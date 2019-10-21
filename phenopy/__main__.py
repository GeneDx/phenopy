import fire
import itertools
import sys

import pandas as pd

from configparser import NoOptionError, NoSectionError
from multiprocessing import Manager, Pool
from phenopy.config import config, logger
from phenopy.network import _load_hpo_network
from phenopy.d2p import load as load_d2p
from phenopy.score import Scorer
from phenopy.util import remove_parents, read_records_file
from phenopy.weights import make_age_distributions


def score(query_hpo_file, records_file=None, query_name='SAMPLE', obo_file=None, disease_to_phenotype_file=None, threads=1,
          agg_score='BMA', no_parents=False, custom_annotations_file=None, output_file=None):
    """
    Scores a case HPO terms against all diseases associated HPO.

    :param query_hpo_file: File with case HPO terms, one per line.
    :param records_file: One record per line, tab delimited. First column record unique identifier, second column
        pipe separated list of HPO identifier (HP:0000001).
    :param query_name: Unique identifier for the query file.
    :param obo_file: OBO file from https://hpo.jax.org/app/download/ontology.
    :param disease_to_phenotype_file: Disease to phenotype annoations from http://compbio.charite.de/jenkins/job/hpo.annotations.2018/
    :param threads: Number of parallel process to use.
    :param agg_score: The aggregation method to use for summarizing the similarity matrix between two term sets
        Must be one of {'BMA', 'maximum'}
    :param no_parents: If provided, scoring is done by only using the most informative nodes. All parent nodes are removed.
    :param custom_annotations_file: A custom entity-to-phenotype annotation file in the same format as tests/data/test.score-product.txt
    :param output_file: filepath where to store the results.
    """

    if agg_score not in {'BMA', 'maximum', 'BMWA'}:
        logger.critical(
            'agg_score must be one of {BMA, maximum, BMWA}.')
        exit(1)

    if obo_file is None:
        try:
            obo_file = config.get('hpo', 'obo_file')
        except (NoSectionError, NoOptionError):
            logger.critical(
                'No HPO OBO file provided and no "hpo:obo_file" found in the configuration file.')
            exit(1)

    if disease_to_phenotype_file is None:
        try:
            disease_to_phenotype_file = config.get('hpo', 'disease_to_phenotype_file')
        except (NoSectionError, NoOptionError):
            logger.critical(
                'No HPO disease_to_phenotype_file file provided and no "hpo:disease_to_phenotype_file" found in the configuration file.'
            )
            exit(1)

    try:
        with open(query_hpo_file, 'r') as case_fh:
            case_hpo = case_fh.read().splitlines()
    except (FileNotFoundError, PermissionError) as e:
        logger.critical(e)
        exit(1)

    # load phenotypes to diseases associations
    disease_to_phenotypes, phenotype_to_diseases, phenotype_disease_frequencies = load_d2p(disease_to_phenotype_file, logger=logger)

    # load hpo network
    hpo_network = _load_hpo_network(
        obo_file, phenotype_to_diseases, len(disease_to_phenotypes), custom_annotations_file,
        phenotype_disease_frequencies=phenotype_disease_frequencies)

    # create instance the scorer class
    scorer = Scorer(hpo_network, agg_score=agg_score)

    # multiprocessing objects
    manager = Manager()
    lock = manager.Lock()

    if no_parents is True:
        case_hpo = remove_parents(case_hpo, hpo_network)

    # convert and filter disease to phenotypes records terms
    for record_id, phenotypes in disease_to_phenotypes.items():
        disease_to_phenotypes[record_id] = scorer.convert_alternate_ids(phenotypes)
        disease_to_phenotypes[record_id] = scorer.filter_and_sort_hpo_ids(phenotypes)

    # convert and filter the query hpo ids
    case_hpo = scorer.convert_alternate_ids(case_hpo)
    case_hpo = scorer.filter_and_sort_hpo_ids(case_hpo)

    if records_file:
        # score and output case hpo terms against all diseases associated set of hpo terms
        logger.info(
            f'Scoring HPO terms from file: {query_hpo_file} against entities in: {records_file}')

        # build the records dictionary
        records_list = read_records_file(records_file, no_parents, hpo_network, logger=logger)
        records = {item['sample']: item['terms'] for item in records_list}
        # clean the records dictionary
        for record_id, phenotypes in records.items():
            records[record_id] = scorer.convert_alternate_ids(phenotypes)
            records[record_id] = scorer.filter_and_sort_hpo_ids(phenotypes)

        records[query_name] = case_hpo
        if not output_file:
            sys.stdout.write('\t'.join(['#query', 'entity_id', 'score']))
            sys.stdout.write('\n')
            with Pool(threads) as p:
                p.starmap(scorer.score_records, [(records, [
                          (query_name, record) for record in records], lock, i, threads, True, False) for i in range(threads)])
        else:
            with Pool(threads) as p:
                scored_results = p.starmap(scorer.score_records, [(records, [(query_name, record) for record in records],
                                                                 lock, i, threads, False, False) for i in range(threads)])
            scored_results = [item for sublist in scored_results for item in sublist]
            scored_results_df = pd.DataFrame(data=scored_results, columns='#query,entity_id,score'.split(','))
            scored_results_df = scored_results_df.sort_values(by='score', ascending=False)
            scored_results_df.to_csv(output_file, sep='\t', index=False)
            logger.info(f'Scoring completed')
            logger.info(f'Writing results to file: {output_file}')

    else:
        # score and output case hpo terms against all disease associated set of hpo terms
        logger.info(f'Scoring case HPO terms from file: {query_hpo_file}')

        # set min_score_mask
        scorer.min_score_mask = None

        # include the case - to - iteslf
        disease_to_phenotypes[query_name] = case_hpo
        # arbitrarily set the query sample as a custom disease and set weights to 1.0 for self-self scoring.
        for hpo_id in case_hpo:
            hpo_network.node[hpo_id]['weights']['disease_frequency'][query_name] = 1.0
        if not output_file:
            sys.stdout.write('\t'.join(['#query', 'omim_id', 'score']))
            sys.stdout.write('\n')
            # iterate over each cross-product and score the pair of records
            with Pool(threads) as p:
                p.starmap(scorer.score_records, [(disease_to_phenotypes, [
                          (query_name, disease_id) for disease_id in disease_to_phenotypes], lock,  i, threads, True, True) for i in range(threads)])
        else:

            with Pool(threads) as p:
                scored_results = p.starmap(scorer.score_records, [(disease_to_phenotypes,
                                     [(query_name, disease_id) for disease_id in disease_to_phenotypes], lock,  i, threads, False, True)
                                                                for i in range(threads)])
            scored_results = [item for sublist in scored_results for item in sublist]
            scored_results_df = pd.DataFrame(data=scored_results, columns='#query,omim_id,score'.split(','))
            scored_results_df = scored_results_df.sort_values(by='score', ascending=False)
            scored_results_df.to_csv(output_file, sep='\t', index=False)
            logger.info(f'Scoring completed')
            logger.info(f'Writing results to file: {output_file}')


def score_product(records_file, obo_file=None, disease_to_phenotype_file=None, pheno_ages_file=None,
                  threads=1, agg_score='BMA', no_parents=False, custom_annotations_file=None):
    """
    Scores the cartesian product of HPO terms from a list of unique records (cases, diseases, diseases, etc).

    :param records_file: One record per line, tab delimited. 1st column record unique identifier, 2nd column contains
        optional patient age and gender(age=11.0;sex=male). 3d column contains pipe separated list of HPO identifier (HP:0000001).
    :param obo_file: OBO file from https://hpo.jax.org/app/download/ontology.
    :param disease_to_phenotype_file: Disease to phenotype annoations from http://compbio.charite.de/jenkins/job/hpo.annotations.2018/
    :param pheno_ages_file: Phenotypes age summary stats file containing phenotype HPO id, mean_age, and std.
    :param threads: Multiprocessing threads to use [default: 1].
    :param agg_score: The aggregation method to use for summarizing the similarity matrix between two term sets
        Must be one of {'BMA', 'maximum', 'BMWA'}. If BMWA is passed phenotype ages file is expected the first time its run.
    :param no_parents: If provided, scoring is done by only using the most informative nodes. All parent nodes are removed.
    :param custom_annotations_file: A custom entity-to-phenotype annotation file in the same format as tests/data/test.score-product.txt
    """
    if agg_score not in {'BMA', 'maximum', 'BMWA'}:
        logger.critical(
            'agg_score must be one of {BMA, maximum, BMWA}.')
        exit(1)

    if obo_file is None:
        try:
            obo_file = config.get('hpo', 'obo_file')
        except (NoSectionError, NoOptionError):
            logger.critical(
                'No HPO OBO file provided and no "hpo:obo_file" found in the configuration file.')
            exit(1)

    if disease_to_phenotype_file is None:
        try:
            disease_to_phenotype_file = config.get('hpo', 'disease_to_phenotype_file')
        except (NoSectionError, NoOptionError):
            logger.critical(
                'No HPO phenotype.hpoa file provided and no "hpo:disease_to_phenotype_file" found in the configuration file.'
            )
            exit(1)
    if pheno_ages_file is not None:
        try:
            ages = make_age_distributions(pheno_ages_file)
            logger.info(
                'Added custom phenotype age distributions to HPO nodes.'
            )
        except (FileNotFoundError, PermissionError) as e:
            logger.critical(e)
            logger.critical(
                'Specified phenotype ages file could not be loaded or does not exist'
            )
            exit(1)
    elif agg_score == 'BMWA':
        try:
            ages = make_age_distributions(config.get('age', 'open_access_phenotype_age'))
            logger.info(
                'Added default phenotype age distributions to HPO nodes.'
            )
        except (FileNotFoundError, PermissionError) as e:
            logger.critical(e)
            logger.critical(
                'Default phenotype ages file could not be loaded or does not exist'
            )
            exit(1)

    else:
        ages = None

    # load phenotype to diseases associations
    disease_to_phenotypes, phenotype_to_diseases, phenotype_disease_frequencies = load_d2p(disease_to_phenotype_file, logger=logger)

    # load hpo network
    hpo_network = _load_hpo_network(
        obo_file, phenotype_to_diseases, len(disease_to_phenotypes), custom_annotations_file, ages=ages, phenotype_disease_frequencies=phenotype_disease_frequencies)

    # try except
    records_list = read_records_file(records_file, no_parents, hpo_network, logger=logger)

    logger.info(f'Scoring product of records from file: {records_file}')

    # create instance the scorer class
    scorer = Scorer(hpo_network, agg_score=agg_score)

    # # clean HPO ids: convert from alternate primary then filter and sort
    # records = {item['sample']: item['terms'] for item in records}
    # # clean the records dictionary
    # for record_id, phenotypes in records.items():
    #     records[record_id] = scorer.convert_alternate_ids(phenotypes)
    #     records[record_id] = scorer.filter_and_sort_hpo_ids(phenotypes)

    # iterate over each cross-product and score the pair of records
    manager = Manager()
    lock = manager.Lock()
    with Pool(threads) as p:
        p.starmap(scorer.score_pairs, [(records_list,
                                        lock, i, threads) for i in range(threads)])


def main():

    fire.Fire({
        'score': score,
        'score-product': score_product,
    })


if __name__ == '__main__':
    main()
