import fire
import itertools
import sys
from configparser import NoOptionError, NoSectionError
from multiprocessing import Manager, Pool
import pandas as pd
from phenosim.config import config, logger
from phenosim.network import _load_hpo_network
from phenosim.p2g import load as load_p2g
from phenosim.score import Scorer
from phenosim.util import remove_parents, read_records_file


def score(query_hpo_file, records_file=None, query_name='SAMPLE', obo_file=None, pheno2genes_file=None, threads=1,
          agg_score='BMA', no_parents=False, custom_annotations_file=None, output_file=None):
    """
    Scores a case HPO terms against all genes associated HPO.

    :param query_hpo_file: File with case HPO terms, one per line.
    :param records_file: One record per line, tab delimited. First column record unique identifier, second column
        pipe separated list of HPO identifier (HP:0000001).
    :param query_name: Unique identifier for the query file.
    :param obo_file: OBO file from https://hpo.jax.org/app/download/ontology.
    :param pheno2genes_file: Phenotypes to genes from https://hpo.jax.org/app/download/annotation.
    :param threads: Number of parallel process to use.
    :param agg_score: The aggregation method to use for summarizing the similarity matrix between two term sets
        Must be one of {'BMA', 'maximum'}
    :param no_parents: If provided, scoring is done by only using the most informative nodes. All parent nodes are removed.
    :param custom_annotations_file: A comma-separated list of custom annotation files in the same format as tests/data/test.score-product.txt
    """

    if agg_score not in {'BMA', 'maximum', }:
        logger.critical(
            'agg_score must be one of {BMA, maximum}.')
        exit(1)

    if obo_file is None:
        try:
            obo_file = config.get('hpo', 'obo_file')
        except (NoSectionError, NoOptionError):
            logger.critical(
                'No HPO OBO file provided and no "hpo:obo_file" found in the configuration file.')
            exit(1)

    if pheno2genes_file is None:
        try:
            pheno2genes_file = config.get('hpo', 'pheno2genes_file')
        except (NoSectionError, NoOptionError):
            logger.critical(
                'No HPO pheno2genes_file file provided and no "hpo:pheno2genes_file" found in the configuration file.'
            )
            exit(1)

    try:
        with open(query_hpo_file, 'r') as case_fh:
            case_hpo = case_fh.read().splitlines()
    except (FileNotFoundError, PermissionError) as e:
        logger.critical(e)
        exit(1)

    # load phenotypes to genes associations
    terms_to_genes, genes_to_terms, annotations_count = load_p2g(
        pheno2genes_file, logger=logger)

    # load hpo network
    hpo_network = _load_hpo_network(
        obo_file, terms_to_genes, annotations_count, custom_annotations_file)

    # create instance the scorer class
    scorer = Scorer(hpo_network)

    # multiprocessing objects
    manager = Manager()
    lock = manager.Lock()

    if no_parents is True:
        case_hpo = remove_parents(case_hpo, hpo_network)

    if records_file:
        # score and output case hpo terms against all genes associated set of hpo terms
        logger.info(
            f'Scoring case HPO terms from file: {query_hpo_file} against cases in: {records_file}')

        records = read_records_file(records_file, no_parents, hpo_network, logger=logger)

        # include the case-to-iteslf
        records[query_name] = case_hpo
        with Pool(threads) as p:
            p.starmap(scorer.score_pairs, [(records, [
                      (query_name, record) for record in records], lock, agg_score, i, threads) for i in range(threads)])

    else:
        # score and output case hpo terms against all genes associated set of hpo terms
        logger.info(f'Scoring case HPO terms from file: {query_hpo_file}')

        # add the case terms to the genes_to_terms dict
        genes_to_terms[query_name] = case_hpo
        if not output_file:
            sys.stdout.write('\t'.join(['#Query', 'Gene', 'Score']))
            sys.stdout.write('\n')
            # iterate over each cross-product and score the pair of records
            with Pool(threads) as p:
                p.starmap(scorer.score_pairs, [(genes_to_terms, [
                          (query_name, gene) for gene in genes_to_terms], lock, agg_score, i, threads) for i in range(threads)])
        else:

            with Pool(threads) as p:
                scored_results = p.starmap(scorer.score_pairs, [(genes_to_terms,
                                     [(query_name, gene) for gene in genes_to_terms], lock, agg_score, i, threads, False)
                                                                for i in range(threads)])
            scored_results_df = pd.DataFrame(data=scored_results, columns='query,gene,score'.split(','))
            scored_results_df = scored_results_df.sort_values(by='score')
            print(scored_results_df.head())
            scored_results_df.to_csv(output_file, sep='\t')
            logger.info(f'Scoring completed')
            logger.info(f'Writing results to file: {output_file}')


def score_product(records_file, obo_file=None, pheno2genes_file=None, threads=1, agg_score='BMA', no_parents=False,
                  custom_annotations_file=None):
    """
    Scores the cartesian product of HPO terms from a list of unique records (cases, genes, diseases, etc).

    :param records_file: One record per line, tab delimited. First column record unique identifier, second column
        pipe separated list of HPO identifier (HP:0000001).
    :param obo_file: OBO file from https://hpo.jax.org/app/download/ontology.
    :param pheno2genes_file: Phenotypes to genes from https://hpo.jax.org/app/download/annotation.
    :param threads: Multiprocessing threads to use [default: 1].
    :param agg_score: The aggregation method to use for summarizing the similarity matrix between two term sets
        Must be one of {'BMA', 'maximum'}
    :param no_parents: If provided, scoring is done by only using the most informative nodes. All parent nodes are removed.
    :param custom_annotations_file: A comma-separated list of custom annotation files in the same format as tests/data/test.score-product.txt
    """
    if agg_score not in {'BMA', 'maximum', }:
        logger.critical(
            'agg_score must be one of {BMA, maximum}.')
        exit(1)

    if obo_file is None:
        try:
            obo_file = config.get('hpo', 'obo_file')
        except (NoSectionError, NoOptionError):
            logger.critical(
                'No HPO OBO file provided and no "hpo:obo_file" found in the configuration file.')
            exit(1)

    if pheno2genes_file is None:
        try:
            pheno2genes_file = config.get('hpo', 'pheno2genes_file')
        except (NoSectionError, NoOptionError):
            logger.critical(
                'No HPO pheno2genes_file file provided and no "hpo:pheno2genes_file" found in the configuration file.'
            )
            exit(1)

    # load phenotypes to genes associations
    terms_to_genes, _, annotations_count = load_p2g(
        pheno2genes_file, logger=logger)

    # load hpo network
    hpo_network = _load_hpo_network(
        obo_file, terms_to_genes, annotations_count, custom_annotations_file)

    # try except
    records = read_records_file(records_file, no_parents, hpo_network, logger=logger)

    logger.info(f'Scoring product of records from file: {records_file}')

    # create instance the scorer class
    scorer = Scorer(hpo_network)

    # create records product generator
    records_product = itertools.product(records.keys(), repeat=2)

    # iterate over each cross-product and score the pair of records
    manager = Manager()
    lock = manager.Lock()
    with Pool(threads) as p:
        p.starmap(scorer.score_pairs, [(records, records_product,
                                        lock, agg_score, i, threads) for i in range(threads)])


def main():

    fire.Fire({
        'score': score,
        'score-product': score_product,
    })


if __name__ == '__main__':
    main()
