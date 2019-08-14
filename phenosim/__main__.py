import fire
import itertools
import os

from configparser import NoOptionError, NoSectionError
from multiprocessing import Manager, Pool

from phenosim.config import config, data_directory, logger
from phenosim.obo import cache, process, restore
from phenosim.obo import load as load_obo
from phenosim.p2g import load as load_p2g
from phenosim.score import Scorer
from phenosim.util import remove_parents, read_records_file


def _load_hpo_network(obo_file, hpo_network_file, terms_to_genes, annotations_count, custom_annotations_file):
    """
    Load and process phenotypes to genes and obo files if we don't have a processed network already.
    """
    # if the hpo_network_file is provided and it exists, load it, this ignores obo_file and custom_annotations_file
    if hpo_network_file is not None:
        if os.path.exists(hpo_network_file):
            logger.info(f'restoring {hpo_network_file} from disk. Ignoring obo_file and custom_annotations_file')
            hpo_network = restore(hpo_network_file)
        # when the hpo_network_file is provided but does not exist, process the obo file
        else:
            logger.info(f'Loading {obo_file} and writing the cache to: {hpo_network_file}')
            hpo_network = load_obo(obo_file, logger=logger)
            hpo_network = process(hpo_network, terms_to_genes, annotations_count, custom_annotations_file,
                                  logger=logger)
            cache(hpo_network, hpo_network_file)
    else:
        # if the user doesn't provide hpo_network_file, we assume to use the default.
        hpo_network_file = os.path.join(data_directory, 'hpo_network.pickle')
        # restore default hpo_network.pickle
        if os.path.exists(hpo_network_file):
            logger.info(f'Loading the default {hpo_network_file}')
            hpo_network = restore(hpo_network_file)
        # process obo from scratch and write the default pickle.
        else:
            logger.warning(f'processing {obo_file} and writing the default {hpo_network_file}')
            hpo_network = load_obo(obo_file, logger=logger)
            hpo_network = process(hpo_network, terms_to_genes, annotations_count, custom_annotations_file,
                                  logger=logger)
            cache(hpo_network, hpo_network_file)
    return hpo_network


def score(query_hpo_file, records_file=None, query_name='query', obo_file=None, pheno2genes_file=None, threads=1,
          agg_score='BMA', no_parents=False, hpo_network_file=None, custom_annotations_file=None):
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
    :param hpo_network_file: If provided, phenosim will try to load a cached hpo_network obejct from file.
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
        obo_file, hpo_network_file, terms_to_genes, annotations_count, custom_annotations_file)

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
        # iterate over each cross-product and score the pair of records
        with Pool(threads) as p:
            p.starmap(scorer.score_pairs, [(genes_to_terms, [
                      (query_name, gene) for gene in genes_to_terms], lock, agg_score, i, threads) for i in range(threads)])


def score_product(records_file, obo_file=None, pheno2genes_file=None, threads=1, agg_score='BMA', no_parents=False,
                  hpo_network_file=None, custom_annotations_file=None):
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
    :param hpo_network_file: If provided, phenosim will try to load a cached hpo_network obejct from file.
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
        obo_file, hpo_network_file, terms_to_genes, annotations_count, custom_annotations_file)

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
