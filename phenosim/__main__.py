import csv
import fire
import itertools
import os
import sys

import pandas as pd

from configparser import NoOptionError, NoSectionError
from multiprocessing import Manager, Pool

from phenosim.cluster import clustering_grid_search, clustering_assign
from phenosim.config import config, data_directory, logger
from phenosim.obo import cache, process, restore
from phenosim.obo import load as load_obo
from phenosim.p2g import load as load_p2g
from phenosim.score import Scorer


def _load_hpo_network(obo_file, terms_to_genes, annotations_count):
    """
    Load and process phenotypes to genes and obo files if we don't have a processed network already.
    """
    hpo_network_file = os.path.join(data_directory, 'hpo_network.pickle')
    if not os.path.exists(hpo_network_file):
        # load and process hpo network
        logger.info(f'Loading HPO OBO file: {obo_file}')
        hpo_network = load_obo(obo_file, logger=logger)
        hpo_network = process(hpo_network, terms_to_genes, annotations_count)

        # save a cache of the processed network
        cache(hpo_network, hpo_network_file)
    else:
        hpo_network = restore(hpo_network_file)

    return hpo_network


def score(query_hpo_file, records_file=None, query_name='query', obo_file=None, pheno2genes_file=None, threads=1,):
    """
    Scores a case HPO terms against all genes associated HPO.

    :param query_hpo_file: File with case HPO terms, one per line.
    :param records_file: One record per line, tab delimited. First column record unique identifier, second column
        pipe separated list of HPO identifier (HP:0000001).
    :param query_name: Unique identifier for the query file.
    :param obo_file: OBO file from https://hpo.jax.org/app/download/ontology.
    :param pheno2genes_file: Phenotypes to genes from https://hpo.jax.org/app/download/annotation.
    :param threads: Number of parallel process to use.

    """
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
        obo_file, terms_to_genes, annotations_count)

    # create instance the scorer class
    scorer = Scorer(hpo_network)

    # multiprocessing objects
    manager = Manager()
    lock = manager.Lock()

    if records_file:
        # score and output case hpo terms against all genes associated set of hpo terms
        logger.info(
            f'Scoring case HPO terms from file: {query_hpo_file} against cases in: {records_file}')
        try:
            # read records_file
            with open(records_file) as records_fh:
                reader = csv.reader(records_fh, delimiter='\t')
                records = {}
                for line in reader:
                    if line[0].startswith('#'):
                        continue
                    records[line[0]] = line[1].split('|')
        except (FileNotFoundError, PermissionError) as e:
            logger.critical(e)
            exit(1)

        # include the case-to-iteslf
        records[query_name] = case_hpo
        with Pool(threads) as p:
            p.starmap(scorer.score_pairs, [(records, [
                      (query_name, record) for record in records], lock, i, threads) for i in range(threads)])

    else:
        # score and output case hpo terms against all genes associated set of hpo terms
        logger.info(f'Scoring case HPO terms from file: {query_hpo_file}')

        # add the case terms to the genes_to_terms dict
        genes_to_terms[query_name] = case_hpo
        # iterate over each cross-product and score the pair of records
        with Pool(threads) as p:
            p.starmap(scorer.score_pairs, [(genes_to_terms, [
                      (query_name, gene) for gene in genes_to_terms], lock, i, threads) for i in range(threads)])


def score_product(records_file, obo_file=None, pheno2genes_file=None, threads=1):
    """
    Scores the cartesian product of HPO terms from a list of unique records (cases, genes, diseases, etc).

    :param records_file: One record per line, tab delimited. First column record unique identifier, second column
        pipe separated list of HPO identifier (HP:0000001).
    :param obo_file: OBO file from https://hpo.jax.org/app/download/ontology.
    :param pheno2genes_file: Phenotypes to genes from https://hpo.jax.org/app/download/annotation.
    :param threads: Multiprocessing threads to use [default: 1].
    """
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
        # read records_file
        with open(records_file) as records_fh:
            reader = csv.reader(records_fh, delimiter='\t')
            records = {}
            for line in reader:
                if line[0].startswith('#'):
                    continue
                records[line[0]] = line[1].split('|')
    except (FileNotFoundError, PermissionError) as e:
        logger.critical(e)
        exit(1)

    # load phenotypes to genes associations
    terms_to_genes, _, annotations_count = load_p2g(
        pheno2genes_file, logger=logger)

    # load hpo network
    hpo_network = _load_hpo_network(
        obo_file, terms_to_genes, annotations_count)

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
                                        lock, i, threads) for i in range(threads)])


def cluster_grid_search(score_product_result_file, max_clusters=2):
    """Runs clustering algorithms in parallel on the output of phenosim `score-product`
    :param score_product_result_file: path to file
    :type score_product_result_file: str
    :param max_clusters: The maximum number of clusters to output silhouette score for.
    :type max_clusters: int
    """

    if max_clusters < 2:
        logger.critical('max_clusters must be >=2')
        exit(1)

    LINKAGE = {"complete", "average", "single"}
    linkage_k_combos = itertools.product(LINKAGE, range(2, max_clusters + 1))

    try:
        # read phenosim_result_file
        df = pd.read_csv(score_product_result_file,
                         sep='\t',
                         header=None,
                         names=['id_pairs', 'score'])
    except (FileNotFoundError, PermissionError) as e:
        logger.critical(e)
        exit(1)

    # process the DataFrame
    df['id_pairs'] = df['id_pairs'].str.split('-')
    df[['record1', 'record2']] = pd.DataFrame(
        df['id_pairs'].values.tolist(), index=df.index)
    df.drop('id_pairs', axis=1, inplace=True)
    df = df.set_index(['record1', 'record2']).unstack()
    X = df.values

    for link_method, k in linkage_k_combos:
        clustering_grid_search(X, link_method, k)


def assign_clusters(score_product_result_file, linkage='average', k=2):
    """Runs agglomerative clustering algorithms in parallel on the output of phenosim `score-product`
    :param score_product_result_file: path to file
    :type score_product_result_file: str
    :param linkage: The type of linkage to perform {single, average, complete}
    :type linkage: str
    :param k:
    :type max_clusters: int
    :param *kwargs: Placeholder for arguements to pass to pyclustering or scikit-learn
    """
    if not any(_ == linkage for _ in ['single', 'average', 'complete']):
        sys.stderr('Please select one of the allowed clustering methods')
        exit(1)

    if k <= 1:
        sys.stderr('Please select k to be an int >= 2')
        exit(1)

    try:
        # read phenosim_result_file
        df = pd.read_csv(score_product_result_file,
                         sep='\t',
                         header=None,
                         names=['id_pairs', 'score'])
    except (FileNotFoundError, PermissionError) as e:
        logger.critical(e)
        exit(1)

    # process the DataFrame
    df['id_pairs'] = df['id_pairs'].str.split('-')
    df[['record1', 'record2']] = pd.DataFrame(
        df['id_pairs'].values.tolist(), index=df.index)
    df.drop('id_pairs', axis=1, inplace=True)
    df = df.set_index(['record1', 'record2']).unstack()
    X = df.values
    samples = df.index.tolist()

    clustering_assign(X, linkage, k, samples)


def main():
    fire.Fire({
        'score': score,
        'score-product': score_product,
        'cluster-grid-search': cluster_grid_search,
        'cluster-assign': assign_clusters,
    })


if __name__ == '__main__':
    main()
