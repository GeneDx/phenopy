import csv
import fire
import itertools
import os
import sys

from configparser import NoOptionError, NoSectionError

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


def score_case_to_genes(case_hpo_file, obo_file=None, pheno2genes_file=None):
    """
    Scores a case HPO terms against all genes associated HPO.

    :param case_hpo_file: File with case HPO terms, one per line.
    :param obo_file: OBO file from https://hpo.jax.org/app/download/ontology.
    :param pheno2genes_file: Phenotypes to genes from https://hpo.jax.org/app/download/annotation.
    """
    if obo_file is None:
        try:
            obo_file = config.get('hpo', 'obo_file')
        except (NoSectionError, NoOptionError):
            logger.critical('No HPO OBO file provided and no "hpo:obo_file" found in the configuration file.')
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
        with open(case_hpo_file, 'r') as case_fh:
            case_hpo = case_fh.read().splitlines()
    except (FileNotFoundError, PermissionError) as e:
        logger.critical(e)
        exit(1)

    # load phenotypes to genes associations
    terms_to_genes, genes_to_terms, annotations_count = load_p2g(pheno2genes_file, logger=logger)

    # load hpo network
    hpo_network = _load_hpo_network(obo_file, terms_to_genes, annotations_count)

    # score and output case hpo terms against all genes associated set of hpo terms
    logger.info(f'Scoring case HPO terms from file: {case_hpo_file}')

    # create instance the scorer class
    scorer = Scorer(hpo_network)

    for gene in genes_to_terms.keys():
        gene_hpo = genes_to_terms[gene]

        # skip genes with not hpo terms in the network
        if not gene_hpo:
            continue

        sys.stdout.write('\t'.join([
            gene,
            str(scorer.score(case_hpo, gene_hpo)),
        ]))
        sys.stdout.write('\n')


def score_all(records_file, obo_file=None, pheno2genes_file=None):
    """
    Scores the cross-product of HPO terms from a list of unique records (cases, genes, diseases, etc).

    :param records_file: One record per line, tab delimited. First column record unique identifier, second column
        pipe separated list of HPO identifier (HP:0000001).
    :param obo_file: OBO file from https://hpo.jax.org/app/download/ontology.
    :param pheno2genes_file: Phenotypes to genes from https://hpo.jax.org/app/download/annotation.
    """
    if obo_file is None:
        try:
            obo_file = config.get('hpo', 'obo_file')
        except (NoSectionError, NoOptionError):
            logger.critical('No HPO OBO file provided and no "hpo:obo_file" found in the configuration file.')
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
        # read records_file and convert to DataFrame #record_id | hpo_ids
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
    terms_to_genes, _, annotations_count = load_p2g(pheno2genes_file, logger=logger)

    # load hpo network
    hpo_network = _load_hpo_network(obo_file, terms_to_genes, annotations_count)

    # create instance the scorer class
    scorer = Scorer(hpo_network)

    # iterate over each cross-product and score the pair of records
    for record_a, record_b in itertools.product(records, repeat=2):
        sys.stdout.write('\t'.join([
            f'{record_a}-{record_b}',
            str(scorer.score(records[record_a], records[record_b])),
        ]))
        sys.stdout.write('\n')


def main():
    fire.Fire({
        'score': score_case_to_genes,
        'score-all': score_all,
    })


if __name__ == '__main__':
    main()
