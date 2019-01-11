import fire
import os
import sys

from configparser import NoOptionError, NoSectionError

from phenosim.config import config, data_directory, logger
from phenosim.obo import cache, process, restore
from phenosim.obo import load as load_obo
from phenosim.p2g import load as load_p2g
from phenosim.score import score as score_phenotypes
from phenosim.scoreterms import score_case_genes


def score(case_hpo_file, obo_file=None, pheno2genes_file=None):
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
        if logger is not None:
            logger.critical(e)
        else:
            sys.stderr.write(str(e))
        exit(1)

    # load phenotypes to genes associations
    terms_to_genes, genes_to_terms, annotations_count = load_p2g(pheno2genes_file, logger=logger)

    # load/process phenotypes to genes and obo files if we don't have a processed network already
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

    # score and output case hpo terms against all genes associated set of hpo terms
    logger.info(f'Scoring case HPO terms from file: {case_hpo_file}')
    for gene, gene_hpo in genes_to_terms.items():
        # filter out gene hpo terms not in the network
        gene_hpo = list(filter(lambda x: x in hpo_network.node, gene_hpo))

        # skip genes with not hpo terms in the network
        if not gene_hpo:
            continue

        sys.stdout.write('\t'.join([
            gene,
            str(score_phenotypes(case_hpo, gene_hpo, hpo_network)),
        ]))
        sys.stdout.write('\n')

    # # score case hpo terms against all genes associated set of hpo terms
    # logger.info(f'Scoring case HPO terms from file: {case_hpo_file}')
    # results = score_case_genes(hpo_network, case_hpo, genes_to_terms)
    # # write results to file or something, it's currently a dictionary


def main():
    fire.Fire({
        'score': score,
    })


if __name__ == '__main__':
    main()
