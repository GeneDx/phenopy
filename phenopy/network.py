import os

from phenopy.config import config, logger
from phenopy.obo import cache, process, restore
from phenopy.obo import load as load_obo


def _load_hpo_network(obo_file, terms_to_genes, num_genes_annotated, custom_annotations_file, hpo_network_file=None):
    """
    :param obo_file: path to obo file.
    :param terms_to_genes: Dictionary of HPO terms as keys and list of genes as values.
    :param num_genes_annotated: An integer representing the number of unique genes in the annotation corpus.
    :param custom_annotations_file: Path to a custom annotations file, if it exists.
    :param hpo_network_file: Path to the hpop_network_file.
    Load and process phenotypes to genes and obo files if we don't have a processed network already.
    """
    # We instruct the user that they can set hpo_network_file in .phenopy/phenopy.ini
    # The default value is empty string, so check for that first.

    if not isinstance(terms_to_genes, dict):
        logger.critical(f'terms_to_genes was not a dictionary, please use the terms_to_genes variable returned from '
                        f'load_p2g')
        raise ValueError

    if hpo_network_file is None:
        hpo_network_file = config.get('hpo', 'hpo_network_file')

    if not os.path.exists(hpo_network_file):
        # load and process hpo network
        logger.info(f'Loading HPO OBO file: {obo_file}')
        hpo_network = load_obo(obo_file, logger=logger)
        hpo_network = process(hpo_network, terms_to_genes, num_genes_annotated, custom_annotations_file,
                              logger=logger)

        # save a cache of the processed network
        cache(hpo_network, hpo_network_file)
    # the default hpo_network.pickle file was found
    else:
        try:
            hpo_network = restore(hpo_network_file)
        except (FileNotFoundError, PermissionError, IsADirectoryError) as e:
            logger.critical(f'{hpo_network_file} is not a valid path to a pickled hpo_network file.\n'
                            f'In your $HOME/.phenopy/phenopy.ini, please set hpo_network_file'
                            f'=/path/to/hpo_netowrk.pickle OR leave it empty, which is the default. ')
            raise e
    return hpo_network