import os

from phenosim.config import config, data_directory, logger
from phenosim.obo import cache, process, restore
from phenosim.obo import load as load_obo


def _load_hpo_network(obo_file, terms_to_genes, annotations_count, custom_annotations_file):
    """
    Load and process phenotypes to genes and obo files if we don't have a processed network already.
    """
    # We instruct the user that they can set hpo_network_file in .phenosim/phenosim.ini
    # The default value is empty string, so check for that first.
    hpo_network_file = config.get('hpo', 'hpo_network_file')
    hpo_network_file_empty = hpo_network_file == ''
    if hpo_network_file_empty:
        # The user didn't set hpo_network_file in the config so try to restore hpo_network.pickle from a previous run.
        found_data_directory = config.get('hpo', 'data_directory', fallback=data_directory)
        # the data_directory variable only exists in tests/test_network.py so the fallback is the default behavior.
        hpo_network_file = os.path.join(found_data_directory, 'hpo_network.pickle')
        # if the default hpo_network.pickle file doesn't exist, process a new hpo_network from scratch
        if not os.path.exists(hpo_network_file):
            # load and process hpo network
            logger.info(f'Loading HPO OBO file: {obo_file}')
            hpo_network = load_obo(obo_file, logger=logger)
            hpo_network = process(hpo_network, terms_to_genes, annotations_count, custom_annotations_file,
                                  logger=logger)

            # save a cache of the processed network
            cache(hpo_network, hpo_network_file)
        # the default hpo_network.pickle file was found
        else:
            hpo_network = restore(hpo_network_file)
    # If the user provides a valid path to the hpo_network file in phenosim.ini, restore that.
    else:
        try:
            hpo_network = restore(hpo_network_file)
        except (FileNotFoundError, PermissionError, IsADirectoryError):
            logger.critical(f'{hpo_network_file} is not a valid path to a pickled hpo_network file.\n'
                            f'In your $HOME/.phenosim/phenosim.ini, please set hpo_network_file'
                            f'=/path/to/hpo_netowrk.pickle OR leave it empty, which is the default. ')
            raise
    return hpo_network