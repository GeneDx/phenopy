import os
from configparser import NoOptionError, NoSectionError

from phenopy.config import config, logger
from phenopy.obo import cache, process, restore
from phenopy.obo import load as load_obo
from phenopy.weights import make_age_distributions


def load(obo_file, phenotype_to_diseases, num_diseases_annotated, annotations_file=None, ages_distribution_file=None,
         phenotype_disease_frequencies=None, hpo_network_file=None):
    """
    Load and process phenotypes to diseases and obo files if we don't have a processed network already.

    :param obo_file: path to obo file.
    :param phenotype_to_diseases: Dictionary of HPO terms as keys and list of diseases as values.
    :param num_diseases_annotated: An integer representing the number of unique diseases in the annotation corpus.
    :param annotations_file: Path to annotations file, optional.
    :param ages_distribution_file: Path to phenotypes ages distribution file.
    :param hpo_network_file: Path to pre-processed HPO network file.
    """
    # We instruct the user that they can set hpo_network_file in .phenopy/phenopy.ini
    # The default value is ~/.phenopy/data/hpo_network.pickle
    if hpo_network_file is None:
        try:
            hpo_network_file = config.get('hpo', 'hpo_network_file')
        except (NoSectionError, NoOptionError):
            logger.critical(
                'No HPO network file found in the configuration file. See "hpo:hpo_network_file" parameter.')
            exit(1)

    # make ages distributions
    ages = None
    if ages_distribution_file is not None:
        try:
            ages = make_age_distributions(ages_distribution_file)
            logger.info(
                f'Adding custom phenotype age distributions to HPO nodes from file: {ages_distribution_file}'
            )
        except (FileNotFoundError, PermissionError) as e:
            logger.critical(e)
            logger.critical(
                f'Specified phenotype ages file could not be loaded or does not exist: {e}'
            )
            exit(1)

    # load and process hpo network if we receive an annotation file or if we don't have a pre-processed one
    if None not in [annotations_file, ages_distribution_file] or not os.path.exists(hpo_network_file):
        logger.info(f'Loading HPO OBO file: {obo_file}')
        hpo_network = load_obo(obo_file, logger=logger)
        hpo_network = process(hpo_network, phenotype_to_diseases, num_diseases_annotated, annotations_file,
                              ages=ages, phenotype_disease_frequencies=phenotype_disease_frequencies, logger=logger)

        # save a cache of the processed network only if no optional files were provided
        if None not in [annotations_file, ages_distribution_file]:
            cache(hpo_network, hpo_network_file)

    # the default hpo_network pickle file was found
    else:
        try:
            hpo_network = restore(hpo_network_file)
        except (FileNotFoundError, PermissionError, IsADirectoryError) as e:
            logger.critical(f'{hpo_network_file} is not a valid path to a pickled hpo_network file.\n'
                            f'In your $HOME/.phenopy/phenopy.ini, please set hpo_network_file'
                            f'=/path/to/hpo_netowrk.pickle OR leave it empty to use the default. ')
            raise e
    return hpo_network
