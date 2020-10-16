import os

from configparser import NoOptionError, NoSectionError

import joblib
# import lightgbm as lgb

from phenopy import generate_annotated_hpo_network
from phenopy.config import config, logger, project_data_dir
from phenopy.util import encode_phenotypes, read_phenotype_groups

def predict_likelihood_moldx(phenotypes, phenotype_groups=None, hpo_network=None, alt2prim=None):
    """
    Predicts the likelihood of molecular diagnosis given a set of phenotypes.
    :param phenotypes: A list of phenotypes or a list of lists of phenotypes.
    :param phenotype_groups: <optionnal> A dictionary of phenotype to phenotype group mappings.
    :param hpo_network: <optional> The hpo networkx object.
    :param alt2prim: <optional> A dictionary of alternate phenotype ids to primary phenotype ids. (must be given if hpo_network is provided)
    :return: An array of probabilities for the positive class.
    """
    # detect if phenotypes is 1d or 2d
    if hpo_network is None or alt2prim is None:
        try:
            obo_file = config.get('hpo', 'obo_file')
        except (NoSectionError, NoOptionError):
            logger.critical(
                'No HPO OBO file found in the configuration file. See "hpo:obo_file" parameter.')
            exit(1)
        try:
            disease_to_phenotype_file = config.get('hpo', 'disease_to_phenotype_file')
        except (NoSectionError, NoOptionError):
            logger.critical(
                'No HPO annotated dataset file found in the configuration file.'
                ' See "hpo:disease_to_phenotype_file" parameter.'
            )
            exit(1)

        logger.info(f'Loading HPO OBO file: {obo_file}')
        hpo_network, alt2prim, _ = \
            generate_annotated_hpo_network(obo_file,
                                        disease_to_phenotype_file,
                                        )
    if phenotype_groups is None:
        phenotype_groups = read_phenotype_groups()

    encoded_phenotypes = encode_phenotypes(phenotypes, phenotype_groups, hpo_network, alt2prim)
    model_file = os.path.join(project_data_dir, 'lgb.model.pkl')
    model = joblib.load(model_file)
    probabilities = model.predict_proba(encoded_phenotypes)
    return probabilities[:, 1]
