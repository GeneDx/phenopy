import configparser
import logging
import os
import urllib.request
import shutil

from phenosim import __project__, __version__


def download_hpo_files():
    """
    Check if HPO files exist, if not download them
    :return: None
    """
    def download(url, file_path):
        """
        Download and save a file
        :param url: where to get it from
        :param file_path: where to put it
        :return: None
        """
        try:
            response = urllib.request.urlopen(url)

        except ValueError:
            logger.info(f'Incorrect url specified for HPO files: {url}')
            raise

        except urllib.error.URLError as e:
            if hasattr(e, 'reason'):
                logger.info(f'Incorrect url specified for HPO files: {url}')
                raise
                logger.info('Reason: ', e.reason)
            elif hasattr(e, 'code'):
                logger.info('The server could not fulfill the request')
                logger.info('Reason: ', e.code)
                raise


        try:
            with open(file_path, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)

        except PermissionError:
            logger.info(f'No permission accessing data directory: {file_path}')
            raise

    # read the config file to get file paths and urls
    phen_to_genes_path = config.get('hpo', 'pheno2genes_file')
    phen_to_genes_url = config.get('hpo', 'pheno2genes_file_url')

    obo_path = config.get('hpo', 'obo_file')
    obo_url = config.get('hpo', 'obo_file_url')

    if not os.path.isfile(phen_to_genes_path):
        logger.info(f'Downloading HPO phenotype to: {phen_to_genes_path}')
        download(phen_to_genes_url, phen_to_genes_path)

    if not os.path.isfile(obo_path):
        logger.info(f'Downloading HPO obo file to: {obo_path}')
        download(obo_url, obo_path)


# create logger
logger = logging.getLogger(__project__)
logger.setLevel(logging.DEBUG)

# create console handler
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter and add it to the handler
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)

# add the handler to the logger
logger.addHandler(ch)

# create config
config = configparser.ConfigParser()

# create config directory if it doesn't exist
config_directory = os.path.join(os.environ.get('HOME'), f'.{__project__}')
try:
    os.makedirs(config_directory)
except FileExistsError:
    pass

# create data directory if it doesn't exist
data_directory = os.path.join(config_directory, 'data')
try:
    os.makedirs(data_directory)
except FileExistsError:
    pass

# if phenosim.ini doesnt exist make one
logger.info(f'checking if config file exists: {config_directory}')
if not os.path.isfile(os.path.join(config_directory, 'phenosim.ini')):
    config = configparser.ConfigParser()
    config['hpo'] = {
            'obo_file': os.path.join(
                data_directory,
                'hp.obo',
            ),
            'obo_file_url':'http://purl.obolibrary.org/obo/hp.obo',
            'pheno2genes_file': os.path.join(
                data_directory,
                'phenotype_to_genes.txt',
            ),
            'pheno2genes_file_url': 'http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt',
            'hpo_network_file': os.path.join(
                data_directory,
                'hpo_network.pickle',
            ),
        }

    with open(os.path.join(config_directory, 'phenosim.ini'), 'w') as configfile:
        logger.info('writing config file to: %s '%config_directory)
        config.write(configfile)

# log project and version
logger.info(f'{__project__} {__version__}')

# read config
config_file = os.environ.get(
    f'{__project__.upper()}_CONFIG',
    os.path.join(config_directory, f'{__project__}.ini',)
)
config.read(config_file)
logger.info(f'Using configuration file: {config_file}')

logger.info('Checking if HPO files exist')
download_hpo_files()

