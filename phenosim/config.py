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
            with urllib.request.urlopen(url) as response, open(file_path, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)

        except PermissionError:
            logger.info('No permission accessing data directory: %s' % file_path)
            exit(1)

        except ValueError:
            logger.info('Incorrect url specified for HPO files: %s' % url)
            raise

        except urllib.error.URLError:
            logger.info('Incorrect url specified for HPO files: %s' % url)
            raise

        except Exception as unhandled_exception:
            logger.info('Could not download HPO files')
            logger.info('Try downloading them manually')
            raise unhandled_exception


        # read the config file to get file paths and urls
        phen_to_genes_path = config.get('hpo','pheno2genes_file')
        phen_to_genes_url = config.get('hpo','pheno2genes_file_url')

        obo_path = config.get('hpo','obo_file')
        obo_url = config.get('hpo','obo_file_url')

        if not os.path.isfile(phen_to_genes_path):
            logger.info('Downloading HPO phenotype to: %s'%phen_to_genes_path)
            download(phen_to_genes_url, phen_to_genes_path)

        if not os.path.isfile(obo_path):
            logger.info('Downloading HPO obo file to: %s'%obo_path)
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


#set defaults
config.read_dict({
    'hpo': {
        'obo_file': os.path.join(
            data_directory,
            'hp.obo',
        ),
        'obo_file_url':'http://purl.obolibrary.org/obo/hp.obo',

        'pheno2genes_file': os.path.join(
            data_directory,
            'phenotype_to_genes.txt',
        ),
        'pheno2genes_file_url':'http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt'

    },
})

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

