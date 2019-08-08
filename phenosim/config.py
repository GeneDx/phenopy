import configparser
import logging
import os

from phenosim import __project__, __version__

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


# set defaults
config.read_dict({
    'hpo': {
        'obo_file': os.path.join(
            data_directory,
            'hp.obo',
        ),
        'pheno2genes_file': os.path.join(
            data_directory,
            'phenotype_to_genes.txt',
        ),
        'hpo_pickle': os.path.join(
            data_directory,
            'hpo_network.pickle',
        ),
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
