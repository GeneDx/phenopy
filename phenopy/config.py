import configparser
import logging
import os
import urllib.request
import shutil
from pathlib import Path
from gensim.models import KeyedVectors
from phenopy import __project__, __version__


def download_resource_files():
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
            logger.info(f"Incorrect url specified for HPO files: {url}")
            raise

        except urllib.error.URLError as e:
            if hasattr(e, "reason"):
                logger.info(f"Incorrect url specified for HPO files: {url}")
                logger.info("Reason: ", e.reason)
                raise
            elif hasattr(e, "code"):
                logger.info("The server could not fulfill the request")
                logger.info("Reason: ", e.code)
                raise

        try:
            with open(file_path, "wb") as out_file:
                shutil.copyfileobj(response, out_file)

        except PermissionError:
            logger.info(f"No permission accessing data directory: {file_path}")
            raise

    # read the config file to get file paths and urls
    obo_path = config.get("hpo", "obo_file")
    obo_url = config.get("hpo", "obo_file_url")

    hpoa_path = config.get("hpo", "disease_to_phenotype_file")
    hpoa_url = config.get("hpo", "disease_to_phenotype_file_url")

    if not os.path.isfile(obo_path):
        logger.info(f"Downloading HPO obo file to: {obo_path}")
        download(obo_url, obo_path)

    if not os.path.isfile(hpoa_path):
        logger.info(f"Downloading phenotype to disease annotations to {hpoa_path}")
        download(hpoa_url, hpoa_path)


# create logger
logger = logging.getLogger(__project__)
logger.setLevel(logging.DEBUG)

# create console handler
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter and add it to the handler
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
ch.setFormatter(formatter)

# add the handler to the logger
logger.addHandler(ch)

# create config
config = configparser.ConfigParser()

# create config directory if it doesn't exist
config_directory = os.path.join(os.environ.get("HOME"), f".{__project__}")
project_directory = Path(__file__).parent
project_data_dir = os.path.join(project_directory, "data")
try:
    os.makedirs(config_directory)
except FileExistsError:
    pass

# create data directory if it doesn't exist
data_directory = os.path.join(config_directory, "data")
try:
    os.makedirs(data_directory)
except FileExistsError:
    pass

# if phenopy.ini doesnt exist make one
logger.info(f"checking if config file exists: {config_directory}")
if not os.path.isfile(os.path.join(config_directory, "phenopy.ini")):
    config = configparser.ConfigParser()
    w2v_path = os.path.join(os.path.dirname(__file__), "data/phenopy.wv.model.txt.gz")

    w2v_vw_path = os.path.join(data_directory, "phenopy.w2v.model")

    wv = KeyedVectors.load_word2vec_format(w2v_path)
    # save model in faster to load format in users directory
    wv.save(w2v_vw_path)

    # copy the lmd model to the data directory
    lmd_path = os.path.join(os.path.dirname(__file__), "data/lgb.model.pkl")
    lmd_data_path = os.path.join(data_directory, "lgb.model.pkl")
    shutil.copyfile(lmd_path, lmd_data_path)

    config["hpo"] = {
        "obo_file": os.path.join(
            data_directory,
            "hp.obo",
        ),
        "obo_file_url": "http://purl.obolibrary.org/obo/hp.obo",
        "hpo_network_file": os.path.join(
            data_directory,
            "hpo_network.pickle",
        ),
        "disease_to_phenotype_file_url": "http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa",
        "disease_to_phenotype_file": os.path.join(
            data_directory,
            "phenotype.hpoa",
        ),
    }

    config["models"] = {
        "phenopy.wv.model": w2v_vw_path,
    }
    config["age"] = {
        "open_access_phenotype_age": os.path.join(
            project_data_dir,
            "oa_phenotype_age.tsv",
        )
    }
    config["omim"] = {
        "omim_api_key": "",
    }
    config["phenotype_groups"] = {
        "phenotype_groups_file": os.path.join(project_data_dir, "phenotype_groups.txt")
    }

    with open(os.path.join(config_directory, "phenopy.ini"), "w") as configfile:
        logger.info("writing config file to: %s " % config_directory)
        config.write(configfile)

# log project and version
logger.info(f"{__project__} {__version__}")

# read config
config_file = os.environ.get(
    f"{__project__.upper()}_CONFIG",
    os.path.join(
        config_directory,
        f"{__project__}.ini",
    ),
)
config.read(config_file)
logger.info(f"Using configuration file: {config_file}")

# download resource files if necessary
download_resource_files()
