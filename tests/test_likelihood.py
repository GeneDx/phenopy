import os
import pytest
from configparser import NoOptionError, NoSectionError

from phenopy.config import config
from phenopy.util import read_phenotype_groups
from phenopy.build_hpo import generate_annotated_hpo_network
from phenopy.likelihood import predict_likelihood_moldx


@pytest.fixture(scope="module")
def test_data():
    data = {}
    data["parent_dir"] = os.path.dirname(os.path.realpath(__file__))

    if "hpo" not in config.sections():
        config.add_section("hpo")

    config.set("hpo", "obo_file", os.path.join(data["parent_dir"], "data/hp.obo"))
    config.set("hpo", "disease_to_phenotype_file", os.path.join(
        data["parent_dir"], "data/phenotype.hpoa")
               )

    data["obo_file"] = config.get("hpo", "obo_file")
    data["disease_to_phenotype_file"] = config.get("hpo", "disease_to_phenotype_file")

    data["hpo_network"], data["alt2prim"], data["disease_records"] = \
        generate_annotated_hpo_network(
            data["obo_file"],
            data["disease_to_phenotype_file"],
        )
    data["phenotype_groups"] = read_phenotype_groups()
    return data


def test_read_phenotype_groups_1000(test_data):
    hp_to_pg = read_phenotype_groups()
    assert hp_to_pg["HP:0012759"]["k1000"] == 62


def test_read_phenotype_groups_1500(test_data):
    hp_to_pg = read_phenotype_groups()
    assert hp_to_pg["HP:0012759"]["k1500"] == 739


def test_predict_likelihood(test_data):
    phenotypes = [
        ["HP:0012759", "HP:0003011", "HP:0011442"],
        ["HP:0012759", "HP:0003011"],
    ]
    probabilities = predict_likelihood_moldx(
        phenotypes, test_data["phenotype_groups"],
        test_data["hpo_network"], test_data["alt2prim"]
    )
    assert round(probabilities[0], 2) == 0.33


def test_predict_likelihood_phenotypes_only(test_data):
    phenotypes = [
        ["HP:0012759", "HP:0003011", "HP:0011442"],
        ["HP:0012759", "HP:0003011"],
    ]
    probabilities = predict_likelihood_moldx(phenotypes)
    assert round(probabilities[0], 2) == 0.33


def test_no_hpo_config_section(test_data):
    config.remove_section("hpo")
    phenotypes = [
        ["HP:0012759", "HP:0003011", "HP:0011442"],
        ["HP:0012759", "HP:0003011"],
    ]
    with pytest.raises(NoSectionError):
        predict_likelihood_moldx(phenotypes)


def test_no_hpo_config_option(test_data):
    config.remove_option("hpo", "disease_to_phenotype_file")
    phenotypes = [
        ["HP:0012759", "HP:0003011", "HP:0011442"],
        ["HP:0012759", "HP:0003011"],
    ]
    with pytest.raises(NoOptionError):
        predict_likelihood_moldx(phenotypes)


def test_bad_k(test_data):
    phenotypes = [
        ["HP:0012759", "HP:0003011", "HP:0011442"],
        ["HP:0012759", "HP:0003011"],
    ]
    with pytest.raises(KeyError):
        predict_likelihood_moldx(
            phenotypes,
            test_data["phenotype_groups"],
            test_data["hpo_network"],
            test_data["alt2prim"],
            k_phenotype_groups=500
        )
