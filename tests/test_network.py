import os
import pytest

from phenopy.config import config
from phenopy.d2p import load as load_d2p
from phenopy.network import load as load_network
from phenopy.network import annotate
from phenopy.util import generate_alternate_ids


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
    data["obo_file"] = os.path.join(data["parent_dir"], "data/hp.obo")
    return data


def test_load_network(test_data):
    hpo_network = load_network(test_data["obo_file"])
    assert len(hpo_network) == 16861


def test_annotate_network(test_data):
    hpo_network = load_network(test_data["obo_file"])
    alt2prim = generate_alternate_ids(hpo_network)

    # load phenotypes to diseases associations
    disease_to_phenotype_file = os.path.join(
        test_data["parent_dir"], "data/phenotype.hpoa"
    )
    disease_records, phenotype_to_diseases = load_d2p(
        disease_to_phenotype_file, hpo_network, alt2prim
    )

    num_diseases_annotated = len(disease_records)
    hpo_network = annotate(
        hpo_network, phenotype_to_diseases, num_diseases_annotated, alt2prim
    )

    assert round(hpo_network.nodes["HP:0010863"]["ic"], 2) == 7.21
    assert round(hpo_network.nodes["HP:0001263"]["ic"], 2) == 1.55
