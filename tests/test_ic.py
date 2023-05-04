import os
import pytest

from phenopy.network import annotate
from phenopy.network import load as load_network
from phenopy.d2p import load as load_d2p
from phenopy.util import export_phenotype_hpoa_with_no_parents, generate_alternate_ids


@pytest.fixture(scope="module")
def test_data():
    data = {}
    data["parent_dir"] = os.path.dirname(os.path.realpath(__file__))
    data["obo_file"] = os.path.join(data["parent_dir"], "data/hp.obo")
    data["hpo_network"] = load_network(data["obo_file"])
    data["alt2prim"] = generate_alternate_ids(data["hpo_network"])
    data["disease_to_phenotype_file"] = os.path.join(
        data["parent_dir"], "data/phenotype.hpoa"
    )
    data["disease_records"], data["phenotype_to_diseases"] = load_d2p(
        data["disease_to_phenotype_file"], data["hpo_network"], data["alt2prim"]
    )
    data["num_diseases_annotated"] = len(data["disease_records"])
    data["hpo_network"] = annotate(
        data["hpo_network"], data["phenotype_to_diseases"],
        data["num_diseases_annotated"], data["alt2prim"]
    )

    data["hpo_id"] = "HP:0010863"
    data["disease_to_phenotype_output_file"] = os.path.join(
        data["parent_dir"], "data/phenotype.noparents.hpoa"
    )
    return data


def test_ic_d2p(test_data):
    """Calculate the information content of a phenotype"""
    assert round(test_data["hpo_network"].nodes[test_data["hpo_id"]]["ic"], 2) == 7.21


def test_ic_custom(test_data):
    """
    Calculate the information content of a phenotype when multiple
    annotations are present
    """
    custom_annotation_file = os.path.join(
        test_data["parent_dir"], "data/test.score-long.txt"
    )
    hpo_network = load_network(test_data["obo_file"])
    hpo_network = annotate(
        hpo_network, test_data["phenotype_to_diseases"],
        test_data["num_diseases_annotated"],
        test_data["alt2prim"], annotations_file=custom_annotation_file)

    assert round(hpo_network.nodes[test_data["hpo_id"]]["ic"], 2) == 8.11


def test_ic_d2p_no_parents(test_data):
    export_phenotype_hpoa_with_no_parents(
        test_data["disease_to_phenotype_file"],
        test_data["disease_to_phenotype_output_file"],
        test_data["hpo_network"])
    assert os.path.exists(test_data["disease_to_phenotype_output_file"])
    os.remove(test_data["disease_to_phenotype_output_file"])

