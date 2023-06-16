import os
import pytest

from phenopy.network import annotate
from phenopy.network import load as load_network
from phenopy.util import export_phenotype_hpoa_with_no_parents


def test_ic_d2p(test_data):
    """Calculate the information content of a phenotype"""
    assert (
        pytest.approx(test_data["hpo_network"].nodes["HP:0010863"]["ic"], 0.01) == 7.21
    )


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
        hpo_network,
        test_data["phenotype_to_diseases"],
        test_data["num_diseases_annotated"],
        test_data["alt2prim"],
        annotations_file=custom_annotation_file,
    )

    assert pytest.approx(hpo_network.nodes["HP:0010863"]["ic"], 0.01) == 8.11


def test_ic_d2p_no_parents(test_data):
    export_phenotype_hpoa_with_no_parents(
        test_data["disease_to_phenotype_file"],
        test_data["disease_to_phenotype_output_file"],
        test_data["hpo_network"],
    )
    assert os.path.exists(test_data["disease_to_phenotype_output_file"])
    os.remove(test_data["disease_to_phenotype_output_file"])
