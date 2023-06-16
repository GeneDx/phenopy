import os
import pytest

from phenopy.d2p import load as load_d2p
from phenopy.network import load as load_network
from phenopy.network import annotate
from phenopy.util import generate_alternate_ids


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

    assert pytest.approx(hpo_network.nodes["HP:0010863"]["ic"], 0.01) == 7.21
    assert pytest.approx(hpo_network.nodes["HP:0001263"]["ic"], 0.01) == 1.55
