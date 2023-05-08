import os
import pytest

from phenopy.d2p import load as load_d2p
from phenopy.network import annotate
from phenopy.network import load as load_network
from phenopy.score import Scorer
from phenopy.util import generate_alternate_ids, read_phenotype_groups


@pytest.fixture(scope="function")
def test_data():
    data = {}
    data["parent_dir"] = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    data["obo_file"] = os.path.join(data["parent_dir"], "data/hp.obo")
    data["hpo_network"] = load_network(data["obo_file"])
    data["alt2prim"] = generate_alternate_ids(data["hpo_network"])
    data["ages_distribution_file"] = os.path.join(
        data["parent_dir"], "data/phenotype_age.tsv"
    )

    data["disease_to_phenotype_file"] = os.path.join(
        data["parent_dir"], "data/phenotype.hpoa"
    )
    data["disease_records"], data["phenotype_to_diseases"] = load_d2p(
        data["disease_to_phenotype_file"], data["hpo_network"], data["alt2prim"]
    )

    data["num_diseases_annotated"] = len(data["disease_records"])
    data["hpo_network"] = annotate(
        data["hpo_network"],
        data["phenotype_to_diseases"],
        data["num_diseases_annotated"],
        data["alt2prim"],
    )

    data["scorer"] = Scorer(data["hpo_network"], min_score_mask=None)
    data["disease_to_phenotype_output_file"] = os.path.join(
        data["parent_dir"], "data/phenotype.noparents.hpoa"
    )

    data["phenotype_groups"] = read_phenotype_groups()

    return data
