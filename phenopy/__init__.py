__project__ = 'phenopy'
__version__ = '0.4.0'

import sys
from contextlib import contextmanager
from phenopy.d2p import load as load_d2p
from phenopy.network import load as load_network
from phenopy.network import annotate
from phenopy.util import generate_alternate_ids


@contextmanager
def open_or_stdout(filename):
    if filename != '-':
        with open(filename, 'w') as f:
            yield f
    else:
        yield sys.stdout


def generate_annotated_hpo_network(obo_file, disease_to_phenotype_file, annotations_file=None, ages_distribution_file=None):
    hpo_network = load_network(obo_file)

    alt2prim = generate_alternate_ids(hpo_network)

    # load phenotypes to diseases associations
    (
        disease_records,
        phenotype_to_diseases,
    ) = load_d2p(disease_to_phenotype_file, hpo_network, alt2prim)

    # load hpo network
    hpo_network = annotate(
        hpo_network,
        phenotype_to_diseases,
        len(disease_records),
        alt2prim,
        annotations_file=annotations_file,
        ages_distribution_file=ages_distribution_file,
    )

    return hpo_network, alt2prim, disease_records