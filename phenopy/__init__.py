__project__ = 'phenopy'
__version__ = '0.2.1'

import csv
import sys
from contextlib import contextmanager

from phenopy.config import logger
from phenopy.d2p import load as load_d2p
from phenopy.network import load as load_network
from phenopy.network import annotate
from phenopy.util import remove_parents
from phenopy.weights import calculate_age_weights


@contextmanager
def open_or_stdout(filename):
    if filename != '-':
        with open(filename, 'w') as f:
            yield f
    else:
        yield sys.stdout


def generate_alternate_ids(hpo_network):
    """Create a key, value store of alternate terms to canonical terms."""
    alt2prim = {}
    for n in hpo_network.nodes(data=True):
        n = n[0]
        try:
            for alt in hpo_network.node[n]['alt_id']:
                alt2prim[alt] = n
        except KeyError:
            # no alternate HPO ids for this term
            continue
    return alt2prim


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


def parse_input(input_file, hpo_network, alt2prim):
    """
    Parse input file.
    """
    try:
        with open(input_file, 'r') as input_fh:
            reader = csv.reader(filter(lambda l: not l.startswith('#'), input_fh), delimiter='\t')
            records = []
            for line in reader:
                # prcoess terms with convert and filter first
                terms = []
                for term_id in line[2].split('|'):
                    # convert alternate ids to primary
                    if term_id in alt2prim:
                        term_id = alt2prim[term_id]
                    # filtering terms not in the hpo network
                    if term_id not in hpo_network.nodes():
                        continue
                    terms.append(term_id)

                record = {
                    'record_id': line[0],
                    'terms': remove_parents(terms, hpo_network),
                    'weights': {},
                    **dict(item.split('=') for item in line[1].split(';') if line[1] != '.')
                }

                # set weights
                if 'age' in record:
                    record['weights']['age'] = calculate_age_weights(record['terms'], record['age'], hpo_network)
                # assign new weights here ex. Sex weights (similar to the age weights).
                records.append(record)

    except (FileNotFoundError, PermissionError) as e:
        logger.critical(f'Provided input file could not be loaded or does not exist: {e}')
        exit(1)
    except ValueError:
        logger.critical(f'Unable to parse input file, invalid line number: {reader.line_num}:{input_file}')
        exit(1)

    return records
