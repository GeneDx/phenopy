import networkx as nx
import obonet
import re
import sys

from phenopy.config import logger
from phenopy.ic import calculate_information_content
from phenopy.weights import make_age_distributions
from phenopy.util import parse_input


def load(obo_file):
    """
    Load OBO file into a networkx graph.

    :param obo_file: OBO definition file.
    :param logger: Python `logging` logger instance.
    :return: `networkx.MultiDiGraph`
    """
    try:
        hpo_network = obonet.read_obo(obo_file)
        #return nx.MultiDiGraph(hpo_network.subgraph(['HP:0000118'] + list(nx.ancestors(hpo_network, 'HP:0000118'))))
    except (FileNotFoundError, PermissionError) as e:
        if logger is not None:
            logger.critical(e)
        else:
            sys.stderr.write(str(e))
        exit(1)

    # roots for non-phenotype nodes
    non_phenotypes = {
        'mortality_aging': 'HP:0040006',
        'mode_of_inheritance': 'HP:0000005',
        'clinical_modifier': 'HP:0012823',
        'frequency': 'HP:0040279',
        'clinical_course': 'HP:0031797',
    }

    # remove non-phenotype branches
    for _, hpo_id in non_phenotypes.items():
        if hpo_id in hpo_network.nodes:
            children = nx.ancestors(hpo_network, hpo_id)
            hpo_network.remove_nodes_from([hpo_id] + list(children))

    return hpo_network


def annotate(hpo_network, phenotype_to_diseases, num_diseases_annotated, alt2prim, annotations_file=None, ages_distribution_file=None,
            phenotype_disease_frequencies=None):
    """
    Cleans the HPO network.

    Removes non-phenotype branches of the network, and merges all synonyms into one tag.

    :param hpo_network: `networkx.MultiDiGraph` to clean.
    :param phenotype_to_diseases: Dictionary mapping HPO terms to diseases.
    :param num_diseases_annotated: Number of diseases with HPO annotations.
    :param alt2prim: The dict of alternate terms to canonical terms.
    :param annotations_file: A list of custom annotation files, in the same format as tests/data/test.score-long.txt
    :param ages: age distributions object
    :param phenotype_disease_frequencies: dictionary of phenotype to disease frequencies
    :param logger: Python `logging` logger instance.
    :param ages_distribution_file: Path to phenotypes ages distribution file.
    :return: `networkx.MultiDiGraph`
    """

    # Before calculating information content, check for custom_annotations_file and load
    custom_annos = None
    if annotations_file is not None:
        custom_annos = {}
        for record in parse_input(annotations_file, hpo_network, alt2prim):
            for term_id in record['terms']:
                if term_id not in custom_annos:
                    custom_annos[term_id] = []
                custom_annos[term_id].append(record['record_id'])

    # make ages distributions
    ages = None
    if ages_distribution_file is not None:
        try:
            ages = make_age_distributions(ages_distribution_file)
            logger.info(
                f'Adding custom phenotype age distributions to HPO nodes from file: {ages_distribution_file}'
            )
        except (FileNotFoundError, PermissionError) as e:
            logger.critical(e)
            logger.critical(
                f'Specified phenotype ages file could not be loaded or does not exist: {e}'
            )
            exit(1)

    for node_id, data in hpo_network.nodes(data=True):
        # annotate with information content value
        hpo_network.nodes[node_id]['ic'] = calculate_information_content(
            node_id,
            hpo_network,
            phenotype_to_diseases,
            num_diseases_annotated,
            custom_annos,
        )
        # annotate with phenotype age distribution
        hpo_network.nodes[node_id]['disease_weights'] = {}

        if ages is not None and node_id in ages.index:
            hpo_network.nodes[node_id]['age_dist'] = ages.loc[node_id]['age_dist']

        # add the disease_frequency weights as attributes to the node
        if phenotype_disease_frequencies is not None:
            if node_id in phenotype_disease_frequencies:
                for disease_id, frequency in phenotype_disease_frequencies[node_id].items():
                    hpo_network.nodes[node_id]['weights']['disease_frequency'][disease_id] = frequency

        # annotate with depth value
        # hard-coding origin node for now
        origin = 'HP:0000001'
        hpo_network.nodes[node_id]['depth'] = nx.shortest_path_length(
            hpo_network,
            node_id,
            origin
        )

        # clean synonyms
        synonyms = []
        try:
            for synonym in data['synonym']:
                synonyms.append(synonym)
            hpo_network.nodes[node_id]['synonyms'] = re.findall(r'"(.*?)"', ','.join(synonyms))
        except KeyError:
            # pass if no synonym tags in the node
            pass

    return hpo_network
