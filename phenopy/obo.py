import re

import obonet
import sys

import networkx as nx

from phenopy.ic import calculate_information_content


def load(obo_file, logger=None):
    """
    Load OBO file into a networkx graph.

    :param obo_file: OBO definition file.
    :param logger: Python `logging` logger instance.
    :return: `networkx.MultiDiGraph`
    """
    try:
        return obonet.read_obo(obo_file)
    except (FileNotFoundError, PermissionError) as e:
        if logger is not None:
            logger.critical(e)
        else:
            sys.stderr.write(str(e))
        exit(1)

def process(hpo_network, terms_to_genes, num_genes_annotated, custom_annotations_file=None, ages=None, logger=None):
    """
    Cleans the HPO network.

    Removes non-phenotype branches of the network, and merges all synonyms into one tag.

    :param hpo_network: `networkx.MultiDiGraph` to clean.
    :param terms_to_genes: Dictionary mapping HPO terms to genes.
    :param num_genes_annotated: Number of genes with HPO annotations.
    :param custom_annotations_file: A list of custom annotation files, in the same format as tests/data/test.score-product.txt
    :return: `networkx.MultiDiGraph`
    """

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

    # Before calculating information content, check for custom_annotations_file and load
    if custom_annotations_file is not None:
        custom_annos = {}
        try:
            with open(custom_annotations_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    entity_id, hpo_ids_string = line.split('\t')
                    hpo_ids = hpo_ids_string.split('|')
                    for hpo_id in hpo_ids:
                        if hpo_id in custom_annos:
                            custom_annos[hpo_id].append(entity_id)
                        else:
                            custom_annos[hpo_id] = [entity_id]

        except (FileNotFoundError, PermissionError, IsADirectoryError) as e:
            logger.critical(f'{custom_annotations_file} not found or is not accessible.')
            raise e
    else:
        custom_annos = None

    for node_id, data in hpo_network.nodes(data=True):
        # annotate with information content value
        hpo_network.node[node_id]['ic'] = calculate_information_content(
            node_id,
            hpo_network,
            terms_to_genes,
            num_genes_annotated,
            custom_annos,
        )

        hpo_network.node[node_id]['weights'] = {
            'age_dist': None,
            'age_exists': False,
            'default_weight': 1.0
        }
        if ages is not None:
            if node_id in ages.index:
                hpo_network.node[node_id]['weights']['age_dist'] = ages.loc[node_id]['age_dist']
                hpo_network.node[node_id]['weights']['age_exists'] = True

        # annotate with depth value
        # hard-coding origin node for now
        origin = 'HP:0000001'
        hpo_network.node[node_id]['depth'] = nx.shortest_path_length(
            hpo_network,
            node_id,
            origin
        )

        # clean synonyms
        synonyms = []
        try:
            for synonym in data['synonym']:
                synonyms.append(synonym)
            hpo_network.node[node_id]['synonyms'] = re.findall(r'"(.*?)"', ','.join(synonyms))
        except KeyError:
            # pass if no synonym tags in the node
            pass

    return hpo_network


def cache(hpo_network, cache_file):
    """
    Write a cache of a network as a pickle file.

    :param hpo_network: network to cache.
    :param cache_file: file path to use for cache.
    """
    nx.write_gpickle(hpo_network, cache_file)


def restore(cache_file):
    """
    Restores a network from a pickled cache file.

    :param cache_file: cache file path.
    :return: `networkx.MultiDiGraph`
    """
    return nx.read_gpickle(cache_file)
