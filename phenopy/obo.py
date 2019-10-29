import re
import obonet
import sys
import networkx as nx

from phenopy import parse_input
from phenopy.ic import calculate_information_content





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
