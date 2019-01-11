import itertools

import networkx as nx
import pandas as pd


def calculate_gamma(term_a, term_b, term_lcafn, hpo_network):
    """
    TODO: Documentation?

    :param term_a: HPO term A.
    :param term_b: HPO term B.
    :param term_lcafn: TODO: Documentation?
    :param hpo_network: HPO network.
    :return: `float` (term pair comparison score)
    """
    # calculate gamma?
    if term_a == term_b:
        return 0
    elif term_b in nx.ancestors(hpo_network, term_a):
        if nx.shortest_path_length(hpo_network, term_b, term_a) == 1:
            return 1

    a_to_lcafn = nx.shortest_path_length(hpo_network, term_a, term_lcafn)
    b_to_lcafn = nx.shortest_path_length(hpo_network, term_b, term_lcafn)

    return a_to_lcafn + b_to_lcafn


def score_hpo_pair_hrss(term_a, term_b, hpo_network):
    """
    Scores the comparison of a pair of terms, using HRSS algorithm.
    
    :param term_a: HPO term A.
    :param term_b: HPO term B.
    :param hpo_network: HPO network.
    :return: `float` (term pair comparison score)
    """
    # find mils? information content for each term
    mil_ic = []
    for term in [term_a, term_b]:
        if hpo_network.in_edges(term):
            mil_ic.append(max(
                {hpo_network.node[p]['ic'] for p in hpo_network.predecessors(term) if 'ic' in hpo_network.node[p]}
            ))
        else:
            mil_ic.append(hpo_network.node[term]['ic'])

    # calculate beta_ic?
    beta_ic = ((mil_ic[0] - hpo_network.node[term_a]['ic']) + (mil_ic[1] - hpo_network.node[term_b]['ic'])) / 2.0

    # find common bfs? predecessors
    bfs_predecessors = []
    for term in [term_a, term_b]:
        bfs_predecessors.append({p[0] for p in nx.bfs_predecessors(hpo_network, term)})
    common_bfs_predecessors = bfs_predecessors[0].intersection(bfs_predecessors[1])

    # lcafn? node
    lcafn_node = max(common_bfs_predecessors, key=lambda n: hpo_network.node[n]['depth'])

    # calculate alpha_ic?
    alpha_ic = hpo_network.node[lcafn_node]['ic']

    # calculate gamma?
    gamma = calculate_gamma(term_a, term_b, lcafn_node, hpo_network)

    return (1.0 / float(1.0 + gamma)) * (alpha_ic / float(alpha_ic + beta_ic))


def score(terms_a, terms_b, hpo_network):
    """
    Scores the comparison of terms in list A to terms in list B.

    :param terms_a: List of HPO terms A.
    :param terms_b: List of HPO terms B.
    :param hpo_network: network to cache.
    :return: `float` (comparison score)
    """
    df = pd.DataFrame(
        [(a, b, score_hpo_pair_hrss(a, b, hpo_network)) for a, b in itertools.product(terms_a, terms_b)],
        columns=['a', 'b', 'score']
    ).set_index(
        ['a', 'b']
    ).unstack()

    return round((df.max(axis=1).mean() + df.max(axis=0).mean()) / 2.0, 4)
