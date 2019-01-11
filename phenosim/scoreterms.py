import pandas

import networkx as nx


# where is hpo_network coming from in each of these???


def lcafn(hpo_network, p, q):
    """
    Determine the lowest common ancestor for a two terms

    :param hpo_network: `networkx.MultiDiGraph`
    :param term1: First HPO term, "HP:0000001"
    :param term2: Second HPO term, "HP:0001263"
    :return: Least Common Ancestor for two terms, "HP:0000001"
    """

    preds_1 = dict(nx.bfs_predecessors(hpo_network, p))
    preds_2 = dict(nx.bfs_predecessors(hpo_network, q))

    common_preds = set([n for n in preds_1]).intersection(
        set([n for n in preds_2]))

    return max(common_preds, key=lambda n: hpo_network.node[n]['depth'])


def BetaIC(hpo_network, term1, term2):
    """
    Beta Information Content
    :param hpo_network: `networkx.MultiDiGraph`
    :param term1: First HPO term, "HP:0000001"
    :param term2: Second HPO term, "HP:0001263"
    :return: BetaIC for two terms, <float>
    """
    if list(hpo_network.in_edges(term1)) == []:
        MILi = term1
    else:
        MILi = max({term for term in list(hpo_network.predecessors(
            term1)) if term in hpo_network['ic'].keys()}, key=lambda n: hpo_network['ic'][n])

    if list(hpo_network.in_edges(term2)) == []:
        MILj = term2
    else:
        MILj = max({term for term in list(hpo_network.predecessors(
            term2)) if term in hpo_network['ic'].keys()}, key=lambda n: hpo_network['ic'][n])
    return ((hpo_network['ic'][MILi] - hpo_network['ic'][term1]) + (hpo_network['ic'][MILj] - hpo_network['ic'][term2])) / 2.0


def Gamma(hpo_network, term1, term2):
    """
    Gamma term
    :param hpo_network: `networkx.MultiDiGraph`
    :param term1: First HPO term, "HP:0000001"
    :param term2: Second HPO term, "HP:0001263"
    :return: Gamma for two terms, <float>
    """
    if term1 == term2:
        return 0

    elif term2 in nx.ancestors(hpo_network, term1):
        if nx.shortest_path_length(hpo_network, term2, term1) == 1:
            return 1

    return nx.shortest_path_length(hpo_network, term1,lcafn(term1, term2)) + nx.shortest_path_length(hpo_network, term2, lcafn(term1, term2))


def HRSS(hpo_network, term1, term2):
    """
    Main phenotype similarity scoring function.
    HRSS algorithm from https://www.ncbi.nlm.nih.gov/pubmed/23741529

    :param hpo_network: `networkx.MultiDiGraph`
    :param term1: First HPO term, "HP:0000001"
    :param term2: Second HPO term, "HP:0001263"
    :return: HRSS for two terms, <float>
    """

    betaIC = BetaIC(hpo_network, term1, term2)
    alphaIC = hpo_network['ic'][lcafn(hpo_network, term1, term2)]  # IC of lca
    gamma = Gamma(hpo_network, term1, term2)
    return round((1.0 / float(1.0 + gamma)) * (alphaIC / float(alphaIC + betaIC)), 4)


def HRSS_wrapper(terms):
    return HRSS(terms[0], terms[1])

def score_two_lists(hpo_network, tlist1, tlist2):
    """
    Score two lists of HPO terms by the HRSS alogorithm.

    :param hpo_network: `networkx.MultiDiGraph`
    :param tlist1: [hpo_ids]
    :param tlist2: [hpo_ids]
    """
    df = pandas.DataFrame(index=pandas.MultiIndex.from_product((tlist1, tlist2)))
    df['score'] = df.index.map(HRSS_wrapper).values
    df = df.unstack()
    # return BMA
    return round((df.max(axis=1).mean() + df.max(axis=0).mean()) / 2.0, 4)


def score_case_genes(hpo_network, case_terms, genes_to_terms):
    """
    Score the phenotype similarity of one case to all genes.

    :param hpo_network: `networkx.MultiDiGraph`
    :param case_terms: [hpo_ids]
    :param gene_to_terms {gene_name: [hpo_ids]}
    """
    score_holder = {}
    for gene_name, gene_terms in genes_to_terms.items():
        score_holder[gene_name] = score_two_lists(hpo_network, case_terms, gene_terms)
    return score_holder
