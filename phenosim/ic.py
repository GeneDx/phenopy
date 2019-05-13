import networkx as nx
import numpy as np

SMOOTH = 1


def calculate_information_content(hpo_id, hpo_network, terms_to_genes, a_count):
    """
    Calculates information content for an HPO term.

    :param hpo_id: HPO term to calculate information content for.
    :param hpo_network: HPO network.
    :param terms_to_genes: HPO terms to genes mapping dictionary.
    :param a_count: Total count of HPO genes annotations.
    :return: `float`
    """
    # compile list of HPO terms to include in the calculation, term plus children
    hpo_id_plus_children = [hpo_id] + list(nx.ancestors(hpo_network, hpo_id))

    # count hpo genes annotations related to this term and children
    count_hpo_annos = sum([len(terms_to_genes[hpo_id]) for hpo_id in hpo_id_plus_children if hpo_id in terms_to_genes])

    # negative log of the number of hpo annotations divided by the total number of hpo terms in the
    # phenotypes_to_genes file
    return -np.log((count_hpo_annos + SMOOTH) / float(a_count + SMOOTH))
