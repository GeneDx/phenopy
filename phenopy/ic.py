import networkx as nx
import numpy as np

SMOOTH = 1


def calculate_information_content(hpo_id, hpo_network, phenotype_to_diseases, num_diseases_annotated, custom_annotations):
    """
    Calculates information content for an HPO term.

    :param hpo_id: HPO term to calculate information content for.
    :param hpo_network: HPO network.
    :param phenotype_to_diseases: HPO terms to diseases mapping dictionary.
    :param num_diseases_annotated: Number of diseases with HPO annotations.
    :param custom_annotations: A list of dictionaries mapping HPO terms to custom annotations
    :return: `float`
    """
    # compile list of HPO terms to include in the calculation, term plus children
    hpo_id_plus_children = [hpo_id] + list(nx.ancestors(hpo_network, hpo_id))
    # num_diseases_annotated is the total number of diseases in the annotation corpus.

    def get_ic(hpo_ids, annotations):
        # count the number of unique diseases annotated to the hpo term and it's children
        n_unique_diseases = len({disease_id for hpo_id in hpo_ids if hpo_id in annotations for disease_id in annotations[hpo_id]})
        # negative log of the number of hpo annotations divided by the total number of hpo terms in phenotype.hpoa
        # additive smoothing using Laplace's Rule of succession
        # https://en.wikipedia.org/wiki/Additive_smoothing
        information_content = -np.log((n_unique_diseases + SMOOTH) / float(num_diseases_annotated))

        if information_content < 0.0:
            information_content = 0.0

        return information_content

    annotations_list = [phenotype_to_diseases]
    if custom_annotations is not None:
        annotations_list.append(custom_annotations)

    return np.mean([get_ic(hpo_id_plus_children, annotations=annotations) for annotations in annotations_list])

