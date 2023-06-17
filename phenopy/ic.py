import networkx as nx
import numpy as np
from typing import Dict

SMOOTH = 1


def calculate_information_content(
        hpo_id: str,
        hpo_network: nx.MultiDiGraph,
        phenotype_to_diseases: Dict,
        num_diseases_annotated: int,
        custom_annotations: str = None) -> np.ndarray:
    """
    Calculates information content for an HPO term.
    """
    # compile list of HPO terms to include in the calculation, term plus children
    hpo_id_plus_children = [hpo_id] + list(nx.ancestors(hpo_network, hpo_id))
    # num_diseases_annotated is the total number of diseases in the annotation corpus.

    def get_ic(hpo_ids, annotations):
        # count the # of unique diseases annotated to the hpo term and it's children
        n_unique_diseases = len(
            {g for h in hpo_ids if h in annotations for g in annotations[h]}
        )
        # negative log of the number of hpo annotations divided by the total number
        # of hpo terms in the
        # phenotypes_to_genes file
        information_content = -np.log((n_unique_diseases + SMOOTH) /
                                      float(num_diseases_annotated + SMOOTH))

        return information_content

    annotations_list = [phenotype_to_diseases]
    if custom_annotations is not None:
        annotations_list.append(custom_annotations)
    output_mean = np.mean([get_ic(hpo_id_plus_children, annotations=annotations)
                           for annotations in annotations_list])
    return output_mean
