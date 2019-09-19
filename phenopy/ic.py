import networkx as nx
import numpy as np

SMOOTH = 1


def calculate_information_content(hpo_id, hpo_network, terms_to_genes, num_genes_annotated, custom_annotations):
    """
    Calculates information content for an HPO term.

    :param hpo_id: HPO term to calculate information content for.
    :param hpo_network: HPO network.
    :param terms_to_genes: HPO terms to genes mapping dictionary.
    :param num_genes_annotated: Number of genes with HPO annotations.
    :param custom_annotations: A list of dictionaries mapping HPO terms to custom annotations
    :return: `float`
    """
    # compile list of HPO terms to include in the calculation, term plus children
    hpo_id_plus_children = [hpo_id] + list(nx.ancestors(hpo_network, hpo_id))
    # num_genes_annotated is the total number of genes in the annotation corpus.
    def get_ic(hpo_ids, annotations):
        # count the number of unique genes annotated to the hpo term and it's children
        n_unique_genes = len({g for h in hpo_ids if h in annotations for g in annotations[h]})
        # negative log of the number of hpo annotations divided by the total number of hpo terms in the
        # phenotypes_to_genes file
        information_content = -np.log((n_unique_genes + SMOOTH) / float(num_genes_annotated + SMOOTH))
        # include this to avoid returning zero division warnings.
        # when comparing two different leaves whose LCA is "phenotypic abnormality", the right side of the equation:
        # alphaIC / (alphaIC / betaIC), would be 0.0 / (0.0 / 0.0)
        if np.isnan(information_content):
            return np.nextafter(0, 1)
        elif np.isinf(information_content):
            return -np.log(np.nextafter(0, 1))
        return information_content

    annotations_list = [terms_to_genes]
    if custom_annotations is not None:
        annotations_list.append(custom_annotations)

    return np.mean([get_ic(hpo_id_plus_children, annotations=annotations) for annotations in annotations_list])

