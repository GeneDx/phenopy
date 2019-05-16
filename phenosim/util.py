import sys
import networkx as nx
import pandas as pd


def export_pheno2genes_with_no_parents(pheno2genes_file, pheno2genes_no_parents_file, hpo_network, logger=None):
    """
    Load HPO terms associated to genes as annotated in https://hpo.jax.org/app/download/annotation.
    Filter the parent terms for each gene.
    Dump pheno2genes_no_parents_file

    :param pheno2genes_file: Phenotypes to genes file.
    :param pheno2genes_no_parents_file: Phenotypes to genes file with parents removed.
    :param hpo_network: The HPO networkx object.
    :param logger: Python `logging` logger instance.
    :return: None
    """
    try:
        df = pd.read_csv(
            pheno2genes_file,
            comment='#',
            names=['hpo_id', 'hpo_name', 'gene_id', 'gene_name'],
            sep='\t',
        )
    except (FileNotFoundError, PermissionError) as e:
        if logger is not None:
            logger.critical(e)
        else:
            sys.stderr.write(str(e))
        exit(1)

    no_parents_df = df.copy()
    for gene, annotations in df.groupby('gene_name'):
        termlist = [node for node in annotations['hpo_id'].tolist() if node in hpo_network.nodes()]
        termlist = remove_parent_termlist(termlist, hpo_network)
        parent_idx = annotations.loc[~annotations['hpo_id'].isin(termlist)].index
        no_parents_df.drop(parent_idx, inplace=True)

    try:
        no_parents_df.to_csv(pheno2genes_no_parents_file,
                             header=['#hpo_id', 'hpo_name', 'gene_id', 'gene_name'],
                             sep='\t',
                             index=False)
    except PermissionError as e:
        if logger is not None:
            logger.critical(e)
        else:
            sys.stderr.write(str(e))
        exit(1)


def remove_parent_termlist(termlist, hpo_network):
    """remove parents from termlist
    :param termlist: List of HPO terms.
    :param hpo_network: The HPO networkx object.
    """
    pruned_termlist = []
    for source_term in termlist:
        for target_term in termlist:
            # has_path will evaluate True for a term to itself, include additional check
            same_terms = source_term == target_term
            has_path = nx.has_path(hpo_network, source_term, target_term)
            if not has_path and not same_terms:
                pruned_termlist.append(source_term)
    return pruned_termlist
