import sys

import pandas as pd


def load(pheno2genes_file, logger=None):
    """
    Load HPO terms associated to genes as annotated in https://hpo.jax.org/app/download/annotation.

    :param pheno2genes_file: Phenotypes to genes file.
    :param logger: Python `logging` logger instance.
    :return: `dict` {hpo_term: [genes]}, `dict` {gene: [hpo_ids]}, `int` number of genes  annotations
    """
    try:
        df = pd.read_csv(
            pheno2genes_file,
            comment='#',
            names=['hpo_id','gene_id','gene_name'],
            usecols=[0,2,3],
            sep='\t',
        )
    except (FileNotFoundError, PermissionError) as e:
        if logger is not None:
            logger.critical(e)
        else:
            sys.stderr.write(str(e))
        exit(1)

    # save map of terms to genes, and keep count of total number of annotations
    terms_to_genes = {}
    count = 0

    # convert gene id to string
    df['gene_id'] = df['gene_id'].astype(str)

    # make unique gene identifier
    df['gene'] = 'NCBI:' + df['gene_id'] + '[' + df['gene_name'] + ']'

    for term, annotations in df.groupby('hpo_id'):
        terms_to_genes[term] = annotations['gene'].tolist()
        count += len(terms_to_genes[term])

    # save map of genes to terms
    genes_to_terms = {gene: annotations['hpo_id'].tolist(
    ) for gene, annotations in df.groupby('gene')}

    return terms_to_genes, genes_to_terms, len(genes_to_terms)
