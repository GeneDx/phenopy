import csv
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
        termlist = remove_parents(termlist, hpo_network)
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


def parse(string, what='HPO'):
    """
    Parse patient parameters in the records file
    :param string: string to parse
    :param what: (HP,age,sex) terms to parse
    :return: parsed object, int for age, string for gender, list for terms
    """
    string = string.strip()
    if string == '.':
        return None
    if what == 'HPO':
        result = [x for x in string.split('|') if x.startswith('HP:')]
        return result
    elif f'{what}=' in string:
        result = [x.split(f'{what}=')[1] for x in string.split(';') if what in x]
        if result:
            result = result[0]
            if what == 'age':
                try:
                    result = round(float(result), 1)
                except ValueError:
                    result = None

            if what == 'sex':
                if result.lower().startswith('f'):
                    result = 'Female'
                elif result.lower().startswith('m'):
                    result = 'Male'
                else:
                    result = None
            return result
        else:
            return None


def read_records_file(records_file, no_parents=False, hpo_network=None, logger=None):
    """
    Parse input file for patient descriptions into an array of dictionaries
    :param records_file: path to the records file to parse
    :param no_parents: remove parent nodes
    :param hpo_network: hpo network to use in removing parents
    :param logger: logger object to use in reporting errors
    :return: list of dictionaries
    """
    try:
        with open(records_file) as records_fh:
            reader = csv.reader(records_fh, delimiter='\t')
            records = []
            for line in reader:
                if line[0].startswith('#'):
                    continue
                dict_ = {
                    'sample': line[0],
                    'age': parse(line[1], what='age'),
                    'gender': parse(line[1], what='sex'),
                    'terms': parse(line[2], what='HPO')
                }

                if no_parents is True and hpo_network is not None:
                    dict_['terms'] = remove_parents(dict_['terms'], hpo_network)
                else:
                    pass
                records.append(dict_)
        return records
    except (FileNotFoundError, PermissionError) as e:
        if logger is not None:
            logger.critical(e)
        else:
            sys.stderr.write(str(e))
        exit(1)


def remove_parents(termlist, hpo_network):
    """remove parents from termlist
    :param termlist: List of HPO terms.
    :param hpo_network: The HPO networkx object.
    """
    terms_to_remove = set()
    for source_term in termlist:
        if source_term not in hpo_network.nodes:
            terms_to_remove.add(source_term)
            continue
        for target_term in termlist:
            if target_term not in hpo_network.nodes:
                terms_to_remove.add(target_term)
                continue
            # has_path will evaluate True for a term to itself, include additional check
            same_terms = source_term == target_term
            source_to_target = nx.has_path(hpo_network, source_term, target_term)
            target_to_source = nx.has_path(hpo_network, target_term, source_term)
            if source_to_target is True and not same_terms:
                terms_to_remove.add(target_term)
            if target_to_source is True and not same_terms:
                terms_to_remove.add(source_term)
    return set(termlist) - terms_to_remove
