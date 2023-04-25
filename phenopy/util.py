import csv
import os
import sys
import networkx as nx
import pandas as pd

from collections import Counter
from contextlib import contextmanager

from phenopy.config import config, logger


def half_product(num_rows, num_columns):
    """yield combinations and the diagonal"""
    for m in range(0, num_rows):
        for n in range(m, num_columns):
            yield (m, n)


def export_phenotype_hpoa_with_no_parents(phenotype_hpoa_file, phenotype_hpoa_no_parents_file, hpo_network, logger=None):
    """
    Load HPO terms associated to genes as annotated in https://hpo.jax.org/app/download/annotation.
    Filter the parent terms for each gene.
    Dump pheno2genes_no_parents_file

    :param phenotype_hpoa_file: Phenotypes to diseases file.
    :param phenotype_hpoa_no_parents_file: Phenotypes to diseases file with parents removed.
    :param hpo_network: The HPO networkx object.
    :param logger: Python `logging` logger instance.
    :return: None
    """
    try:
        with open(phenotype_hpoa_file, 'r') as tsv_fh:
            # skip the comment lines
            [next(tsv_fh) for _ in range(4)]
            df = pd.read_csv(
                tsv_fh,
                sep='\t',
            )
    except (FileNotFoundError, PermissionError) as e:
        if logger is not None:
            logger.critical(e)
        else:
            sys.stderr.write(str(e))
        exit(1)

    no_parents_df = df.copy()
    for gene, annotations in df.groupby('database_id'):
        termlist = [node for node in annotations['hpo_id'].tolist() if node in hpo_network.nodes()]
        termlist = remove_parents(termlist, hpo_network)
        parent_idx = annotations.loc[~annotations['hpo_id'].isin(termlist)].index
        no_parents_df.drop(parent_idx, inplace=True)

    try:
        no_parents_df.to_csv(phenotype_hpoa_no_parents_file,
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
    return sorted(set(termlist) - terms_to_remove)


def generate_alternate_ids(hpo_network):
    """Create a key, value store of alternate terms to canonical terms."""
    alt2prim = {}
    for n in hpo_network.nodes(data=True):
        n = n[0]
        try:
            for alt in hpo_network.nodes[n]['alt_id']:
                alt2prim[alt] = n
        except KeyError:
            # no alternate HPO ids for this term
            continue
    return alt2prim


def parse_input(input_file, hpo_network, alt2prim):
    """
    Parse input file.
    """
    try:
        with open(input_file, 'r') as input_fh:
            reader = csv.reader(filter(lambda l: not l.startswith('#'), input_fh), delimiter='\t')
            records = []
            for line in reader:
                # prcoess terms with convert and filter first
                terms = []
                for term_id in line[2].split('|'):
                    # convert alternate ids to primary
                    if term_id in alt2prim:
                        term_id = alt2prim[term_id]
                    # filtering terms not in the hpo network
                    if term_id not in hpo_network.nodes():
                        continue
                    terms.append(term_id)

                record = {
                    'record_id': line[0],
                    'terms': remove_parents(terms, hpo_network),
                    'weights': {},
                    **dict(item.split('=') for item in line[1].split(';') if line[1] != '.')
                }

                # assign new weights here ex. Sex weights (similar to the age weights).
                records.append(record)

    except (FileNotFoundError, PermissionError) as e:
        logger.critical(f'Provided input file could not be loaded or does not exist: {e}')
        exit(1)
    except ValueError:
        logger.critical(f'Unable to parse input file, invalid line number: {reader.line_num}:{input_file}')
        exit(1)

    return records


def read_phenotype_groups(phenotype_group_file=None):
    """
    Reads the phenotype group mappping file into a dictionary.
    :param phenotype_group_file: Filepath to the phenotype group file.
    """
    if phenotype_group_file is None:
        phenotype_group_file = config["phenotype_groups"]["phenotype_groups_file"]
    
    hp_to_pg = {}
    with open(phenotype_group_file, "r") as f:
        header = f.readline()
        for line in f:
            hpid, phenotype_group_1000, phenotype_group_1500 = line.strip("\n").split("\t")
            hp_to_pg[hpid] = {
                'k1000': int(phenotype_group_1000),
                'k1500': int(phenotype_group_1500),
            }
    return hp_to_pg


def standardize_phenotypes(terms, hpo_network, alt2prim):
    """Given a list of HPO ids, first try to convert synonyms to primary ids,
    then filter if terms are not in the ontology"""
    terms = [alt2prim[term] if term in alt2prim else term for term in terms]
    terms = list(filter(lambda term: term in hpo_network.nodes, terms))
    terms = remove_parents(terms, hpo_network)
    return terms


def encode_phenotypes(phenotypes, phenotype_groups, hpo_network, alt2prim, k=1000):
    """
    :param phenotypes:
    :param phenotype_groups: A dictionary 
    :return: a numpy array of one-hot encoded feature counts
    """

    def build_feature_array(cntr, n_features=k):
        a = [0] * (n_features)
        for feature_index, count in cntr.items():
            a[feature_index] = count
        return a

    nested = all(isinstance(element, list) for element in phenotypes)

    encode = lambda hpo_ids: Counter(hpo_ids)
    if nested:
        return [build_feature_array(encode([phenotype_groups[hpoid][f'k{k}'] for hpoid in standardize_phenotypes(phenotypes_, hpo_network, alt2prim)])) for phenotypes_ in phenotypes]
    
    return build_feature_array(encode([phenotype_groups[hpoid][f'k{k}'] for hpoid in standardize_phenotypes(phenotypes, hpo_network, alt2prim)]))


@contextmanager
def open_or_stdout(filename):
    if filename != '-':
        with open(filename, 'w') as f:
            yield f
    else:
        yield sys.stdout
