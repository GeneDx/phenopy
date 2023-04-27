import numpy as np
import csv
import sys
import networkx as nx
from typing import (
    List,
    Tuple,
)

hpo_id_to_float = {
    'HP:0040280': 1.0,
    'HP:0040281': np.mean([0.80, 0.99]),
    'HP:0040282': np.mean([0.30, 0.79]),
    'HP:0040283': np.mean([0.05, 0.29]),
    'HP:0040284': np.mean([0.01, 0.04]),
    'HP:0040285': 0.0,
}


def read_hpo_annotation_file(phenotype_annotations_file: str,
                             hpo_network: nx.MultiDiGraph,
                             logger=None) -> List:
    """
    Reads the annotation files from the HPO website
    """
    try:
        with open(phenotype_annotations_file, 'r') as tsv_fh:
            [next(tsv_fh) for _ in range(4)]
            reader = csv.DictReader(tsv_fh, delimiter='\t')
            # this removes the leading hash
            reader.fieldnames[0] = 'database_id'

            records = []

            for row in reader:
                # phenotype term id
                term_id = row.get('HPO_ID') if 'HPO_ID' in row else row.get('hpo_id')
                if term_id not in hpo_network.nodes():
                    continue
                # parse disease id, currently only supports omim entries
                db, disease_accession = row['database_id'].split(':')
                if db not in ['OMIM']:
                    continue
                # For now, skip negative phenotype annotations
                if row['qualifier'] == 'NOT':
                    continue

                records.append(
                    (term_id, disease_accession, frequency_converter(row['frequency']))
                )

        return records

    except (FileNotFoundError, PermissionError):
        hpoa_file_error_msg = f"{phenotype_annotations_file} " \
                              f"not found or incorrect permissions"
        if logger is not None:
            logger.critical(hpoa_file_error_msg)
        else:
            sys.stderr.write(hpoa_file_error_msg)
        sys.exit(1)


def read_custom_annotation_file(
        custom_annotation_file_path: str,
        hpo_network: nx.MultiDiGraph,
        logger: None = None) -> List:
    try:
        with open(custom_annotation_file_path, 'r') as tsv_fh:
            reader = csv.reader(tsv_fh, delimiter='\t')

            records = []
            for row in reader:
                # phenotype term id
                # convert alternate phenotype id to primary
                term_id, disease_accession, freq = row
                if term_id not in hpo_network.nodes():
                    continue

                records.append((term_id, disease_accession, float(freq)))

        return records

    except (FileNotFoundError, PermissionError):
        hpoa_file_error_msg = f"{custom_annotation_file_path} " \
                              f"not found or incorrect permissions"
        if logger is not None:
            logger.critical(hpoa_file_error_msg)
        else:
            sys.stderr.write(hpoa_file_error_msg)
        sys.exit(1)


def load(
        phenotype_annotations_file: str,
        hpo_network: nx.MultiDiGraph,
        alt2prim,
        default_frequency: float = 0.5) -> Tuple:
    """
    Parse the hpoa file
    """
    if phenotype_annotations_file.endswith("hpoa"):
        records = read_hpo_annotation_file(phenotype_annotations_file, hpo_network)
    else:
        records = read_custom_annotation_file(phenotype_annotations_file, hpo_network)

    disease_to_phenotypes = dict()
    phenotype_to_diseases = dict()

    for r in records:
        term_id, disease_accession, freq = r
        if term_id not in phenotype_to_diseases:
            phenotype_to_diseases[term_id] = {
                disease_accession: {'frequency': default_frequency}
            }
        else:
            if disease_accession not in phenotype_to_diseases[term_id]:
                phenotype_to_diseases[term_id].update(
                    {disease_accession: {'frequency': default_frequency}}
                )

        phenotype_to_diseases[term_id][disease_accession]['frequency'] = freq

        # add the phenotype to the disease in the disease_records dictionary
        if disease_accession not in disease_to_phenotypes:
            disease_to_phenotypes[disease_accession] = {
                'record_id': disease_accession,
                'terms': [],
                'weights': {'disease_frequency': [], }
            }
        disease_to_phenotypes[disease_accession]['terms'].append(term_id)

    # going from dict to a list of disease records and setting weights
    disease_records = list()
    for disease_accession, disease in disease_to_phenotypes.items():
        disease['terms'] = sorted(set(disease['terms']))
        for term_id in disease['terms']:
            # convert alternate phenotype id to primary
            term_id = term_id if term_id not in alt2prim else alt2prim[term_id]
            if term_id not in hpo_network.nodes():
                continue

            frequency_weight = \
                phenotype_to_diseases[term_id][disease_accession]['frequency']
            #
            disease['weights']['disease_frequency'].append(frequency_weight)

        disease_records.append(disease)

    # TODO: do we need phenotype_to_diseases?
    return disease_records, phenotype_to_diseases


def frequency_converter(
        hpoa_frequency: str,
        default_frequency: float = 0.5) -> float:
    """
    convert the frequency column from the hpoa file to a float
    """
    if 'HP:' in hpoa_frequency:
        # TODO discuss the best default
        return hpo_id_to_float.get(hpoa_frequency, default_frequency)

    elif '/' in hpoa_frequency:
        n, d = hpoa_frequency.split('/')
        return float(n) / float(d)

    elif '%' in hpoa_frequency:
        return float(hpoa_frequency.strip('%')) / 100

    # TODO discuss the best default
    return default_frequency
