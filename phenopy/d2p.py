import numpy as np
import csv
import sys

hpo_id_to_float = {
    'HP:0040280': 1.0,
    'HP:0040281': np.mean([80, 99]),
    'HP:0040282': np.mean([0.30, 0.79]),
    'HP:0040283': np.mean([0.05, 0.29]),
    'HP:0040284': np.mean([0.01, 0.04]),
    'HP:0040285': 0.0,
}


def load(phenotype_annotations_file, hpo_network, alt2prim, logger=None):
    """Parse the hpoa file
    :param phenotype_annotations_file: path to the phenotype.hpoa file
    :return: three dictionaries of disease to phenotypes, phenotypes to disease, and phenotypes to disease frequencies
    """
    try:
        with open(phenotype_annotations_file, 'r') as tsv_fh:
            reader = csv.DictReader(filter(lambda line: line[0] != '#', tsv_fh), delimiter='\t')

            disease_to_phenotypes = dict()
            phenotype_to_diseases = dict()
            for row in reader:
                # phenotype term id
                # convert alternate phenotype id to primary
                term_id = row['HPO_ID'] if row['HPO_ID'] not in alt2prim else alt2prim[row['HPO_ID']]

                if term_id not in hpo_network.nodes():
                    continue

                # parse disease id, currently only supports omim entries
                db, disease_accession = row['DatabaseID'].split(':')
                if db not in ['OMIM']:
                    continue

                # For now, skip negative phenotype annotations
                if row['Qualifier'] == 'NOT':
                    continue

                if term_id not in phenotype_to_diseases:
                    phenotype_to_diseases[term_id] = {disease_accession: {'frequency': []}}

                else:
                    if disease_accession not in phenotype_to_diseases[term_id]:
                        phenotype_to_diseases[term_id].update({disease_accession: {'frequency': []}})

                phenotype_to_diseases[term_id][disease_accession]['frequency'].append(frequency_converter(row['Frequency']))

                # add the phenotype to the disease in the disease_records dictionary
                if disease_accession not in disease_to_phenotypes:
                    disease_to_phenotypes[disease_accession] = {'record_id': disease_accession,
                                                                'terms': [],
                                                                'weights': {'disease_frequency': [],
                                                                            },
                                                                }
                disease_to_phenotypes[disease_accession]['terms'].append(term_id)

        # going from dict to a list of disease records and setting weights
        disease_records = list()
        for disease_accession, disease in disease_to_phenotypes.items():
            disease['terms'] = sorted(set(disease['terms']))
            for term_id in disease['terms']:
                try:
                    frequency_weight = np.mean(phenotype_to_diseases[term_id][disease_accession]['frequency'])
                except TypeError:
                    # TODO: discuss with team what is a good default for the unannotated disease frequency.
                    frequency_weight = 0.5
                disease['weights']['disease_frequency'].append(frequency_weight)
            disease_records.append(disease)

        # TODO: do we need phenotype_to_diseases
        return disease_records, phenotype_to_diseases

    except (FileNotFoundError, PermissionError) as e:
        hpoa_file_error_msg = f'{phenotype_annotations_file} not found or incorrect permissions'
        if logger is not None:
            logger.critical(hpoa_file_error_msg)
        else:
            sys.stderr.write(hpoa_file_error_msg)
        exit(1)


def frequency_converter(hpoa_frequency):
    """convert the frequency column from the hpoa file to a float"""
    if 'HP:' in hpoa_frequency:
        return hpo_id_to_float.get(hpoa_frequency, '')

    elif '/' in hpoa_frequency:
        n, d = hpoa_frequency.split('/')
        return float(n) / float(d)

    elif '%' in hpoa_frequency:
        return float(hpoa_frequency.strip('%')) / 100
    # return 0.5 by default
    return None
