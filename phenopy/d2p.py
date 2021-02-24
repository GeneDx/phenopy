import numpy as np
import csv
import sys

hpo_id_to_float = {
    'HP:0040280': 1.0,
    'HP:0040281': np.mean([0.80, 0.99]),
    'HP:0040282': np.mean([0.30, 0.79]),
    'HP:0040283': np.mean([0.05, 0.29]),
    'HP:0040284': np.mean([0.01, 0.04]),
    'HP:0040285': 0.0,
}


def read_hpo_annotation_file(phenotype_annotations_file, hpo_network, logger=None):
    """
    :param phenotype_annotations_file: path to the phenotype.hpoa file
    :return: records
    """
    try:
        with open(phenotype_annotations_file, 'r') as tsv_fh:
            #[next(tsv_fh) for _ in range(3)]
            reader = csv.DictReader(tsv_fh, delimiter='\t')
            # this removes the leading hash
            reader.fieldnames[0] = 'disease-db'

            records = []

            for row in reader:
                # phenotype term id
                term_id = row['HPO-ID']
                if term_id not in hpo_network.nodes():
                    continue
                # parse disease id, currently only supports omim entries
                if row['disease-db'] not in ['OMIM']:
                    continue
                # For now, skip negative phenotype annotations
                if row['negation'] == 'NOT':
                    continue

                records.append((term_id, row["disease-identifier"], frequency_converter(row['frequencyHPO'])))

        return records

    except (FileNotFoundError, PermissionError) as e:
        hpoa_file_error_msg = f'{phenotype_annotations_file} not found or incorrect permissions'
        if logger is not None:
            logger.critical(hpoa_file_error_msg)
        else:
            sys.stderr.write(hpoa_file_error_msg)
        sys.exit(1)


def read_custom_annotation_file(custom_annotation_file_path, hpo_network, logger=None):
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

    except (FileNotFoundError, PermissionError) as e:
        hpoa_file_error_msg = f'{custom_annotation_file_path} not found or incorrect permissions'
        if logger is not None:
            logger.critical(hpoa_file_error_msg)
        else:
            sys.stderr.write(hpoa_file_error_msg)
        sys.exit(1)


def load(phenotype_annotations_file, hpo_network, alt2prim, default_frequency=0.5):
    """Parse the hpoa file
    :param phenotype_annotations_file: path to the phenotype.hpoa file
    :return: three dictionaries of disease to phenotypes, phenotypes to disease, and phenotypes to disease frequencies
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
            phenotype_to_diseases[term_id] = {disease_accession: {'frequency': default_frequency}}
        else:
            if disease_accession not in phenotype_to_diseases[term_id]:
                phenotype_to_diseases[term_id].update({disease_accession: {'frequency': default_frequency}})

        phenotype_to_diseases[term_id][disease_accession]['frequency'] = freq

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
            # convert alternate phenotype id to primary
            term_id = term_id if term_id not in alt2prim else alt2prim[term_id]
            if term_id not in hpo_network.nodes():
                continue

            frequency_weight = phenotype_to_diseases[term_id][disease_accession]['frequency']
            #
            disease['weights']['disease_frequency'].append(frequency_weight)

        disease_records.append(disease)

    # TODO: do we need phenotype_to_diseases?
    return disease_records, phenotype_to_diseases


def frequency_converter(hpoa_frequency, default_frequency=0.5):
    """convert the frequency column from the hpoa file to a float"""
    if 'HP:' in hpoa_frequency:
        #TODO discuss the best default
        return hpo_id_to_float.get(hpoa_frequency, default_frequency)

    elif '/' in hpoa_frequency:
        n, d = hpoa_frequency.split('/')
        return float(n) / float(d)

    elif '%' in hpoa_frequency:
        return float(hpoa_frequency.strip('%')) / 100

    # return 0.5 by default
    #TODO discuss the best default
    return default_frequency
