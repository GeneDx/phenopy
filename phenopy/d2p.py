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


def load(phenotype_annotations_file, logger=None):
    """Parse the hpoa file
    :param phenotype_annotations_file: path to the phenotype.hpoa file
    :return: three dictionaries of disease to phenotypes, phenotypes to disease, and phenotypes to disease frequencies
    """
    try:
        with open(phenotype_annotations_file, 'r') as tsv_fh:
            reader = csv.DictReader(filter(lambda line: line[0] != '#', tsv_fh), delimiter='\t')

            disease_to_phenotypes = dict()
            phenotype_disease_frequencies = dict()
            for row in reader:
                # phenotype term id
                term_id = row['HPO_ID']

                # parse disease id, currently only supports omim entries
                db, disease_accession = row['DatabaseID'].split(':')
                if db not in ['OMIM']:
                    continue

                # For now, skip negative phenotpye annotations
                if row['Qualifier'] == 'NOT':
                    continue

                # annotate the frequency of the phenotype to disease
                # assign a new dictionary for a phenotype term when one doesn't exist yet
                if term_id not in phenotype_disease_frequencies:
                    phenotype_disease_frequencies[term_id] = {disease_accession:
                                                                  [frequency_converter(row['Frequency'])]}
                else:
                    # if the phenotype key already exists, check if the disease accession also already exists.
                    # there can be multiple reported disease phenotype frequencies in the .hpoa file.
                    if disease_accession in phenotype_disease_frequencies[term_id]:
                        phenotype_disease_frequencies[term_id][disease_accession].append(
                            frequency_converter(row['Frequency']))
                    else:
                        phenotype_disease_frequencies[term_id][disease_accession] = \
                            [frequency_converter(row['Frequency'])]

                # add the phenotype to the disease in the disease_to_phenotypes dictionary
                if disease_accession not in disease_to_phenotypes:
                    disease_to_phenotypes[disease_accession] = [term_id]
                else:
                    disease_to_phenotypes[disease_accession].append(term_id)

            for phenotype, diseases in phenotype_disease_frequencies.items():
                for disease in diseases:
                    phenotype_disease_frequencies[phenotype][disease] = \
                        np.mean(phenotype_disease_frequencies[phenotype][disease])

        phenotype_to_diseases = dict()
        for disease_accession, phenotype_ids in disease_to_phenotypes.items():
            for phenotype_id in phenotype_ids:
                if phenotype_id not in phenotype_to_diseases:
                    phenotype_to_diseases[phenotype_id] = [disease_accession]
                else:
                    phenotype_to_diseases[phenotype_id].append(disease_accession)

        return disease_to_phenotypes, phenotype_to_diseases, phenotype_disease_frequencies

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
    return 0.5