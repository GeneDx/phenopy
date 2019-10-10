import csv
import sys


def load(phenotype_annotations_file, logger=None):
    """Parse the hpoa file"""
    try:
        with open(phenotype_annotations_file, 'r') as tsv_fh:
            reader = csv.DictReader(filter(lambda line: line[0] != '#', tsv_fh), delimiter='\t')

            disease_to_phenotypes = dict()
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

                # add the phenotype to the disease in the disease_to_phenotypes dictionary
                if disease_accession not in disease_to_phenotypes:
                    disease_to_phenotypes[disease_accession] = [term_id]
                else:
                    disease_to_phenotypes[disease_accession].append(term_id)

        phenotype_to_diseases = dict()
        for disease_accession, phenotype_ids in disease_to_phenotypes.items():
            for phenotype_id in phenotype_ids:
                if phenotype_id not in phenotype_to_diseases:
                    phenotype_to_diseases[phenotype_id] = [disease_accession]
                else:
                    phenotype_to_diseases[phenotype_id].append(disease_accession)

        return disease_to_phenotypes, phenotype_to_diseases

    except (FileNotFoundError, PermissionError) as e:
        hpoa_file_error_msg = f'{phenotype_annotations_file} not found or incorrect permissions'
        if logger is not None:
            logger.critical(hpoa_file_error_msg)
        else:
            sys.stderr.write(hpoa_file_error_msg)
        exit(1)
