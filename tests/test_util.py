import os
import unittest

from phenopy.util import parse, read_records_file

from phenopy.config import config
from phenopy.util import read_phenotype_groups, encode_phenotypes, parse_input, generate_annotated_hpo_network


class UtilTestCase(unittest.TestCase):
    @classmethod
    def setUp(cls):
        # parent dir
        cls.parent_dir = os.path.dirname(os.path.realpath(__file__))
        
        if 'hpo' not in config.sections():
            config.add_section('hpo')
        
        config.set('hpo', 'obo_file', os.path.join(cls.parent_dir, 'data/hp.obo'))
        config.set('hpo', 'disease_to_phenotype_file', os.path.join(cls.parent_dir, 'data/phenotype.hpoa'))

        cls.obo_file = config.get('hpo', 'obo_file')
        cls.disease_to_phenotype_file = config.get('hpo', 'disease_to_phenotype_file')

        cls.hpo_network, cls.alt2prim, cls.disease_records = generate_annotated_hpo_network(
            cls.obo_file,
            cls.disease_to_phenotype_file,
            )
        cls.phenotype_groups = read_phenotype_groups()


    def test_read_records_file(self):
        with self.assertRaises(SystemExit) as se:
            read_records_file('notafilepath/notafile')

        syserr = se.exception
        self.assertEqual(syserr.code, 1)
        records_truth = [
            {
                'sample': '118200',
                'age': 9.0,
                'gender': 'Female',
                'terms': 'HP:0001263|HP:0001251|HP:0001290|HP:0004322'.split('|')
             },
            {
                'sample': '118210',
                'age': 4.0,
                'gender': None,
                'terms': 'HP:0001249|HP:0001263|HP:0001290'.split('|')
            },
            {
                'sample': '118211',
                'age': None,
                'gender': None,
                'terms': 'HP:0001249|HP:0001263|HP:0001290'.split('|')
            }
        ]
        records_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/test.score-short.txt')
        records = read_records_file(records_path, no_parents=False)
        self.assertEqual(records,records_truth)


    def test_parse(self):
        string = 'age=13;sex=Male'
        self.assertEqual(parse(string, what='sex'), 'Male')
        self.assertEqual(parse(string, what='age'), 13.0)

        string = 'age=13.64;sex=male'
        self.assertEqual(parse(string, what='sex'), 'Male')
        self.assertEqual(parse(string, what='age'), 13.6)

        string = 'age=12.9;sex=female'
        self.assertEqual(parse(string, what='sex'), 'Female')
        self.assertEqual(parse(string, what='age'), 12.9)

        string = 'sex=Female'
        self.assertEqual(parse(string, what='sex'), 'Female')

        string = 'sex=FEMALE'
        self.assertEqual(parse(string, what='sex'), 'Female')

        string = 'sex=F'
        self.assertEqual(parse(string, what='sex'), 'Female')

        string = 'age=1'
        self.assertEqual(parse(string, what='age'), 1.0)

        string = '.'
        self.assertEqual(parse(string, what='age'), None)

        string = '. '
        self.assertEqual(parse(string, what='age'), None)

        string = ' . '
        self.assertEqual(parse(string, what='age'), None)

        string = '13?'
        self.assertEqual(parse(string, what='age'), None)

        string = 'sex=NA'
        self.assertEqual(parse(string, what='sex'), None)

        string = 'sex=Unknown'
        self.assertEqual(parse(string, what='sex'), None)


    def test_encode_phenotypes_file(self):
        input_file = os.path.join(self.parent_dir, "data/test.score-short.txt")
        records = parse_input(input_file, self.hpo_network, self.alt2prim)
        encoded_phenotypes = encode_phenotypes(
            [record["terms"] for record in records],
            self.phenotype_groups,
            self.hpo_network,
            self.alt2prim
        )
        self.assertEqual(sum(encoded_phenotypes[0]), 4)

    
    def test_encode_1d_phenotypes(self):
        phenotypes = ['HP:0012759', 'HP:0003011', 'HP:0011442']
        encoded_phenotypes = encode_phenotypes(phenotypes, self.phenotype_groups, self.hpo_network, self.alt2prim, k=1000)
        self.assertEqual(sum(encoded_phenotypes), 3)


    def test_encode_2d_phenotypes(self):
        phenotypes = [
            ['HP:0012759', 'HP:0003011', 'HP:0011442'], 
            ['HP:0012759', 'HP:0003011'],
        ]
        encoded_phenotypes = encode_phenotypes(phenotypes, self.phenotype_groups, self.hpo_network, self.alt2prim, k=1000)
        self.assertEqual(sum(encoded_phenotypes[1]), 2)