from configparser import NoOptionError, NoSectionError
import os
import unittest

from phenopy.config import config
from phenopy.util import read_phenotype_groups
from phenopy.build_hpo import generate_annotated_hpo_network
from phenopy.likelihood import predict_likelihood_moldx


class LikelihoodTestCase(unittest.TestCase):
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

    def test_read_phenotype_groups_1000(self):
        hp_to_pg = read_phenotype_groups()
        self.assertEqual(hp_to_pg['HP:0012759']['k1000'], 62)

    def test_read_phenotype_groups_1500(self):
        hp_to_pg = read_phenotype_groups()
        self.assertEqual(hp_to_pg['HP:0012759']['k1500'], 739)

    def test_predict_likelihood(self):
        phenotypes = [
            ['HP:0012759', 'HP:0003011', 'HP:0011442'], 
            ['HP:0012759', 'HP:0003011'],
        ]
        probabilities = predict_likelihood_moldx(phenotypes, self.phenotype_groups, self.hpo_network, self.alt2prim)
        self.assertAlmostEqual(probabilities[0], 0.33, places=2)

    def test_predict_likelihood_phenotypes_only(self):
        phenotypes = [
            ['HP:0012759', 'HP:0003011', 'HP:0011442'], 
            ['HP:0012759', 'HP:0003011'],
        ]
        probabilities = predict_likelihood_moldx(phenotypes)
        self.assertAlmostEqual(probabilities[0], 0.33, places=2)

    def test_no_hpo_config_section(self):
        config.remove_section('hpo')
        phenotypes = [
            ['HP:0012759', 'HP:0003011', 'HP:0011442'], 
            ['HP:0012759', 'HP:0003011'],
        ]
        with self.assertRaises(NoSectionError):
            predict_likelihood_moldx(phenotypes)

    def test_no_hpo_config_option(self):
        config.remove_option('hpo', 'disease_to_phenotype_file')
        phenotypes = [
            ['HP:0012759', 'HP:0003011', 'HP:0011442'], 
            ['HP:0012759', 'HP:0003011'],
        ]
        with self.assertRaises(NoOptionError):
            predict_likelihood_moldx(phenotypes)   

    def test_bad_k(self):
        phenotypes = [
            ['HP:0012759', 'HP:0003011', 'HP:0011442'], 
            ['HP:0012759', 'HP:0003011'],
        ]
        with self.assertRaises(KeyError):
            predict_likelihood_moldx(
                phenotypes, 
                self.phenotype_groups,
                self.hpo_network,
                self.alt2prim, 
                k_phenotype_groups=500
                )