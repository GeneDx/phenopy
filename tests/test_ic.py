import numpy as np
import os
import unittest

from phenopy import generate_alternate_ids
from phenopy.ic import calculate_information_content
from phenopy.network import annotate
from phenopy.network import load as load_network
from phenopy.d2p import load as load_d2p
from phenopy.util import export_phenotype_hpoa_with_no_parents


class ICTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # parent dir
        cls.parent_dir = os.path.dirname(os.path.realpath(__file__))

        # load and process the network
        cls.obo_file = os.path.join(cls.parent_dir, 'data/hp.obo')
        cls.hpo_network = load_network(cls.obo_file)
        cls.alt2prim = generate_alternate_ids(cls.hpo_network)

        # load phenotypes to genes associations
        cls.disease_to_phenotype_file = os.path.join(cls.parent_dir, 'data/phenotype.hpoa')
        cls.disease_records, cls.phenotype_to_diseases = load_d2p(cls.disease_to_phenotype_file, cls.hpo_network, cls.alt2prim)

        cls.num_diseases_annotated = len(cls.disease_records)
        cls.hpo_network = annotate(cls.hpo_network, cls.phenotype_to_diseases, cls.num_diseases_annotated, cls.alt2prim)

        cls.hpo_id = 'HP:0010863'
        cls.disease_to_phenotype_output_file = os.path.join(cls.parent_dir, 'data/phenotype.noparents.hpoa')

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(os.path.join(cls.parent_dir, 'data/phenotype.noparents.hpoa')):
            os.remove(os.path.join(cls.parent_dir, 'data/phenotype.noparents.hpoa'))

    def test_ic_d2p(self):
        """Calculate the information content of a phenotype"""
        self.assertAlmostEqual(self.hpo_network.nodes[self.hpo_id]['ic'], 5.69, 2)

    def test_ic_custom(self):
        """Calculate the information content of a phenotype when multiple annotations are present"""
        custom_annotation_file = os.path.join(self.parent_dir, 'data/test.score-long.txt')
        hpo_network = load_network(self.obo_file)
        hpo_network = annotate(hpo_network, self.phenotype_to_diseases, self.num_diseases_annotated, self.alt2prim,
                              annotations_file=custom_annotation_file)

        self.assertAlmostEqual(hpo_network.nodes[self.hpo_id]['ic'], 6.38, 1)

    def test_ic_d2p_no_parents(self):
        export_phenotype_hpoa_with_no_parents(self.disease_to_phenotype_file, self.disease_to_phenotype_output_file, self.hpo_network)
        self.assertTrue(os.path.exists(self.disease_to_phenotype_output_file))
