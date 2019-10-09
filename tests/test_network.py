import os
import unittest

import networkx as nx

from phenopy.config import config
from phenopy.d2p import load as load_d2p
from phenopy.network import _load_hpo_network


class NetworkTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # parent dir
        cls.parent_dir = os.path.dirname(os.path.realpath(__file__))
        cls.obo_file = os.path.join(cls.parent_dir, 'data/hp.obo')
        cls.phenotype_hpoa_file = os.path.join(cls.parent_dir, 'data/phenotype.hpoa')
        cls.hpo_network_file = os.path.join(cls.parent_dir, 'data/hpo_network.pickle')
        config.set('hpo', 'data_directory', os.path.join(cls.parent_dir, 'data'))
        cls.disease_to_phenotypes, cls.phenotype_to_diseases = load_d2p(cls.phenotype_hpoa_file)
        cls.hpo_id = 'HP:0010863'

    def tearDown(cls):
        if os.path.exists(cls.hpo_network_file):
            os.remove(cls.hpo_network_file)
        if os.path.exists(cls.hpo_network_file):
            os.remove(cls.hpo_network_file)

    def test_load_from_hp_obo(self):
        self.hpo_network = _load_hpo_network(self.obo_file, self.phenotype_to_diseases,
                                        len(self.disease_to_phenotypes), custom_annotations_file=None,
                                        hpo_network_file=self.hpo_network_file)

        # this is a cleaned version of the network, so it is not the same as test_obo.py
        self.assertEqual(len(self.hpo_network), 28)
        self.assertAlmostEqual(self.hpo_network.node[self.hpo_id]['ic'], 5.69, 2)

    def test_load_custom(self):
        self.hpo_network = _load_hpo_network(self.obo_file, self.phenotype_to_diseases, len(self.disease_to_phenotypes),
                                        custom_annotations_file=os.path.join(self.parent_dir,
                                                                             'data/test.score-product.txt'),
                                        hpo_network_file=self.hpo_network_file
                                        )

        self.assertAlmostEqual(self.hpo_network.node[self.hpo_id]['ic'], 6.38, 2)

    def test_restore(self):
        hpo_network1 = _load_hpo_network(self.obo_file, self.phenotype_to_diseases, len(self.disease_to_phenotypes),
                                        custom_annotations_file=os.path.join(self.parent_dir,
                                                                             'data/test.score-product.txt'),
                                        hpo_network_file=self.hpo_network_file
                                        )
        # this should get loaded from the pickled hpo_network_file
        hpo_network2 = _load_hpo_network(self.obo_file, self.phenotype_to_diseases, len(self.disease_to_phenotypes),
                                        custom_annotations_file=os.path.join(self.parent_dir,
                                                                             'data/test.score-product.txt'),
                                        hpo_network_file=self.hpo_network_file
                                        )
        self.assertTrue(nx.is_isomorphic(hpo_network1, hpo_network2))

    def test_phenotype_to_diseases(self):
        with self.assertRaises(ValueError) as se:
            self.hpo_network_file = _load_hpo_network(self.obo_file, self.phenotype_hpoa_file, len(self.disease_to_phenotypes),
                                                      custom_annotations_file=None, hpo_network_file=None)