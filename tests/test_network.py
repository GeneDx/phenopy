import os
import unittest

import networkx as nx

from phenopy.config import config
from phenopy.p2g import load as load_p2g
from phenopy.network import _load_hpo_network


class NetworkTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # parent dir
        cls.parent_dir = os.path.dirname(os.path.realpath(__file__))
        cls.obo_file = os.path.join(cls.parent_dir, 'data/hp.obo')
        cls.pheno2genes_file = os.path.join(cls.parent_dir, 'data/phenotype_to_genes.txt')
        cls.hpo_network_file = os.path.join(cls.parent_dir, 'data/hpo_network.pickle')
        config.set('hpo', 'data_directory', os.path.join(cls.parent_dir, 'data'))
        cls.terms_to_genes, cls.genes_to_terms, cls.num_genes_annotated = load_p2g(cls.pheno2genes_file)
        cls.hpo_id = 'HP:0010863'

    def tearDown(cls):
        if os.path.exists((cls.hpo_network_file)):
            os.remove(cls.hpo_network_file)

    def test_load_from_hp_obo(self):
        self.hpo_network = _load_hpo_network(self.obo_file, self.terms_to_genes,
                                        self.num_genes_annotated, custom_annotations_file=None,
                                        hpo_network_file=self.hpo_network_file)

        # this is a cleaned version of the network, so it is not the same as test_obo.py
        self.assertEqual(len(self.hpo_network), 16)
        self.assertAlmostEqual(self.hpo_network.node[self.hpo_id]['ic'], 6.12, 2)

    def test_load_custom(self):
        self.hpo_network = _load_hpo_network(self.obo_file, self.terms_to_genes, self.num_genes_annotated,
                                        custom_annotations_file=os.path.join(self.parent_dir,
                                                                             'data/test.score-product.txt'),
                                        hpo_network_file=self.hpo_network_file
                                        )

        self.assertAlmostEqual(self.hpo_network.node[self.hpo_id]['ic'], 7.09, 2)

    def test_restore(self):
        hpo_network1 = _load_hpo_network(self.obo_file, self.terms_to_genes, self.num_genes_annotated,
                                        custom_annotations_file=os.path.join(self.parent_dir,
                                                                             'data/test.score-product.txt'),
                                        hpo_network_file=self.hpo_network_file
                                        )
        # this should get loaded from the pickled hpo_network_file
        hpo_network2 = _load_hpo_network(self.obo_file, self.terms_to_genes, self.num_genes_annotated,
                                        custom_annotations_file=os.path.join(self.parent_dir,
                                                                             'data/test.score-product.txt'),
                                        hpo_network_file=self.hpo_network_file
                                        )
        self.assertTrue(nx.is_isomorphic(hpo_network1, hpo_network2))

    def test_terms_to_genes(self):
        with self.assertRaises(ValueError) as se:
            self.hpo_network_file = _load_hpo_network(self.obo_file, self.pheno2genes_file, self.num_genes_annotated,
                                                      custom_annotations_file=None, hpo_network_file=None)