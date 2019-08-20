import os
import unittest

from phenosim.config import config
from phenosim.p2g import load as load_p2g
from phenosim.network import _load_hpo_network


class NetworkTestCase(unittest.TestCase):
    def setUp(self):
        config.set('hpo', 'data_directory', 'data')
        terms_to_genes, genes_to_terms, self.annotations_count = load_p2g('data/phenotypes_to_genes.txt')
        # choose an HPO id that is in the custom annotations file, so it should have different information content
        self.hpo_id = 'HP:0000545'

    def test_load_from_hp_obo(self):
        hpo_network = _load_hpo_network('data/hp.obo', 'data/phenotypes_to_genes.txt',
                                        self.annotations_count, custom_annotations_file=None)
        # clean up before asserting
        os.remove('data/hpo_network.pickle')
        # this is a cleaned version of the network, so it is not the same as test_obo.py
        self.assertEqual(len(hpo_network), 20)
        self.assertAlmostEqual(hpo_network.node[self.hpo_id]['ic'], 2.48, 2)

    def test_load_custom(self):
        hpo_network = _load_hpo_network('data/hp.obo', 'data/phenotypes_to_genes.txt',
                                        self.annotations_count, custom_annotations_file='data/test.score-product.txt')
        # clean up before asserting
        os.remove('data/hpo_network.pickle')
        self.assertAlmostEqual(hpo_network.node[self.hpo_id]['ic'], 2.14, 2)
