import os
import unittest

from phenosim.config import config
from phenosim.p2g import load as load_p2g
from phenosim.network import _load_hpo_network


class NetworkTestCase(unittest.TestCase):
    def setUp(cls):
        # parent dir
        cls.parent_dir = os.path.dirname(os.path.realpath(__file__))
        cls.obo_file = os.path.join(cls.parent_dir, 'data/hp.obo')
        cls.pheno2genes_file = os.path.join(cls.parent_dir, 'data/phenotypes_to_genes.txt')
        cls.hpo_network_file = os.path.join(cls.parent_dir, 'data/hpo_network.pickle')

        config.set('hpo', 'data_directory', os.path.join(cls.parent_dir, 'data'))
        terms_to_genes, genes_to_terms, cls.annotations_count = load_p2g(cls.pheno2genes_file)
        # choose an HPO id that is in the custom annotations file, so it should have different information content
        cls.hpo_id = 'HP:0000545'

    def tearDown(cls):
        os.remove(cls.hpo_network_file)

    def test_load_from_hp_obo(self):

        hpo_network = _load_hpo_network(self.obo_file, self.pheno2genes_file,
                                        self.annotations_count, custom_annotations_file=None,
                                        hpo_network_file=self.hpo_network_file)

        # this is a cleaned version of the network, so it is not the same as test_obo.py
        self.assertEqual(len(hpo_network), 20)
        self.assertAlmostEqual(hpo_network.node[self.hpo_id]['ic'], 2.48, 2)

    def test_load_custom(self):
        hpo_network = _load_hpo_network(self.obo_file, self.pheno2genes_file, self.annotations_count,
                                        custom_annotations_file=os.path.join(self.parent_dir,
                                                                             'data/test.score-product.txt'),
                                        hpo_network_file=self.hpo_network_file
                                        )

        self.assertAlmostEqual(hpo_network.node[self.hpo_id]['ic'], 2.14, 2)
