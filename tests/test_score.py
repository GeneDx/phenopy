import os
import unittest

from phenosim.obo import process
from phenosim.obo import load as load_obo
from phenosim.p2g import load as load_p2g
from phenosim.score import Scorer


class ScorerTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # parent dir
        parent_dir = os.path.dirname(os.path.realpath(__file__))

        # load phenotypes to genes associations
        pheno2genes_file = os.path.join(parent_dir, 'data/phenotypes_to_genes.txt')
        terms_to_genes, genes_to_terms, annotations_count = load_p2g(pheno2genes_file)

        # load and process the network
        obo_file = os.path.join(parent_dir, 'data/hp.obo')
        hpo_network = load_obo(obo_file)
        hpo_network = process(hpo_network, terms_to_genes, annotations_count)

        # create instance the scorer class
        cls.scorer = Scorer(hpo_network)

    def test_find_lca(self):
        # find the lowest common ancestor between glaucoma and myopia
        lca = self.scorer.find_lca('HP:0000501', 'HP:0000545')

        self.assertEqual(lca, 'HP:0012373')

        # make sure that when the root node is passed, it's returned as lca.
        root_lca = self.scorer.find_lca('HP:0000501', 'HP:0000001')

        self.assertEqual(root_lca, 'HP:0000001')

    def test_calculate_gamma(self):
        t1 = 'HP:0012372'
        t2 = 'HP:0012373'

        # term to itself should be distance 0
        gamma0 = self.scorer.calculate_gamma(t1, t1, 'HP:0000478')
        self.assertEqual(gamma0, 0)

        gamma2 = self.scorer.calculate_gamma(t1, t2, 'HP:0000478')
        self.assertEqual(gamma2, 2)
