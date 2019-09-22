import os
import unittest

from phenopy.obo import process
from phenopy.obo import load as load_obo
from phenopy.p2g import load as load_p2g
from phenopy.score import Scorer
from phenopy.util import remove_parents


class ScorerTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # parent dir
        parent_dir = os.path.dirname(os.path.realpath(__file__))

        # load phenotypes to genes associations
        pheno2genes_file = os.path.join(
            parent_dir, 'data/phenotype_to_genes.txt')
        terms_to_genes, genes_to_terms, annotations_count = load_p2g(
            pheno2genes_file)

        # load and process the network
        obo_file = os.path.join(parent_dir, 'data/hp.obo')
        hpo_network = load_obo(obo_file)
        hpo_network = process(hpo_network, terms_to_genes, annotations_count)

        # create instance the scorer class
        cls.scorer = Scorer(hpo_network)

    def test_find_lca(self):
        # find the lowest common ancestor between glaucoma and myopia
        lca = self.scorer.find_lca('HP:0001249', 'HP:0012434')
        self.assertEqual(lca, 'HP:0012759')

        # make sure that when the root node is passed, it's returned as lca.
        root_lca = self.scorer.find_lca('HP:0012759', 'HP:0000001')
        self.assertEqual(root_lca, 'HP:0000001')

        # LCA of parent-child is parent
        parent_lca = self.scorer.find_lca('HP:0012759', 'HP:0012758')
        self.assertEqual(parent_lca, 'HP:0012759')

        # LCA of parent-child is parent (inverse)
        parent_lca = self.scorer.find_lca('HP:0012758', 'HP:0012759')
        self.assertEqual(parent_lca, 'HP:0012759')

    def test_calculate_gamma(self):
        t1 = 'HP:0012758'
        t2 = 'HP:0012759'

        # term to itself should be distance 0
        gamma0 = self.scorer.calculate_gamma(t1, t1, t2)
        self.assertEqual(gamma0, 0)

        # term to a parent should be distance 1
        gamma1 = self.scorer.calculate_gamma(t1, t2, t2)
        self.assertEqual(gamma1, 1)

        # term to a neighbor should be 2
        gamma2 = self.scorer.calculate_gamma('HP:0000750', 'HP:0012434', t1)
        self.assertEqual(gamma2, 2)

    def test_calculate_beta(self):
        t1 = 'HP:0012758'
        t2 = 'HP:0012759'
        beta = self.scorer.calculate_beta(t1, t2)
        self.assertAlmostEqual(beta, 6.88, places=2)

    def test_score_hpo_pair_hrss(self):
        t1 = 'HP:0011351'
        t2 = 'HP:0012434'

        # score two terms
        score = self.scorer.score_hpo_pair_hrss((t1, t2))
        self.assertAlmostEqual(score, 0.2, places=2)

        # test that the cache is working
        score = self.scorer.score_hpo_pair_hrss((t1, t2))
        self.assertAlmostEqual(score, 0.2, places=2)

        # and test that the cache is working for inverse comparisons
        score = self.scorer.score_hpo_pair_hrss((t2, t1))
        self.assertAlmostEqual(score, 0.2, places=2)

    def test_score(self):
        terms_a = ['HP:0012433', 'HP:0012434']
        terms_b = ['HP:0001249', 'HP:0012758']

        # if no terms in one set, return 0.0
        score0 = self.scorer.score([], terms_b)
        self.assertEqual(score0, 0.0)

        # test BMA
        score_bma = self.scorer.score(terms_a, terms_b, agg_score='BMA')
        self.assertAlmostEqual(score_bma, 0.1462, places=4)

        score_max = self.scorer.score(terms_a, terms_b, agg_score='maximum')
        self.assertAlmostEqual(score_max, 0.25, places=4)

    def test_no_parents(self):
        terms_a = ['HP:0012433', 'HP:0000708']
        terms_b = ['HP:0001249', 'HP:0012758']

        self.assertEqual('HP:0012433', list(remove_parents(terms_a, self.scorer.hpo_network))[0])
        self.assertEqual(len(remove_parents(terms_b, self.scorer.hpo_network)), 2)

    def test_alt_term(self):
        terms_a = ['HP:0000715', 'HP:0012434']
        self.assertIn('HP:0000708', self.scorer.convert_alternate_ids(terms_a))
