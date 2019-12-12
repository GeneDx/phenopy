import itertools
import numpy as np
import os
import pandas as pd
import unittest

from phenopy.d2p import load as load_d2p
from phenopy.network import annotate
from phenopy.network import load as load_network
from phenopy.score import Scorer
from phenopy.util import remove_parents, parse_input, generate_alternate_ids, half_product
from phenopy.weights import calculate_age_weights


class ScorerTestCase(unittest.TestCase):
    @classmethod
    def setUp(cls):
        # parent dir
        cls.parent_dir = os.path.dirname(os.path.realpath(__file__))

        # load and process the network
        cls.obo_file = os.path.join(cls.parent_dir, 'data/hp.obo')
        cls.hpo_network = load_network(cls.obo_file)
        cls.alt2prim = generate_alternate_ids(cls.hpo_network)
        cls.ages_distribution_file = os.path.join(cls.parent_dir, 'data/phenotype_age.tsv')

        # load phenotypes to genes associations
        cls.disease_to_phenotype_file = os.path.join(cls.parent_dir, 'data/phenotype.hpoa')
        cls.disease_records, cls.phenotype_to_diseases = load_d2p(cls.disease_to_phenotype_file, cls.hpo_network, cls.alt2prim)

        cls.num_diseases_annotated = len(cls.disease_records)
        cls.hpo_network = annotate(cls.hpo_network, cls.phenotype_to_diseases, cls.num_diseases_annotated, cls.alt2prim)

        # create instance the scorer class
        cls.scorer = Scorer(cls.hpo_network, min_score_mask=None)

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

        # LCA of self-self is self
        parent_lca = self.scorer.find_lca('HP:0012759', 'HP:0012759')
        self.assertEqual(parent_lca, 'HP:0012759')

        # LCA of grandparent-child is grandparent
        parent_lca = self.scorer.find_lca('HP:0012759', 'HP:0000750')
        self.assertEqual(parent_lca, 'HP:0012759')

    def test_calculate_gamma(self):
        t1 = 'HP:0012758'
        t2 = 'HP:0012759'

        # term to itself should be distance 0
        gamma0 = self.scorer.calculate_gamma(t1, t1, t2)
        self.assertEqual(gamma0, 0)

        # term to a parent should be distance 1
        gamma1a = self.scorer.calculate_gamma(t1, t2, t2)
        self.assertEqual(gamma1a, 1)

        gamma1b = self.scorer.calculate_gamma(t2, t1, t2)
        self.assertEqual(gamma1b, 1)

        # term to a neighbor should be 2
        gamma2 = self.scorer.calculate_gamma('HP:0000750', 'HP:0012434', t1)
        self.assertEqual(gamma2, 2)

    def test_calculate_beta(self):
        t1 = 'HP:0001344'
        t2 = 'HP:0012759'
        beta = self.scorer.calculate_beta(t1, t2)
        self.assertAlmostEqual(beta, 3.51, places=2)

    def test_score_hpo_pair_hrss(self):
        t1 = 'HP:0011351'
        t2 = 'HP:0012434'

        # score two terms
        score = self.scorer.score_hpo_pair_hrss(t1, t2)
        self.assertAlmostEqual(score, 0.2, places=2)

        # test that the cache is working
        score = self.scorer.score_hpo_pair_hrss(t1, t2)
        self.assertAlmostEqual(score, 0.2, places=2)

        # and test that the cache is working for inverse comparisons
        score = self.scorer.score_hpo_pair_hrss(t2, t1)
        self.assertAlmostEqual(score, 0.2, places=2)

    def test_score(self):
        record_a = {'record_id': 'sample_1',
                    'terms': ['HP:0012433', 'HP:0012434'],
                    'weights': {}
                    }
        record_b = {'record_id': 'sample_2',
                    'terms': [],
                    'weights': {}
                    }

        # if no terms in one set, return 0.0
        score0 = self.scorer.score(record_a, record_b)
        self.assertEqual(score0, 0.0)
        record_b['terms'] = ['HP:0001249', 'HP:0012758']

        # test BMA
        score_bma = self.scorer.score(record_a, record_b)
        self.assertAlmostEqual(score_bma, 0.207, places=2)
        self.scorer.summarization_method = 'maximum'
        score_max = self.scorer.score(record_a, record_b)
        self.assertAlmostEqual(score_max, 0.25, places=4)

        # test wrong method
        self.scorer.summarization_method = 'not_a_method'
        with self.assertRaises(ValueError):
            self.scorer.score(record_a, record_b)

        # test BMWA with age weights
        record_a.update({
            'terms': ['HP:0001251', 'HP:0001263', 'HP:0001290', 'HP:0004322', 'HP:0012433'], # ATAX, DD,  HYP, SS, AbnSocBeh
            'weights': {'age': [0.67, 1., 1., 0.4, 0.4]},
        })
        record_b.update({
            'terms': ['HP:0001249', 'HP:0001263', 'HP:0001290'],  # ID,  DD, HYP
            'weights': {'age': [1., 1., 1.]},
        })

        self.scorer.summarization_method = 'BMWA'
        self.scorer.min_score_mask = 0.05
        score_bmwa = self.scorer.score(record_a, record_b)
        self.assertAlmostEqual(score_bmwa, 0.6239, places=4)

        record_a.update({
            'terms': ['HP:0001251', 'HP:0001263', 'HP:0001290', 'HP:0004322'], # ATAX, DD,  HYP, SS, AbnSocBeh
            'weights': {'age': [0.67, 1., 1., 0.4]},
        })
        record_b.update({
            'terms': ['HP:0001263', 'HP:0001249', 'HP:0001290'],  # ID,  DD, HYP
            'weights': {'age': [1., 1., 0.5]},
        })

        scorer = self.scorer
        scorer.summarization_method = 'BMWA'

        # test with two weights
        score_bwma_both_weights = scorer.score(record_a, record_b)
        self.assertAlmostEqual(score_bwma_both_weights, 0.6901, 4)

        # test with one weight array
        scorer.min_score_mask = None
        record_a['weights'].pop('age', None)
        score_bwma_one_weights = scorer.score(record_a, record_b)
        self.assertAlmostEqual(score_bwma_one_weights, 0.6026, 4)

    def test_score_records(self,):
        query_name = 'SAMPLE'
        query_terms = [
            'HP:0000750',
            'HP:0010863',
        ]
        input_records = [{
            'record_id': query_name,
            'terms': query_terms,
            'weights': {}
        }]
        score_records = self.disease_records

        # using weights
        results = self.scorer.score_records(
            input_records,
            score_records,
            itertools.product(range(len(input_records)), range(len(score_records))),
            thread_index=0,
            threads=1,
        )
        self.assertEqual(1184, len(results))
        self.assertAlmostEqual(float(results[0][2]), 0.1325, 2)

        # without weights -- defaults to best match average (BMA)
        [record['weights'].pop('disease_frequency') for record in score_records]
        results = self.scorer.score_records(
            input_records,
            score_records,
            itertools.product(range(len(input_records)), range(len(score_records))),
            thread_index=0,
            threads=1,
        )
        self.assertEqual(1184, len(results))
        self.assertAlmostEqual(float(results[0][2]), 0.1387, 2)

    def test_no_parents(self):
        terms_a = ['HP:0012433', 'HP:0000708']
        terms_b = ['HP:0001249', 'HP:0012758']

        self.assertEqual('HP:0012433', list(remove_parents(terms_a, self.scorer.hpo_network))[0])
        self.assertEqual(len(remove_parents(terms_b, self.scorer.hpo_network)), 2)

    def test_score_self(self):
        # read in records
        records = parse_input(os.path.join(self.parent_dir, 'data/test.score-long.txt'),
                              self.hpo_network,
                              self.alt2prim
                              )

        # limit to records with HPO terms since many test cases don't have the sub-graph terms from tests/data/hp.obo
        input_records = [x for x in records if x['record_id'] in ['213200', '302801']]

        results = self.scorer.score_records(
            input_records,
            input_records,
            half_product(len(input_records), len(input_records))
        )
        self.assertEqual(len(results), 3)

        # test the score of '213200' - '302801'
        self.assertAlmostEqual(float(results[1][2]), 0.3758, 2)

    def test_bmwa(self):
        # test best match weighted average
        # load and annotate the network

        terms_a = ['HP:0001251', 'HP:0001263', 'HP:0001290', 'HP:0004322']  # ATAX, DD, HYP, SS

        terms_b = ['HP:0001263', 'HP:0001249', 'HP:0001290']  # DD, ID, HYP
        weights_a = {'age': [0.67, 1., 1., 0.4]}
        weights_b = {'age': [1., 1., 1.]}

        df = pd.DataFrame(
            [[4.22595743e-02, 3.92122308e-02, 3.04851573e-04],
             [1.07473687e-01, 5.05101655e-01, 3.78305515e-04],
             [3.69780479e-04, 3.78305515e-04, 4.64651944e-01],
             [4.17139800e-04, 4.12232546e-04, 3.67984322e-04]],
            index=pd.Index(terms_a, name='a'),
            columns=pd.MultiIndex.from_arrays([['score'] * len(terms_b), terms_b],
                                              names=[None, 'b'])
        )

        score_bmwa = self.scorer.best_match_weighted_average(df, weights_a, weights_b)

        self.assertAlmostEqual(score_bmwa, 0.3419, 4)

        # set all weights to 1.0
        weights_a = {'age': [1.] * len(terms_a)}
        score_bmwa = self.scorer.best_match_weighted_average(df, weights_a, weights_b)
        self.assertAlmostEqual(score_bmwa, 0.2985, 4)

        # set all weights to 0.0, result should be the same as all weights being 1.0
        weights_a = {'age': [1.] * len(terms_a)}
        weights_b = {'age': [1.] * len(terms_b)}
        self.min_score_mask = None
        score_bmwa = self.scorer.best_match_weighted_average(df, weights_a, weights_b)
        self.assertAlmostEqual(score_bmwa, 0.2985, 4)

        # Test weight adjustment masking
        # make pairwise scores matrix

        # Patient A age = 9.0 years
        # Patient B age = 4.0 years

        terms_a = ['HP:0001251', 'HP:0001249', 'HP:0001263', 'HP:0001290', 'HP:0004322']  # ATAX, ID, DD, HYP, SS
        terms_b = ['HP:0001263', 'HP:0001249', 'HP:0001290']  # DD, ID, HYP

        df = pd.DataFrame(
            [[4.22595743e-02, 3.92122308e-02, 3.04851573e-04],
             [1.07473687e-01, 5.05101655e-01, 3.78305515e-04],
             [1.07473687e-01, 5.05101655e-01, 3.78305515e-04],
             [3.69780479e-04, 3.78305515e-04, 4.64651944e-01],
             [4.17139800e-04, 4.12232546e-04, 3.67984322e-04]],
            index=pd.Index(terms_a, name='a'),
            columns=pd.MultiIndex.from_arrays([['score'] * len(terms_b), terms_b],
                                              names=[None, 'b'])
        )

        # calculate weights based on patients age

        weights_a = {'age': [0.67, .4, 1., 1., 0.4]}  # patient_b is too young to have ataxia, ID and short stature
        weights_b = {'age': [1., 1., 1.]}

        # compute pairwise best match weighted average
        self.scorer.min_score_mask = None
        score_bmwa = self.scorer.best_match_weighted_average(df, weights_a, weights_b)

        self.assertAlmostEqual(score_bmwa, 0.352, 4)

        # because both patients were described to have ID, but only patient a has ataxia and ss
        # we mask good phenotype matches from being weighted down by default
        # we expect to get a better similarity score
        self.scorer.min_score_mask = 0.05
        score_bmwa = self.scorer.best_match_weighted_average(df, weights_a, weights_b)

        self.assertAlmostEqual(score_bmwa, 0.365, 4)

    def test_age_weight(self):
        # Test age based weight distribution and best_match_weighted_average calculation

        terms_a = ['HP:0001251', 'HP:0001263', 'HP:0001290', 'HP:0004322']  # ATAX, DD, HYP, SS
        terms_b = ['HP:0001263', 'HP:0001249', 'HP:0001290']  # DD, ID, HYP

        self.hpo_network = annotate(self.hpo_network, self.phenotype_to_diseases, self.num_diseases_annotated,
                                    self.alt2prim, ages_distribution_file=self.ages_distribution_file)

        age_a = 9.0
        age_b = 4.0

        # calculate weights based on patients age
        weights_a = {'age': calculate_age_weights(terms_a, age_b, self.hpo_network)}
        weights_b = {'age': calculate_age_weights(terms_b, age_a, self.hpo_network)}

        # make pairwise scores matrix
        df = pd.DataFrame(
            [[4.22595743e-02,   3.92122308e-02, 3.04851573e-04],
             [1.07473687e-01,   5.05101655e-01, 3.78305515e-04],
             [3.69780479e-04,   3.78305515e-04, 4.64651944e-01],
             [4.17139800e-04,   4.12232546e-04, 3.67984322e-04]],
            index=pd.Index(terms_a, name='a'),
            columns=pd.MultiIndex.from_arrays([['score'] * len(terms_b), terms_b],
                                              names=[None, 'b'])
        )
        # compute pairwise best match weighted average
        score_bmwa = self.scorer.best_match_weighted_average(df, weights_a, weights_b)

        self.assertAlmostEqual(score_bmwa, 0.3741, 4)

        # set all weights to 1.0, result should be the same as BMA without weights
        weights_a = {'disease_frequency': [1.] * len(terms_a)}
        weights_b = {'disease_frequency': [1.] * len(terms_b)}
        score_bmwa = self.scorer.best_match_weighted_average(df, weights_a, weights_b)

        self.assertAlmostEqual(score_bmwa, 0.2985, 4)

        # test term not in network
        terms_a = ['HP:Not_a_term']
        weights_a = calculate_age_weights(terms_a, age_b, self.hpo_network)
        self.assertEqual(weights_a, [1.0])

        # term in network no age
        terms_a = ['HP:0000001']
        weights_a = calculate_age_weights(terms_a, age_b, self.hpo_network)
        self.assertEqual(weights_a, [1.0])

    def test_score_pairs_age(self):
        # Test reading in records files and calculating pairwise scores
        # read in records
        self.hpo_network = annotate(self.hpo_network, self.phenotype_to_diseases, self.num_diseases_annotated,
                                    self.alt2prim, ages_distribution_file=self.ages_distribution_file)

        records = parse_input(os.path.join(self.parent_dir, 'data/test.score-short.txt'), self.hpo_network,
                              self.alt2prim)

        # create instance the scorer class
        scorer = Scorer(self.hpo_network, summarization_method='BMWA', min_score_mask=None)

        # select which patients to test in pairwise best_match_weighted_average
        input_records = [x for x in records if x['record_id'] in ['118200', '118210']]

        results = scorer.score_records(
            input_records,
            input_records,
            [(0, 1), ],
        )
        self.assertEqual(len(results), 1)

        # the right answer =
        answer = np.average([0.166, 1.0, 1.0, 0.125, 0.25, 1.0, 1.0], weights=[0.481, 1.0, 1.0, 0.0446, 1.0, 1.0, 1.0])

        self.assertAlmostEqual(float(results[0][2]), answer, 2)

        # Test identical records for which one age exist and one doesn't
        input_records = [x for x in records if x['record_id'] in ['118210', '118211']]

        results = scorer.score_records(
            input_records,
            input_records,
            [(0, 1), ],
        )
        self.assertEqual(len(results), 1)

        self.assertAlmostEqual(float(results[0][2]), 1.0, 1)

    def test_alpha_zero(self):
        """the root term should contain all diseases therefore the IC should be zero"""

        root_term_ic = self.hpo_network.nodes['HP:0000118']['ic']
        self.assertEqual(0.0,root_term_ic)

    def test_leaves_diff_branches_score_zero(self):
        """two leaves in different branches
        two leaves therfore beta is zero
        different branches therefore alpha is zero
        define I = (0.0 / (0.0 + 0.0)) as zero and not nan"""

        # generalized hypotonia
        term_a = 'HP:0001290'

        # moderate receptive langage delay
        term_b = 'HP:0011351'

        score_two_leaves_diff_branches = self.scorer.score_hpo_pair_hrss(term_a,term_b)
        self.assertEqual(0.0,score_two_leaves_diff_branches)
