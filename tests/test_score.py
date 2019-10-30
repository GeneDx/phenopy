import itertools
import numpy as np
import os
import pandas as pd
import unittest

from phenopy import parse_input
from phenopy.obo import process
from phenopy.obo import load as load_obo
from phenopy.d2p import load as load_d2p
from phenopy.score import Scorer
from phenopy.util import remove_parents
from phenopy.weights import get_truncated_normal


class ScorerTestCase(unittest.TestCase):
    @classmethod
    def setUp(cls):
        # parent dir
        cls.parent_dir = os.path.dirname(os.path.realpath(__file__))
        cls.hpo_network_file = os.path.join(cls.parent_dir, 'data/hpo_network.pickle')
        cls.records_file = os.path.join(cls.parent_dir, 'data/test.score-product.txt')

        # load phenotypes to genes associations
        phenotype_hpoa_file = os.path.join(
            cls.parent_dir, 'data/phenotype.hpoa')
        cls.disease_to_phenotypes, cls.phenotype_to_diseases, cls.phenotype_disease_frequencies = load_d2p(
            phenotype_hpoa_file)
        cls.num_diseases_annotated = len(cls.disease_to_phenotypes)

        # load and process the network
        obo_file = os.path.join(cls.parent_dir, 'data/hp.obo')
        hpo_network = load_obo(obo_file)
        cls.hpo_network = process(hpo_network, cls.phenotype_to_diseases, len(cls.disease_to_phenotypes),
                                  phenotype_disease_frequencies=cls.phenotype_disease_frequencies)

        # create instance the scorer class
        cls.scorer = Scorer(hpo_network, min_score_mask=None)

    def tearDown(cls):
        if os.path.exists(cls.hpo_network_file):
            os.remove(cls.hpo_network_file)

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
        score_bma = self.scorer.score(terms_a, terms_b)
        self.assertAlmostEqual(score_bma, 0.207, places=2)
        self.scorer.summarization_method = 'maximum'
        score_max = self.scorer.score(terms_a, terms_b)
        self.assertAlmostEqual(score_max, 0.25, places=4)

        # test wrong method
        self.scorer.summarization_method = 'not_a_method'
        score_not_method = self.scorer.score(terms_a, terms_b)
        self.assertAlmostEqual(score_not_method, 0.0, places=4)

        # test BMWA with age weights
        terms_a = ['HP:0001251', 'HP:0001263', 'HP:0001290', 'HP:0004322', 'HP:0012433'] # ATAX, DD,  HYP, SS, AbnSocBeh
        terms_b = ['HP:0001249', 'HP:0001263', 'HP:0001290']  # ID,  DD, HYP
        weights_a = [0.67, 1., 1., 0.4, 0.4]
        weights_b = [1., 1., 1.]
        weights = [weights_a, weights_b]
        self.scorer.summarization_method = 'BMWA'
        self.scorer.min_score_mask = 0.05
        score_bmwa = self.scorer.score(terms_a, terms_b, weights=weights)
        self.assertAlmostEqual(score_bmwa, 0.5927, places=4)

        terms_a = ['HP:0001251', 'HP:0001263', 'HP:0001290', 'HP:0004322']  # ATAX, DD, HYP, SS
        terms_b = ['HP:0001263', 'HP:0001249', 'HP:0001290']  # DD, ID, HYP
        weights_a = [0.67, 1., 1., 0.4]
        weights_b = [1., 1., 0.5]

        scorer = self.scorer
        scorer.summarization_method = 'BMWA'

        terms_a = scorer.convert_alternate_ids(terms_a)
        terms_b = scorer.convert_alternate_ids(terms_b)

        terms_a = scorer.filter_and_sort_hpo_ids(terms_a)
        terms_b = scorer.filter_and_sort_hpo_ids(terms_b)

        # test with two weights
        score_bwma_both_weights = scorer.score(terms_a, terms_b, [weights_a, weights_b])
        self.assertEqual(score_bwma_both_weights, 0.6488)

        # test with one weight array
        scorer.min_score_mask = None
        score_bwma_one_weights = scorer.score(terms_a, terms_b, [weights_b])
        self.assertEqual(score_bwma_one_weights, 0.6218)

        # make sure phenopy exits if length of weights is not 1 or 2.
        with self.assertRaises(SystemExit):
            score_bwma_zero_weights = scorer.score(terms_a, terms_b, [])
        with self.assertRaises(SystemExit):
            score_bwma_three_weights = scorer.score(terms_a, terms_b, [[0.1], [0.2], [0.3]])

    def test_score_records(self,):
        query_name = 'SAMPLE'
        query_terms = [
            'HP:0000750',
            'HP:0010863',
        ]
        input_records = [{
            'record_id': query_name,
            'terms': query_terms,
        }]
        score_records = self.disease_to_phenotypes
        #
        for hpo_id in query_terms:
            self.hpo_network.node[hpo_id]['weights']['disease_frequency'][query_name] = 1.0
        for record in score_records:
            record['terms'] = self.scorer.convert_alternate_ids(record['terms'])
            record['terms'] = self.scorer.filter_and_sort_hpo_ids(record['terms'])

        # using weights
        results = self.scorer.score_records(
            input_records,
            score_records,
            itertools.product(range(len(input_records)), range(len(score_records))),
            thread_index=0,
            threads=1,
            use_weights=True,
        )
        self.assertEqual(1184, len(results))
        self.assertAlmostEqual(float(results[0][2]), 0.1387, 2)

        # without weights
        results = self.scorer.score_records(
            input_records,
            score_records,
            itertools.product(range(len(input_records)), range(len(score_records))),
            thread_index=0,
            threads=1,
            use_weights=False,
        )
        self.assertEqual(1184, len(results))
        self.assertAlmostEqual(float(results[0][2]), 0.1387, 2)

    def test_no_parents(self):
        terms_a = ['HP:0012433', 'HP:0000708']
        terms_b = ['HP:0001249', 'HP:0012758']

        self.assertEqual('HP:0012433', list(remove_parents(terms_a, self.scorer.hpo_network))[0])
        self.assertEqual(len(remove_parents(terms_b, self.scorer.hpo_network)), 2)

    def test_alt_term(self):
        terms_a = ['HP:0000715', 'HP:0012434']
        self.assertIn('HP:0000708', self.scorer.convert_alternate_ids(terms_a))

    def test_score_self(self):
        # read in records
        records = parse_input(os.path.join(self.parent_dir, 'data/test.score-product.txt'))

        for record in records:
            # prune parent terms
            record['terms'] = remove_parents(record['terms'], self.hpo_network)

            # convert and filter the query hpo ids
            record['terms'] = self.scorer.convert_alternate_ids(record['terms'])
            record['terms'] = self.scorer.filter_and_sort_hpo_ids(record['terms'])

        # limit to records with HPO terms since many test cases don't have the sub-graph terms from tests/data/hp.obo
        input_records = [x for x in records if x['record_id'] in ['213200', '302801']]

        results = self.scorer.score_records(
            input_records,
            input_records,
            itertools.combinations(range(len(input_records)), 2),
        )
        self.assertEqual(len(results), 1)

        # test the score of '213200' - '302801'
        self.assertAlmostEqual(float(results[0][2]), 0.415, 2)

    def test_bmwa(self):
        # test best match weighted average
        # load and process the network

        terms_a = ['HP:0001251', 'HP:0001263', 'HP:0001290', 'HP:0004322']  # ATAX, DD, HYP, SS

        terms_b = ['HP:0001263', 'HP:0001249', 'HP:0001290']  # DD, ID, HYP
        weights_a = [0.67, 1., 1., 0.4]
        weights_b = [1., 1., 1.]

        df = pd.DataFrame(
            [[4.22595743e-02, 3.92122308e-02, 3.04851573e-04],
             [1.07473687e-01, 5.05101655e-01, 3.78305515e-04],
             [3.69780479e-04, 3.78305515e-04, 4.64651944e-01],
             [4.17139800e-04, 4.12232546e-04, 3.67984322e-04]],
            index=pd.Index(terms_a, name='a'),
            columns=pd.MultiIndex.from_arrays([['score'] * len(terms_b), terms_b],
                                              names=[None, 'b'])
        )

        score_bmwa = self.scorer.bmwa(df, weights_a, weights_b)

        self.assertEqual(score_bmwa, 0.3419)

        # set all weights to 1.0
        weights_a = np.ones(len(weights_a))
        score_bmwa = self.scorer.bmwa(df, weights_a, weights_b)
        self.assertEqual(score_bmwa, 0.2985)

        # set all weights to 0.0, result should be the same as all weights being 1.0
        weights_a = np.zeros(len(weights_a))
        weights_b = np.zeros(len(weights_b))
        self.min_score_mask = None
        score_bmwa = self.scorer.bmwa(df, weights_a, weights_b)
        self.assertEqual(score_bmwa, 0.2985)

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

        weights_a = [0.67, .4, 1., 1., 0.4]  # patient_b is too young to have ataxia, ID and short stature
        weights_b = [1., 1., 1.]

        # compute pairwise best match weighted average
        self.scorer.min_score_mask = None
        score_bmwa = self.scorer.bmwa(df, weights_a, weights_b)

        self.assertEqual(score_bmwa, 0.352)

        # because both patients were described to have ID, but only patient a has ataxia and ss
        # we mask good phenotype matches from being weighted down by default
        # we expect to get a better similarity score
        self.scorer.min_score_mask = 0.05
        score_bmwa = self.scorer.bmwa(df, weights_a, weights_b)

        self.assertEqual(score_bmwa, 0.365)

    def test_age_weight(self):
        # Test age based weight distribution and bmwa calculation

        terms_a = ['HP:0001251', 'HP:0001263', 'HP:0001290', 'HP:0004322']  # ATAX, DD, HYP, SS
        terms_b = ['HP:0001263', 'HP:0001249', 'HP:0001290']  # DD, ID, HYP

        # model age using truncated normal
        ages = pd.DataFrame([
                {'hpid': 'HP:0001251', 'age_dist': get_truncated_normal(6.0, 3.0, 0.0, 6.0)},
                {'hpid': 'HP:0001263', 'age_dist': get_truncated_normal(1.0, 1.0, 0.0, 1.0)},
                {'hpid': 'HP:0001290', 'age_dist': get_truncated_normal(1.0, 1.0, 0.0, 1.0)},
                {'hpid': 'HP:0004322', 'age_dist': get_truncated_normal(10.0, 3.0, 0.0, 10.0)},
                {'hpid': 'HP:0001249', 'age_dist': get_truncated_normal(6.0, 3.0, 0.0, 6.0)},
                ]).set_index('hpid')

        self.hpo_network = process(self.hpo_network, self.phenotype_to_diseases, self.num_diseases_annotated, ages=ages)

        age_a = 9.0
        age_b = 4.0

        # calculate weights based on patients age
        weights_a = self.scorer.calculate_age_weights(terms_a, age_b)
        weights_b = self.scorer.calculate_age_weights(terms_b, age_a)

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
        score_bmwa = self.scorer.bmwa(df, weights_a, weights_b)

        self.assertEqual(score_bmwa, 0.3741)

        # set all weights to 1.0, result should be the same as BMA without weights
        weights_a = [1.] * len(terms_a)
        weights_b = [1.] * len(terms_b)
        score_bmwa = self.scorer.bmwa(df, weights_a, weights_b)

        self.assertEqual(score_bmwa, 0.2985)

        # test term not in network
        terms_a = ['HP:Not_a_term']
        weights_a = self.scorer.calculate_age_weights(terms_a, age_b)
        self.assertEqual(weights_a, [1.0])

        # term in network no age

        terms_a = ['HP:0000001']
        weights_a = self.scorer.calculate_age_weights(terms_a, age_b)
        self.assertEqual(weights_a, [1.0])

    def test_score_pairs_age(self):
        # Test reading in records files and calculating pairwise scores
        # read in records
        records = parse_input(os.path.join(self.parent_dir, 'data/test.score-product-age.txt'))

        for record in records:
            # prune parent terms
            record['terms'] = remove_parents(record['terms'], self.hpo_network)

            # convert and filter the query hpo ids
            record['terms'] = self.scorer.convert_alternate_ids(record['terms'])
            record['terms'] = self.scorer.filter_and_sort_hpo_ids(record['terms'])

        # model age using truncated normal
        ages = pd.DataFrame([
            {'hpid': 'HP:0001251', 'age_dist': get_truncated_normal(6.0, 3.0, 0.0, 6.0)},
            {'hpid': 'HP:0001263', 'age_dist': get_truncated_normal(1.0, 1.0, 0.0, 1.0)},
            {'hpid': 'HP:0001290', 'age_dist': get_truncated_normal(1.0, 1.0, 0.0, 1.0)},
            {'hpid': 'HP:0004322', 'age_dist': get_truncated_normal(10.0, 3.0, 0.0, 10.0)},
            {'hpid': 'HP:0001249', 'age_dist': get_truncated_normal(6.0, 3.0, 0.0, 6.0)},
        ]).set_index('hpid')

        self.hpo_network = process(self.hpo_network, self.phenotype_to_diseases, self.num_diseases_annotated, ages=ages)

        # create instance the scorer class
        scorer = Scorer(self.hpo_network, summarization_method='BMWA', min_score_mask=None)

        # select which patients to test in pairwise bmwa
        input_records = [x for x in records if x['record_id'] in ['118200', '118210']]

        # sort terms in records
        sorted_records = []
        for record in input_records:
            record['terms'] = sorted(record['terms'])
            sorted_records.append(record)

        results = scorer.score_records(
            sorted_records,
            sorted_records,
            itertools.combinations(range(len(sorted_records)), 2),
        )
        self.assertEqual(len(results), 4)

        # the right answer =
        answer = np.average([0.166, 1.0, 1.0, 0.125, 0.25, 1.0, 1.0], weights=[0.481, 1.0, 1.0, 0.0446, 1.0, 1.0, 1.0])

        self.assertAlmostEqual(float(results[1][2]), answer, 4)

        # select  patients as before but without pre-sorting hpids
        input_records = [x for x in records if x['record_id'] in ['118200', '118210']]

        results = scorer.score_records(
            input_records,
            input_records,
            itertools.combinations(range(len(input_records)), 2)
        )
        self.assertEqual(len(results), 4)

        self.assertAlmostEqual(float(results[1][2]), answer, 4)

        # Test identical records for which one age exist and one doesn't
        input_records = [x for x in records if x['record_id'] in ['118210', '118211']]

        results = scorer.score_records(
            input_records,
            input_records,
            itertools.combinations(range(len(input_records)), 2)
        )
        self.assertEqual(len(results), 4)

        self.assertAlmostEqual(float(results[1][2]), 1.0, 1)
