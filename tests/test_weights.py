import pandas as pd

import os
import unittest

from phenopy import generate_annotated_hpo_network
from phenopy.config import logger
from phenopy.weights import get_truncated_normal, hpo_age_to_weight, make_age_distributions


class WeightsTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # parent dir
        cls.parent_dir = os.path.dirname(os.path.realpath(__file__))

        # load and process the network
        cls.obo_file = os.path.join(cls.parent_dir, 'data/hp.obo')
        cls.disease_to_phenotype_file = os.path.join(cls.parent_dir, 'data/phenotype.hpoa')
        cls.ages_distribution_file = os.path.join(cls.parent_dir, 'data/phenotype_age.tsv')
        cls.hpo_network, alt2prim, disease_records = \
            generate_annotated_hpo_network(cls.obo_file,
                                           cls.disease_to_phenotype_file,
                                           ages_distribution_file=cls.ages_distribution_file
                                           )

    def test_age_weights(self):
        self.assertEqual(hpo_age_to_weight(self.hpo_network, 'HP:0001251', 9.0), 1.0)
        self.assertAlmostEqual(hpo_age_to_weight(self.hpo_network, 'HP:0001251', 5.0), 0.726, 2)

    def test_make_age_distributions(self):
        with self.assertRaises(SystemExit) as se:
            make_age_distributions('notafilepath/notafile')

        syserr = se.exception
        self.assertEqual(syserr.code, 1)

        with self.assertRaises(SystemExit) as se:
            make_age_distributions('notafilepath/notafile', logger=logger)

        syserr = se.exception
        self.assertEqual(syserr.code, 1)

        ages_truth = pd.DataFrame([
            {'hpid': 'HP:0001251', 'age_dist': get_truncated_normal(6.0, 3.0, 0.0, 6.0)},
            {'hpid': 'HP:0001263', 'age_dist': get_truncated_normal(1.0, 1.0, 0.0, 1.0)},
            {'hpid': 'HP:0001290', 'age_dist': get_truncated_normal(1.0, 1.0, 0.0, 1.0)},
            {'hpid': 'HP:0004322', 'age_dist': get_truncated_normal(10.0, 3.0, 0.0, 10.0)},
            {'hpid': 'HP:0001249', 'age_dist': get_truncated_normal(6.0, 3.0, 0.0, 6.0)},
        ]).set_index('hpid')

        phenotype_ages_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/phenotype_age.tsv')
        df = make_age_distributions(phenotype_ages_file)
        self.assertEqual(set(ages_truth.index), set(df.index))

        for hpid in ages_truth.index:
            self.assertAlmostEqual(ages_truth.loc[hpid]['age_dist'].mean(), df.loc[hpid]['age_dist'].mean(), 1)

    def test_get_truncated_normal(self):
        self.assertAlmostEqual(get_truncated_normal(6.0, 1.0, 0.0, 6.0).mean(), 5.20, 2)
        self.assertAlmostEqual(get_truncated_normal(6.0, 1.0, 0.0, 6.0).cdf(3.0), 0.0026, 2)
        self.assertAlmostEqual(get_truncated_normal(6.0, 1.0, 0.0, 6.0).cdf(12.0), 1.0, 2)






