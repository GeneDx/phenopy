import pandas as pd

import os
import unittest

from phenopy.weights import get_truncated_normal, age_to_weights, make_age_distributions


class UtilTestCase(unittest.TestCase):

    def test_age_weights(self):
        self.assertEqual(age_to_weights(get_truncated_normal(6.0, 1.0, 0.0, 6.0), 9.0), 1.0)
        self.assertEqual(age_to_weights(get_truncated_normal(9.0, 1.0, 0.0, 9.0), 9.0), 1.0)
        self.assertEqual(age_to_weights(get_truncated_normal(9.0, 1.0, 0.0, 9.0), 9.0), 1.0)
        self.assertAlmostEqual(age_to_weights(get_truncated_normal(9.0, 1.0, 0.0, 9.0), 8.0), 0.317, 2)

    def test_make_age_distributions(self):
        with self.assertRaises(SystemExit) as se:
            make_age_distributions('notafilepath/notafile')

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






