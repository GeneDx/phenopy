
import os
import unittest
import pandas as pd
from phenopy.util import parse_input
from phenopy import generate_annotated_hpo_network
from phenopy.cluster import prep_cluster_data, process_kfile, prep_feature_array, apply_umap, dbscan, compute_dx_yield


class ClusterTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # parent dir
        cls.parent_dir = os.path.dirname(os.path.realpath(__file__))

        # load and process the network
        cls.obo_file = os.path.join(cls.parent_dir, 'data/hp.obo')
        cls.disease_to_phenotype_file = os.path.join(cls.parent_dir, 'data/phenotype.hpoa')
        cls.ages_distribution_file = os.path.join(cls.parent_dir, 'data/phenotype_age.tsv')
        cls.hpo_network, cls.alt2prim, cls.disease_records = \
            generate_annotated_hpo_network(cls.obo_file,
                                           cls.disease_to_phenotype_file,
                                           ages_distribution_file=cls.ages_distribution_file
                                           )
        cls.input_file = os.path.join(cls.parent_dir, 'data/cluster_two_groups.txt')
        cls.kfile = os.path.join(cls.parent_dir, 'data/phenotype_groups_4.txt')

    def test_process_kfile(self):
        feature_to_hps, hp_to_feature, n_features = process_kfile(self.kfile, k=1000)
        self.assertEqual(n_features, 4)
        self.assertEqual(feature_to_hps, {0: 'HP:0003011', 1: 'HP:0001263', 2: 'HP:0001290', 3: 'HP:0001251'})
        self.assertEqual(hp_to_feature, {'HP:0001251': 3, 'HP:0001263': 1, 'HP:0001290': 2, 'HP:0003011': 0})

    def test_prep_cluster_data(self):
        feature_to_hps, hp_to_feature, n_features = process_kfile(self.kfile, k=1000)
        records = parse_input(self.input_file, self.hpo_network, self.alt2prim)
        df_input = pd.DataFrame.from_dict(records)
        result = prep_cluster_data(df_input, hp_to_feature)
        self.assertEqual(result.shape[0], 398)

    def test_prep_feature_array(self):
        feature_to_hps, hp_to_feature, n_features = process_kfile(self.kfile, k=1000)
        records = parse_input(self.input_file, self.hpo_network, self.alt2prim)
        df_input = pd.DataFrame.from_dict(records)
        result = prep_cluster_data(df_input, hp_to_feature)
        X_vect = prep_feature_array(result, n_features)
        self.assertEqual(X_vect.shape[0], 398)

    def test_apply_umap(self):
        feature_to_hps, hp_to_feature, n_features = process_kfile(self.kfile, k=1000)
        records = parse_input(self.input_file, self.hpo_network, self.alt2prim)
        df_input = pd.DataFrame.from_dict(records)
        result = prep_cluster_data(df_input, hp_to_feature)
        X_vect = prep_feature_array(result, n_features)
        umap_result = apply_umap(X_vect)
        self.assertEqual(umap_result.shape[0], 398)

    def test_dbscan(self):
        feature_to_hps, hp_to_feature, n_features = process_kfile(self.kfile, k=1000)
        records = parse_input(self.input_file, self.hpo_network, self.alt2prim)
        df_input = pd.DataFrame.from_dict(records)
        result = prep_cluster_data(df_input, hp_to_feature)
        X_vect = prep_feature_array(result, n_features)
        umap_result = apply_umap(X_vect)
        labels, core_samples_mask, stats = dbscan(umap_result)
        self.assertEqual(len(labels), 398)
        self.assertEqual(stats['n_clusters'], 2)
        self.assertLess(stats['n_noise'], 10)
        self.assertAlmostEqual(stats['silhouette_score'], 0.619, 1)

    def test_compute_dx_yield(self):
        records = parse_input(self.input_file, self.hpo_network, self.alt2prim)
        df_input = pd.DataFrame.from_dict(records)
        df_input['cluster_id'] = df_input['is_diagnosed']
        dx_yield = compute_dx_yield(df_input)
        self.assertEqual(dx_yield.loc[0, 'dx_yield'], 100.0)
        self.assertEqual(dx_yield.loc[1, 'dx_yield'],   0.0)





