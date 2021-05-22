
import os
import unittest
import pandas as pd
from phenopy.build_hpo import generate_annotated_hpo_network
from phenopy.cluster import Cluster


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
        data = pd.read_table(cls.input_file, names=['record_id', 'info', 'hpo_terms'])
        data['hpo_terms'] = data['hpo_terms'].str.split("|")
        cls.cluster = Cluster(data, scoring_method='Jaccard', kfile=cls.kfile, )

    def test_process_kfile(self):

        self.cluster.process_kfile(k=1000)
        self.assertEqual(self.cluster.n_features, 4)
        self.assertEqual(self.cluster.feature_to_hps, {0: ['HP:0003011'], 1: ['HP:0001263'], 2: ['HP:0001290'], 3: ['HP:0001251']})
        self.assertEqual(self.cluster.hp_to_feature, {'HP:0001251': 3, 'HP:0001263': 1, 'HP:0001290': 2, 'HP:0003011': 0})

    def test_prep_cluster_data(self):
        self.cluster.process_kfile(k=1000)
        self.cluster.prep_cluster_data()
        self.assertEqual(self.cluster.data.shape[0], 398)

    def test_prep_feature_array(self):
        self.cluster.process_kfile(k=1000)
        self.cluster.prep_cluster_data()
        self.cluster.prep_feature_array()
        self.assertEqual(self.cluster.feature_array.shape[0], 398)

    def test_apply_umap(self):
        self.cluster.process_kfile(k=1000)
        self.cluster.prep_cluster_data()
        self.cluster.prep_feature_array()
        self.cluster.umap()
        self.assertEqual(self.cluster.data[['umap1','umap2']].values.shape[0], 398)

    def test_dbscan(self):
        self.cluster.process_kfile(k=1000)
        self.cluster.prep_cluster_data()
        self.cluster.prep_feature_array()
        self.cluster.umap()
        self.cluster.dbscan()

        self.assertEqual(self.cluster.data['cluster_id'].shape[0], 398)
        self.assertEqual(self.cluster.dbscan_stats['n_clusters'], 2)
        self.assertLess(self.cluster.dbscan_stats['n_noise'], 10)
        self.assertAlmostEqual(self.cluster.dbscan_stats['silhouette_score'], 0.615, 1)

    def test_hdbscan(self):

        self.cluster.umap(metric='precomputed')
        self.cluster.hdbscan()
        self.assertEqual(self.cluster.data['cluster_id'].shape[0], 398)








