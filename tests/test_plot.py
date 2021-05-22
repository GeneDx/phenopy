
import os
import unittest
import pandas as pd

from phenopy.build_hpo import generate_annotated_hpo_network
from phenopy.cluster import Cluster


class PlotTestCase(unittest.TestCase):
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
        cls.cluster = Cluster(data, scoring_method='Jaccard', kfile=cls.kfile, k=1000)
        cls.save_basic_dbscan_plot = os.path.join(cls.parent_dir, 'data/basic_dbscan_plot.png')

    def test_plot_basic_dbscan(self):
        self.cluster.process_kfile()
        self.cluster.prep_cluster_data()
        self.cluster.prep_feature_array()
        self.cluster.umap()
        self.cluster.dbscan()
        plot_ = self.cluster.visualize(color_by="cluster_id")
        plot_.savefig(self.save_basic_dbscan_plot)
        # TODO assert file exist