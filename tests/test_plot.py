
import os
import unittest

from phenopy import generate_annotated_hpo_network
from phenopy.cluster import prep_cluster_data, process_kfile, prep_feature_array, apply_umap, dbscan
from phenopy.plot import plot_basic_dbscan


class PlotTestCase(unittest.TestCase):
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
        cls.input_file = os.path.join(cls.parent_dir, 'data/cluster_two_groups.txt')
        cls.kfile = os.path.join(cls.parent_dir, 'data/phenotype_groups_4.txt')
        cls.save_basic_dbscan_plot = os.path.join(cls.parent_dir, 'data/basic_dbscan_plot.png')

    def test_plot_basic_dbscan(self):
        feature_to_hps, hp_to_feature, n_features = process_kfile(self.kfile, k=1000)
        test_df = prep_cluster_data(self.input_file, self.kfile)
        X_vect = prep_feature_array(test_df, n_features)
        umap_result = apply_umap(X_vect)
        labels, core_samples_mask, stats = dbscan(umap_result)
        result = plot_basic_dbscan(umap_result, core_samples_mask, labels)
        result.savefig(self.save_basic_dbscan_plot)
        # TODO assert file exist