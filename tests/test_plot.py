
import os
import unittest
import pandas as pd

from phenopy import generate_annotated_hpo_network
from phenopy.cluster import prep_cluster_data, process_kfile, prep_feature_array, apply_umap, dbscan
from phenopy.plot import plot_basic_dbscan
from phenopy.util import parse_input


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
        cls.save_basic_dbscan_plot = os.path.join(cls.parent_dir, 'data/basic_dbscan_plot.png')

    def test_plot_basic_dbscan(self):
        feature_to_hps, hp_to_feature, n_features = process_kfile(self.kfile, k=1000)
        records = parse_input(self.input_file, self.hpo_network, self.alt2prim)
        df_input = pd.DataFrame.from_dict(records)
        result = prep_cluster_data(df_input, hp_to_feature)
        X_vect = prep_feature_array(result, n_features)
        umap_result = apply_umap(X_vect)
        labels, core_samples_mask, stats = dbscan(umap_result)
        plot_ = plot_basic_dbscan(umap_result, core_samples_mask, labels)
        plot_.savefig(self.save_basic_dbscan_plot)
        # TODO assert file exist