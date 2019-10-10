import numpy as np
import os
import unittest

from phenopy.config import config
from phenopy.ic import calculate_information_content
from phenopy.obo import process
from phenopy.obo import load as load_obo
from phenopy.d2p import load as load_d2p
from phenopy.util import export_phenotype_hpoa_with_no_parents


class ICTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # parent dir
        cls.parent_dir = os.path.dirname(os.path.realpath(__file__))
        cls.hpo_network_file = os.path.join(cls.parent_dir, 'data/hpo_network.pickle')

        # load phenotypes to genes associations
        cls.disease_to_phenotype_file = os.path.join(cls.parent_dir, 'data/phenotype.hpoa')
        cls.disease_to_phenotypes, cls.phenotype_to_diseases = load_d2p(cls.disease_to_phenotype_file)

        # load and process the network
        config.set('hpo', 'data_directory', os.path.join(cls.parent_dir, 'data'))
        cls.obo_file = os.path.join(cls.parent_dir, 'data/hp.obo')
        cls.hpo_network = load_obo(cls.obo_file)
        cls.num_diseases_annotated = len(cls.disease_to_phenotypes)
        cls.hpo_network = process(cls.hpo_network, cls.phenotype_to_diseases, cls.num_diseases_annotated,
                                  custom_annotations_file=None)
        cls.hpo_id = 'HP:0010863'
        cls.disease_to_phenotype_output_file = os.path.join(cls.parent_dir, 'data/phenotype.noparents.hpoa')

    def tearDown(cls):
        if os.path.exists(cls.disease_to_phenotype_output_file):
            os.remove(cls.disease_to_phenotype_output_file)
        if os.path.exists(cls.hpo_network_file):
            os.remove(cls.hpo_network_file)

    def test_ic_d2p(self):
        """Calculate the information content of a phenotype"""
        self.assertAlmostEqual(self.hpo_network.node[self.hpo_id]['ic'], 5.69, 2)

    def test_ic_custom(self):
        """Calculate the information content of a phenotype when multiple annotations are present"""
        custom_annotation_file = os.path.join(self.parent_dir, 'data/test.score-product.txt')
        hpo_network = load_obo(self.obo_file)
        hpo_network = process(hpo_network, self.phenotype_to_diseases, self.num_diseases_annotated,
                              custom_annotations_file=custom_annotation_file)

        self.assertAlmostEqual(hpo_network.node[self.hpo_id]['ic'], 6.38, 1)

    def test_inf_ic(self):
        inf_ic = calculate_information_content(
            self.hpo_id,
            self.hpo_network,
            self.phenotype_to_diseases,
            1e310,
            None,
        )
        self.assertAlmostEqual(inf_ic, -np.log(np.nextafter(0, 1)))

    def test_ic_d2p_no_parents(self):
        export_phenotype_hpoa_with_no_parents(self.disease_to_phenotype_file, self.disease_to_phenotype_output_file, self.hpo_network)
        self.assertTrue(os.path.exists(self.disease_to_phenotype_output_file))
