import numpy as np
import os
import unittest

from phenopy.config import config
from phenopy.ic import calculate_information_content
from phenopy.obo import process
from phenopy.obo import load as load_obo
from phenopy.p2g import load as load_p2g
from phenopy.util import export_pheno2genes_with_no_parents


class ScorerTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # parent dir
        cls.parent_dir = os.path.dirname(os.path.realpath(__file__))

        # load phenotypes to genes associations
        cls.pheno2genes_file = os.path.join(cls.parent_dir, 'data/phenotype_to_genes.txt')
        cls.terms_to_genes, cls.genes_to_terms, cls.num_genes_annotated = load_p2g(cls.pheno2genes_file)

        # load and process the network
        config.set('hpo', 'data_directory', os.path.join(cls.parent_dir, 'data'))
        cls.obo_file = os.path.join(cls.parent_dir, 'data/hp.obo')
        cls.hpo_network = load_obo(cls.obo_file)
        cls.hpo_network = process(cls.hpo_network, cls.terms_to_genes, cls.num_genes_annotated,
                                  custom_annotations_file=None)
        cls.hpo_id = 'HP:0010863'
        cls.pheno2genes_output_file = os.path.join(cls.parent_dir, 'data/phenotype_to_genes.noparents.txt')

    def tearDown(cls):
        if os.path.exists(cls.pheno2genes_output_file):
            os.remove(cls.pheno2genes_output_file)

    def test_ic_p2g(self):
        """Calculate the information content of a phenotype"""
        self.assertAlmostEqual(self.hpo_network.node[self.hpo_id]['ic'], 6.20, 1)

    def test_ic_custom(self):
        """Calculate the information content of a phenotype when multiple annotations are present"""
        custom_annotation_file = os.path.join(self.parent_dir, 'data/test.score-product.txt')
        hpo_network = load_obo(self.obo_file)
        hpo_network = process(hpo_network, self.terms_to_genes, self.num_genes_annotated,
                              custom_annotations_file=custom_annotation_file)

        self.assertAlmostEqual(hpo_network.node[self.hpo_id]['ic'], 7.17, 1)

    def test_inf_ic(self):
        inf_ic = calculate_information_content(
            self.hpo_id,
            self.hpo_network,
            self.terms_to_genes,
            1e310,
            None,
        )
        self.assertAlmostEqual(inf_ic, -np.log(np.nextafter(0, 1)))

    def test_ic_p2g_no_parents(self):
        export_pheno2genes_with_no_parents(self.pheno2genes_file, self.pheno2genes_output_file, self.hpo_network)
        self.assertTrue(os.path.exists(self.pheno2genes_output_file))

