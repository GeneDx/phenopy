import os
import unittest

from phenosim.obo import process
from phenosim.obo import load as load_obo
from phenosim.p2g import load as load_p2g

class ScorerTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # parent dir
        cls.parent_dir = os.path.dirname(os.path.realpath(__file__))

        # load phenotypes to genes associations
        pheno2genes_file = os.path.join(
            cls.parent_dir, 'data/phenotypes_to_genes.txt')
        cls.terms_to_genes, cls.genes_to_terms, cls.annotations_count = load_p2g(
            pheno2genes_file)

        # load and process the network
        cls.obo_file = os.path.join(cls.parent_dir, 'data/hp.obo')
        cls.hpo_network = load_obo(cls.obo_file)
        cls.hpo_network = process(cls.hpo_network, cls.terms_to_genes, cls.annotations_count,
                                  custom_annotation_files=None)

    def test_ic_p2g(self):
        """Calculate the information content of a phenotype"""
        self.assertEqual(self.hpo_network.node['HP:0012372']['ic'], 1.38, 2)


    def test_ic_custom(self):
        """Calculate the information content of a phenotype when multiple annotations are present"""
        custom_annotation_file = os.path.join(self.parent_dir, 'data/test.score-product.txt')
        hpo_network = load_obo(self.obo_file)
        hpo_network = process(hpo_network, self.terms_to_genes, self.annotations_count,
                              custom_annotation_files=[custom_annotation_file])

        self.assertAlmostEqual(hpo_network.node['HP:0012372']['ic'], 1.93, 2)

