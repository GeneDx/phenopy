import os
import unittest

from phenopy.config import config
from phenopy.d2p import load as load_d2p
from phenopy.network import load as load_network
from phenopy.network import annotate
from phenopy.util import generate_alternate_ids


class NetworkTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parent_dir = os.path.dirname(os.path.realpath(__file__))
        config.set('hpo', 'data_directory', os.path.join(cls.parent_dir, 'data'))
        cls.obo_file = os.path.join(cls.parent_dir, 'data/hp.obo')

    def test_load_network(self):
        hpo_network = load_network(self.obo_file)
        self.assertEqual(len(hpo_network), 28)

    def test_annotate_network(self):
        hpo_network = load_network(self.obo_file)
        alt2prim = generate_alternate_ids(hpo_network)

        # load phenotypes to diseases associations
        disease_to_phenotype_file = os.path.join(self.parent_dir, 'data/phenotype.hpoa')
        disease_records, phenotype_to_diseases = load_d2p(disease_to_phenotype_file, hpo_network, alt2prim)

        num_diseases_annotated = len(disease_records)
        hpo_network = annotate(hpo_network, phenotype_to_diseases, num_diseases_annotated, alt2prim)

        self.assertAlmostEqual(hpo_network.nodes['HP:0010863']['ic'], 5.69, 2)
