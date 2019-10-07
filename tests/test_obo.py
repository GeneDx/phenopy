import os
import unittest

from phenopy.obo import load


class OboTestCase(unittest.TestCase):
    def test_load(self):
        with self.assertRaises(SystemExit) as se:
            load('notafilepath/notafile')

        syserr = se.exception
        self.assertEqual(syserr.code, 1)

        hpo_network = load(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/hp.obo'))
        self.assertEqual(len(hpo_network), 28)
        self.assertEqual(hpo_network.node['HP:0012434']['name'], 'Delayed social development')
