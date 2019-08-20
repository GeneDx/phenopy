import unittest

from phenosim.obo import load


class OboTestCase(unittest.TestCase):
    def test_load(self):
        with self.assertRaises(SystemExit) as se:
            load('notafilepath/notafile')

        syserr = se.exception
        self.assertEqual(syserr.code, 1)

        hpo_network = load('data/hp.obo')
        self.assertEqual(len(hpo_network), 22)
        self.assertEqual(hpo_network.node['HP:0000009']['name'], 'Functional abnormality of the bladder')
