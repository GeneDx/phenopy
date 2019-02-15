import unittest

from phenosim.obo import load


class OboTestCase(unittest.TestCase):
    def test_load(self):
        with self.assertRaises(SystemExit) as se:
            load('notafilepath/notafile')

        syserr = se.exception
        self.assertEqual(syserr.code, 1)
