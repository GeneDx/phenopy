import unittest

from phenopy.d2p import load


class DPGTestCase(unittest.TestCase):
    def test_load(self):
        with self.assertRaises(SystemExit) as se:
            load('notafilepath/notafile')

        syserr = se.exception
        self.assertEqual(syserr.code, 1)
