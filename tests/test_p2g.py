import unittest

from phenosim.p2g import load


class P2GTestCase(unittest.TestCase):
    def test_load(self):
        with self.assertRaises(SystemExit) as se:
            load('notafilepath/notafile')

        syserr = se.exception
        self.assertEqual(syserr.code, 1)
