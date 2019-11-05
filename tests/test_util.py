import os
import unittest

from phenopy.util import parse, read_records_file


class UtilTestCase(unittest.TestCase):
    def test_read_records_file(self):
        with self.assertRaises(SystemExit) as se:
            read_records_file('notafilepath/notafile')

        syserr = se.exception
        self.assertEqual(syserr.code, 1)
        records_truth = [
            {
                'sample': '118200',
                'age': 9.0,
                'gender': 'Female',
                'terms': 'HP:0001263|HP:0001251|HP:0001290|HP:0004322'.split('|')
             },
            {
                'sample': '118210',
                'age': 4.0,
                'gender': None,
                'terms': 'HP:0001249|HP:0001263|HP:0001290'.split('|')
            },
            {
                'sample': '118211',
                'age': None,
                'gender': None,
                'terms': 'HP:0001249|HP:0001263|HP:0001290'.split('|')
            }
        ]
        records_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/test.score-short.txt')
        records = read_records_file(records_path, no_parents=False)
        self.assertEqual(records,records_truth)


    def test_parse(self):
        string = 'age=13;sex=Male'
        self.assertEqual(parse(string, what='sex'), 'Male')
        self.assertEqual(parse(string, what='age'), 13.0)

        string = 'age=13.64;sex=male'
        self.assertEqual(parse(string, what='sex'), 'Male')
        self.assertEqual(parse(string, what='age'), 13.6)

        string = 'age=12.9;sex=female'
        self.assertEqual(parse(string, what='sex'), 'Female')
        self.assertEqual(parse(string, what='age'), 12.9)

        string = 'sex=Female'
        self.assertEqual(parse(string, what='sex'), 'Female')

        string = 'sex=FEMALE'
        self.assertEqual(parse(string, what='sex'), 'Female')

        string = 'sex=F'
        self.assertEqual(parse(string, what='sex'), 'Female')

        string = 'age=1'
        self.assertEqual(parse(string, what='age'), 1.0)

        string = '.'
        self.assertEqual(parse(string, what='age'), None)

        string = '. '
        self.assertEqual(parse(string, what='age'), None)

        string = ' . '
        self.assertEqual(parse(string, what='age'), None)

        string = '13?'
        self.assertEqual(parse(string, what='age'), None)

        string = 'sex=NA'
        self.assertEqual(parse(string, what='sex'), None)

        string = 'sex=Unknown'
        self.assertEqual(parse(string, what='sex'), None)





