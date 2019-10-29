__project__ = 'phenopy'
__version__ = '0.2.1'

import csv
import sys
from contextlib import contextmanager

from phenopy.config import logger


@contextmanager
def open_or_stdout(filename):
    if filename != '-':
        with open(filename, 'w') as f:
            yield f
    else:
        yield sys.stdout


def parse_input(input_file):
    """
    Parse input file.
    """
    try:
        with open(input_file, 'r') as input_fh:
            reader = csv.reader(filter(lambda l: not l.startswith('#'), input_fh), delimiter='\t')
            records = []
            for line in reader:
                record = {
                    'record_id': line[0],
                    'terms': line[2].split('|'),
                    'weights': {},
                    **dict(item.split('=') for item in line[1].split(';') if line[1] != '.')
                }


    except (FileNotFoundError, PermissionError) as e:
        logger.critical(f'Provided input file could not be loaded or does not exist: {e}')
        exit(1)
    except ValueError:
        logger.critical(f'Unable to parse input file, invalid line number: {reader.line_num}:{input_file}')
        exit(1)

    return records
