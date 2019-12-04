import fire
import itertools
import os
import pickle
import sys

import numpy as np

from configparser import NoOptionError, NoSectionError
from multiprocessing import Pool

from phenopy import open_or_stdout, generate_annotated_hpo_network
from phenopy.config import config, data_directory, logger
from phenopy.hpo_array import generate_hpo_array

from phenopy.score import Scorer
from phenopy.util import parse_input, half_product


def embed(input_file, output_file='-', records_file=None, annotations_file=None, ages_distribution_file=None,
          threads=1):
    """
    Given a list of samples in the input file, embed the samples as phenotype similarity vectors. The output is a
    binary numpy array file. The array has as many columns as there are nodes in the HPO and the number of rows is equal
    to the number of samples in the input file.
    :param input_file:
    :param output_file:
    :param records_file:
    :param annotations_file:
    :param ages_distribution_file:
    :param threads:
    :return: A numpy array if output_
    """
    pass


def score_all_phenotypes(outdir=None, algorithm='wup', threads=1):
    """
    Calculates the HRSS score for every term in the HPO and stores the scores as a pairwise numpy array. Also
    outputs a dictionary mapping the integer index of the HPO term to the HPO name. This runs in ~1 hour on 48 cpus.
    :param outdir: Path to output the pairwise scores array.
    :param algorithm: The semantic similarity algorithm used to score two HPO terms. {'wup', 'hrss'}
    :param threads: Number of threads to calculate the pairwise HPO HRSS array.
    :return: A numpy array of pairwise HRSS scores between every term in HPO. The integer index of the array can be
    mapped back to HPO terms with the int2hpo.pkl file.
    """
    if outdir is None:
        outdir = data_directory

    scores_arr, int2hpo = generate_hpo_array(threads, algorithm=algorithm)

    try:
        # store the mapping so we can recover the indices
        with open(os.path.join(outdir, 'int2hpo.pkl'), 'wb') as f:
            pickle.dump(int2hpo, f)

        # store the array
        np.save(os.path.join(outdir, f'pairwise.{algorithm}.npy'), scores_arr)

    except NotADirectoryError:
        logger.critical('Please choose a valid directory to store results in or set outdir=None.')


def score(input_file, output_file='-', records_file=None, annotations_file=None, ages_distribution_file=None,
          self=False, summarization_method='BMWA', threads=1):
    """
    Scores similarity of provided HPO annotated entries (see format below) against a set of HPO annotated dataset. By
    default scoring happens against diseases annotated by the HPO group. See https://hpo.jax.org/app/download/annotation.

    Phenopy also supports scoring the product of provided entries (see "--product") or scoring against a custom records
    dataset (see "--records-file).

    :param input_file: File with HPO annotated entries, one per line (see format below).
    :param output_file: File path where to store the results. [default: - (stdout)]
    :param records_file: An entity-to-phenotype annotation file in the same format as "input_file". This file, if
     provided, is used to score entries in the "input_file" against entries here. [default: None]
    :param annotations_file: An entity-to-phenotype annotation file in the same format as "input_file". This file, if
     provided, is used to add information content to the network. [default: None]
    :param ages_distribution_file: Phenotypes age summary stats file containing phenotype HPO id, mean_age, and std.
     [default: None]
    :param self: Score entries in the "input_file" against itself.
    :param summarization_method: The method used to summarize the HRSS matrix. Supported Values are best match average
    (BMA), best match weighted average (BMWA), maximum (maximum), and vector (vector). [default: BMWA]
    :param threads: Number of parallel processes to use. [default: 1]
    """

    try:
        obo_file = config.get('hpo', 'obo_file')
    except (NoSectionError, NoOptionError):
        logger.critical(
            'No HPO OBO file found in the configuration file. See "hpo:obo_file" parameter.')
        exit(1)

    try:
        disease_to_phenotype_file = config.get('hpo', 'disease_to_phenotype_file')
    except (NoSectionError, NoOptionError):
        logger.critical(
            'No HPO annotated dataset file found in the configuration file.'
            ' See "hpo:disease_to_phenotype_file" parameter.'
        )
        exit(1)

    logger.info(f'Loading HPO OBO file: {obo_file}')
    hpo_network, alt2prim, disease_records = \
        generate_annotated_hpo_network(obo_file,
                                       disease_to_phenotype_file,
                                       annotations_file=annotations_file,
                                       ages_distribution_file=ages_distribution_file
                                       )

    # parse input records
    input_records = parse_input(input_file, hpo_network, alt2prim)

    # create instance the scorer class
    try:
        scorer = Scorer(hpo_network, summarization_method=summarization_method)
    except ValueError as e:
        logger.critical(f'Failed to initialize scoring class: {e}')
        sys.exit(1)

    if self:
        score_records = input_records

        scoring_pairs = list(half_product(len(score_records), len(score_records)))
    else:
        if records_file:
            score_records = parse_input(records_file, hpo_network, alt2prim)
        else:
            score_records = disease_records

        scoring_pairs = itertools.product(
            range(len(input_records)),
            range(len(score_records)),
        )

    # launch as many scoring process as requested
    with Pool(threads) as p:
        results = p.starmap(
            scorer.score_records,
            [
                (
                    input_records,  # a records
                    score_records,  # b records
                    scoring_pairs,  # pairs
                    i,  # thread_index
                    threads,  # threads
                ) for i in range(threads)
            ]
        )

    with open_or_stdout(output_file) as output_fh:
        output_fh.write('\t'.join(['#query', 'entity_id', 'score']))
        output_fh.write('\n')
        for r in results:
            for s in r:
                output_fh.write('\t'.join(s))
                output_fh.write('\n')


def main():
    fire.Fire({
        'embed': embed,
        'score': score,
        'score-all-phenotypes': score_all_phenotypes,
    })


if __name__ == '__main__':
    main()
