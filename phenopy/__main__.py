import fire
import itertools

from configparser import NoOptionError, NoSectionError
from multiprocessing import Pool

from phenopy import open_or_stdout, parse_input
from phenopy.config import config, logger
from phenopy.network import load as load_network
from phenopy.d2p import load as load_d2p
from phenopy.score import Scorer
from phenopy.util import remove_parents


def score(input_file, output_file='-', records_file=None, annotations_file=None, ages_distribution_file=None,
          self=False, threads=1):
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

    # parse input records
    input_records = parse_input(input_file)

    # load phenotypes to diseases associations
    (
        score_records,
        phenotype_to_diseases,
        phenotype_disease_frequencies,
    ) = load_d2p(disease_to_phenotype_file)

    # load hpo network
    hpo_network = load_network(
        obo_file,
        phenotype_to_diseases,
        len(score_records),
        annotations_file=annotations_file,
        ages_distribution_file=ages_distribution_file,
        phenotype_disease_frequencies=phenotype_disease_frequencies,
    )

    # create instance the scorer class
    scorer = Scorer(hpo_network)

    # TODO: refactor cleanup into a generic function so we can called for all records dataset
    for record in input_records:
        # TODO: can we combine these cleanup operations so we don't have to iterate over all terms multiple times
        # prune parent terms
        record['terms'] = remove_parents(record['terms'], hpo_network)

        # convert and filter the query hpo ids
        record['terms'] = scorer.convert_alternate_ids(record['terms'])
        record['terms'] = scorer.filter_and_sort_hpo_ids(record['terms'])

    if self:
        score_records = input_records

        scoring_pairs = itertools.combinations(
            range(len(input_records)),
            2,
        )
    else:
        if records_file:
            score_records = parse_input(records_file)

            for record in score_records:
                # prune parent terms if requested
                record['terms'] = remove_parents(record['terms'], hpo_network)

                # convert and filter the query hpo ids
                record['terms'] = scorer.convert_alternate_ids(record['terms'])
                record['terms'] = scorer.filter_and_sort_hpo_ids(record['terms'])
        else:
            # convert and filter disease to phenotypes records terms
            # TODO: we should probably do this during d2p
            for record in score_records:
                record['terms'] = scorer.convert_alternate_ids(record['terms'])
                record['terms'] = scorer.filter_and_sort_hpo_ids(record['terms'])

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
                    False,  # use weights
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

    return results


def main():
    fire.Fire({
        'score': score,
    })


if __name__ == '__main__':
    main()
