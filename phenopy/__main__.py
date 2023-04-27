import fire
import itertools
import sys
from configparser import NoOptionError, NoSectionError
from phenopy.util import open_or_stdout
from phenopy.build_hpo import generate_annotated_hpo_network
from phenopy.config import config, logger
from phenopy.likelihood import predict_likelihood_moldx
from phenopy.score import Scorer
from phenopy.util import parse_input, half_product
from phenoseries.experiment import run_phenoseries_experiment


def score(input_file, output_file='-', records_file=None, annotations_file=None,
          custom_disease_file=None, ages_distribution_file=None, self=False,
          summarization_method='BMWA', scoring_method='HRSS', threads=1):
    """
    Scores similarity of provided HPO annotated entries (see format below) against a
    set of HPO annotated dataset. By default scoring happens against diseases
    annotated by the HPO group. See https://hpo.jax.org/app/download/annotation.

    Phenopy also supports scoring the product of provided entries (see "--product") or
    scoring against a custom records dataset (see "--records-file).

    :param input_file: File with HPO annotated entries, one per line (see format below).
    :param output_file: File path where to store the results. [default: - (stdout)]
    :param records_file: An entity-to-phenotype annotation file in the same format as
        "input_file". This file, if
     provided, is used to score entries in the "input_file" against entries here.
        [default: None]
    :param annotations_file: An entity-to-phenotype annotation file in the same format
        as "input_file". This file, if
     provided, is used to add information content to the network. [default: None]
    :param custom_disease_file: entity Annotation for ranking diseases/genes
    :param ages_distribution_file: Phenotypes age summary stats file containing
        phenotype HPO id, mean_age, and std. [default: None]
    :param self: Score entries in the "input_file" against itself.
    :param summarization_method: The method used to summarize the HRSS matrix.
        Supported Values are best match average
    (BMA), best match weighted average (BMWA), and maximum (maximum). [default: BMWA]
    :param scoring_method: Either HRSS or Resnik
    :param threads: Number of parallel processes to use. [default: 1]
    """

    try:
        obo_file = config.get('hpo', 'obo_file')
    except (NoSectionError, NoOptionError):
        logger.critical(
            'No HPO OBO file found in the configuration file. See "hpo:obo_file" '
            'parameter.')
        sys.exit(1)
    if custom_disease_file is None:
        try:
            disease_to_phenotype_file = config.get('hpo', 'disease_to_phenotype_file')
        except (NoSectionError, NoOptionError):
            logger.critical(
                'No HPO annotated dataset file found in the configuration file.'
                ' See "hpo:disease_to_phenotype_file" parameter.'
            )
            sys.exit(1)
    else:
        logger.info(f"using custom disease annotation file: {custom_disease_file}")
        disease_to_phenotype_file = custom_disease_file

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
        scorer = Scorer(hpo_network, summarization_method=summarization_method,
                        scoring_method=scoring_method)
    except ValueError as e:
        logger.critical(f'Failed to initialize scoring class: {e}')
        sys.exit(1)

    if self:
        score_records = input_records

        scoring_pairs = half_product(len(score_records), len(score_records))
    else:
        if records_file:
            score_records = parse_input(records_file, hpo_network, alt2prim)
        else:
            score_records = disease_records

        scoring_pairs = itertools.product(
            range(len(input_records)),
            range(len(score_records)),
        )

    results = scorer.score_records(input_records, score_records, scoring_pairs, threads)

    with open_or_stdout(output_file) as output_fh:
        output_fh.write('\t'.join(['#query', 'entity_id', 'score']))
        output_fh.write('\n')
        for result in results:
            output_fh.write('\t'.join(str(column) for column in result))
            output_fh.write('\n')


def validate_phenoseries(phenotypic_series_filepath, outdir=None, min_hpos=4,
                         min_entities=2, phenoseries_fraction=1.0,
                         scoring_method='HRSS', threads=1, omim_phenotypes_file="",
                         pairwise_mim_scores_file=""):
    """
    This runs the phenoseries experiment for a fraction of the OMIM phenoseries
    (PSid's). It Outputs a file with each row containing: PSid, MIMid, Python list
    of integers (ranks), and the length of the list.

    :param phenotypic_series_filepath: The phenotypicSeries.txt file from OMIM API.
        This is required to run validation.
    :param outdir: Directory where output files will be written.
    :param min_hpos: The minimum number of HPO ids annotated to a MIM id for the
        MIM id to be included in the experiment.
    :param min_entities: The minimum number of MIM ids for a phenoseries id to be
        included in the experiment.
    :param phenoseries_fraction: The fraction of total phenoseries to evaluate.
    :param scoring_method: Either HRSS, Resnik, Jaccard, or word2vec
    :param threads: Number of parallel processes to use. [default: 1]
    :param omim_phenotypes_file: <Optional> Path to the file containing OMIM id in
        the first column and a Python
     list of hpo ids in the second column.
    :param pairwise_mim_scores_file: <Optional> Path to the file containing similarity
        scores for each of the
    """
    run_phenoseries_experiment(
        outdir=outdir,
        phenotypic_series_filepath=phenotypic_series_filepath,
        min_hpos=min_hpos,
        min_entities=min_entities,
        phenoseries_fraction=phenoseries_fraction,
        scoring_method=scoring_method,
        threads=threads,
        omim_phenotypes_file=omim_phenotypes_file,
        pairwise_mim_scores_file=pairwise_mim_scores_file,
        )


def likelihood_moldx(input_file, output_file=None, k_phenotype_groups=1000):
    """
    :param input_file: The file path to a file containing three columns.
        [ID\tkey=value\thpodid,hpoid,hpoid]
    :param output_file: The file path to an output file containing the predicted
        probabilities
    :param k_phenotype_groups: The number of phenotype groups to use for encoding
        phenotypes. The CLI version of phenopy allows for one of [1000, 1500]
    """
    try:
        obo_file = config.get('hpo', 'obo_file')
    except (NoSectionError, NoOptionError):
        logger.critical(
            'No HPO OBO file found in the configuration file. See "hpo:obo_file" '
            'parameter.')
        sys.exit(1)
    try:
        disease_to_phenotype_file = config.get('hpo', 'disease_to_phenotype_file')
    except (NoSectionError, NoOptionError):
        logger.critical(
            'No HPO annotated dataset file found in the configuration file.'
            ' See "hpo:disease_to_phenotype_file" parameter.'
        )
        sys.exit(1)

    logger.info(f'Loading HPO OBO file: {obo_file}')
    hpo_network, alt2prim, _ = \
        generate_annotated_hpo_network(obo_file,
                                       disease_to_phenotype_file,
                                       )

    # parse input records
    input_records = parse_input(input_file, hpo_network, alt2prim)
    record_ids = [record["record_id"] for record in input_records]
    phenotypes = [record["terms"] for record in input_records]

    # predict likelihood of molecular diagnosis
    positive_probabilities = predict_likelihood_moldx(
        phenotypes,
        phenotype_groups=None, 
        hpo_network=hpo_network, 
        alt2prim=alt2prim,
        k_phenotype_groups=k_phenotype_groups,
        )

    if output_file is None:
        output_file = "phenopy.likelihood_moldx.txt"
    try:
        with open(output_file, "w") as f:
            for sample_id, probability in zip(record_ids, positive_probabilities):
                f.write(f"{sample_id}\t{probability}\n")
    except IOError:
        sys.exit("Something went wrong writing the probabilities to file")


def main():
    fire.Fire({
        'score': score,
        'likelihood': likelihood_moldx,
        'validate-phenoseries': validate_phenoseries,
    })


if __name__ == '__main__':
    main()
