import itertools
import numpy as np
import os
import pandas as pd
import pytest

from phenopy.network import annotate
from phenopy.score import Scorer
from phenopy.util import (
    remove_parents,
    parse_input,
    half_product,
)
from phenopy.weights import calculate_age_weights


def test_find_lca(test_data):
    lca = test_data["scorer"].find_lca("HP:0001249", "HP:0012434")
    assert lca == "HP:0012759"

    root_lca = test_data["scorer"].find_lca("HP:0012759", "HP:0000001")
    assert root_lca == "HP:0000001"

    parent_lca = test_data["scorer"].find_lca("HP:0012759", "HP:0012758")
    assert parent_lca == "HP:0012759"

    parent_lca = test_data["scorer"].find_lca("HP:0012758", "HP:0012759")
    assert parent_lca == "HP:0012759"

    parent_lca = test_data["scorer"].find_lca("HP:0012759", "HP:0012759")
    assert parent_lca == "HP:0012759"

    parent_lca = test_data["scorer"].find_lca("HP:0012759", "HP:0000750")
    assert parent_lca == "HP:0012759"


def test_calculate_gamma(test_data):
    t1 = "HP:0012758"
    t2 = "HP:0012759"

    gamma0 = test_data["scorer"].calculate_gamma(t1, t1, t2)
    assert gamma0 == 0

    gamma1a = test_data["scorer"].calculate_gamma(t1, t2, t2)
    assert gamma1a == 1

    gamma1b = test_data["scorer"].calculate_gamma(t2, t1, t2)
    assert gamma1b == 1

    gamma2 = test_data["scorer"].calculate_gamma("HP:0000750", "HP:0012434", t1)
    assert gamma2 == 2


def test_calculate_beta(test_data):
    t1 = "HP:0001344"
    t2 = "HP:0012759"
    beta = test_data["scorer"].calculate_beta(t1, t2)
    assert round(beta, 2) == 3.99


def test_score_hpo_pair_hrss(test_data):
    t1 = "HP:0011351"
    t2 = "HP:0012434"

    score = test_data["scorer"].score_hpo_pair_hrss(t1, t2)
    assert round(score, 2) == 0.14

    score = test_data["scorer"].score_hpo_pair_hrss(t1, t2)
    assert round(score, 2) == 0.14

    score = test_data["scorer"].score_hpo_pair_hrss(t2, t1)
    assert round(score, 2) == 0.14


def test_score(test_data):
    record_a = {
        "record_id": "sample_1",
        "terms": ["HP:0012433", "HP:0012434"],
        "weights": {},
    }
    record_b = {"record_id": "sample_2", "terms": [], "weights": {}}

    score0 = test_data["scorer"].score(record_a, record_b)
    assert score0[2] == 0.0
    record_b["terms"] = ["HP:0001249", "HP:0012758"]

    score_bma = test_data["scorer"].score(record_a, record_b)
    assert round(score_bma[2], 2) == 0.09
    test_data["scorer"].summarization_method = "maximum"
    score_max = test_data["scorer"].score(record_a, record_b)
    assert round(score_max[2], 4) == 0.1251

    test_data["scorer"].summarization_method = "not_a_method"
    with pytest.raises(ValueError):
        test_data["scorer"].score(record_a, record_b)

    record_a.update(
        {
            "terms": [
                "HP:0001251",
                "HP:0001263",
                "HP:0001290",
                "HP:0004322",
                "HP:0012433",
            ],
            "weights": {"age": [0.67, 1.0, 1.0, 0.4, 0.4]},
        }
    )
    record_b.update(
        {
            "terms": ["HP:0001249", "HP:0001263", "HP:0001290"],
            "weights": {"age": [1.0, 1.0, 1.0]},
        }
    )

    test_data["scorer"].summarization_method = "BMWA"
    test_data["scorer"].min_score_mask = 0.05
    score_bmwa = test_data["scorer"].score(record_a, record_b)
    assert round(score_bmwa[2], 4) == 0.1822

    record_a.update(
        {
            "terms": ["HP:0001251", "HP:0001263", "HP:0001290", "HP:0004322"],
            "weights": {"age": [0.67, 1.0, 1.0, 0.4]},
        }
    )
    record_b.update(
        {
            "terms": ["HP:0001263", "HP:0001249", "HP:0001290"],
            "weights": {"age": [1.0, 1.0, 0.5]},
        }
    )

    scorer = test_data["scorer"]
    scorer.summarization_method = "BMWA"

    score_bwma_both_weights = scorer.score(record_a, record_b)
    assert round(score_bwma_both_weights[2], 4) == 0.1918

    scorer.min_score_mask = None
    record_a["weights"].pop("age", None)
    score_bwma_one_weights = scorer.score(record_a, record_b)
    assert round(score_bwma_one_weights[2], 4) == 0.155


def test_score_records(test_data):
    query_name = "SAMPLE"
    query_terms = [
        "HP:0000750",
        "HP:0010863",
    ]
    input_records = [{"record_id": query_name, "terms": query_terms, "weights": {}}]
    score_records = test_data["disease_records"]

    results = test_data["scorer"].score_records(
        input_records,
        score_records,
        itertools.product(range(len(input_records)), range(len(score_records))),
        threads=1,
    )
    assert len(results) == 8118
    assert round(float(results[0][2]), 2) == 0.04

    [record["weights"].pop("disease_frequency") for record in score_records]
    results = test_data["scorer"].score_records(
        input_records,
        score_records,
        itertools.product(range(len(input_records)), range(len(score_records))),
        threads=1,
    )
    assert len(results) == 8118


def test_no_parents(test_data):
    terms_a = ["HP:0012433", "HP:0000708"]
    terms_b = ["HP:0001249", "HP:0012758"]

    assert (
        list(remove_parents(terms_a, test_data["scorer"].hpo_network))[0]
        == "HP:0012433"
    )
    assert len(remove_parents(terms_b, test_data["scorer"].hpo_network)) == 2


def test_score_self(test_data):
    records = parse_input(
        os.path.join(test_data["parent_dir"], "data/test.score-long.txt"),
        test_data["hpo_network"],
        test_data["alt2prim"],
    )

    input_records = [x for x in records if x["record_id"] in ["213200", "302801"]]

    results = test_data["scorer"].score_records(
        input_records,
        input_records,
        half_product(len(input_records), len(input_records)),
    )
    assert len(results) == 3

    assert round(float(results[1][2]), 2) == 0.1


def test_bmwa(test_data):
    terms_a = ["HP:0001251", "HP:0001263", "HP:0001290", "HP:0004322"]

    terms_b = ["HP:0001263", "HP:0001249", "HP:0001290"]
    weights_a = {"age": [0.67, 1.0, 1.0, 0.4]}
    weights_b = {"age": [1.0, 1.0, 1.0]}

    df = pd.DataFrame(
        [
            [4.22595743e-02, 3.92122308e-02, 3.04851573e-04],
            [1.07473687e-01, 5.05101655e-01, 3.78305515e-04],
            [3.69780479e-04, 3.78305515e-04, 4.64651944e-01],
            [4.17139800e-04, 4.12232546e-04, 3.67984322e-04],
        ],
        index=pd.Index(terms_a, name="a"),
        columns=pd.MultiIndex.from_arrays(
            [["score"] * len(terms_b), terms_b], names=[None, "b"]
        ),
    )

    score_bmwa = test_data["scorer"].best_match_weighted_average(
        df, weights_a, weights_b
    )

    assert round(score_bmwa, 4) == 0.3419

    weights_a = {"age": [1.0] * len(terms_a)}
    score_bmwa = test_data["scorer"].best_match_weighted_average(
        df, weights_a, weights_b
    )
    assert round(score_bmwa, 4) == 0.2985

    weights_a = {"age": [1.0] * len(terms_a)}
    weights_b = {"age": [1.0] * len(terms_b)}
    test_data["scorer"].min_score_mask = None
    score_bmwa = test_data["scorer"].best_match_weighted_average(
        df, weights_a, weights_b
    )
    assert round(score_bmwa, 4) == 0.2985

    terms_a = ["HP:0001251", "HP:0001249", "HP:0001263", "HP:0001290", "HP:0004322"]
    terms_b = ["HP:0001263", "HP:0001249", "HP:0001290"]

    df = pd.DataFrame(
        [
            [4.22595743e-02, 3.92122308e-02, 3.04851573e-04],
            [1.07473687e-01, 5.05101655e-01, 3.78305515e-04],
            [1.07473687e-01, 5.05101655e-01, 3.78305515e-04],
            [3.69780479e-04, 3.78305515e-04, 4.64651944e-01],
            [4.17139800e-04, 4.12232546e-04, 3.67984322e-04],
        ],
        index=pd.Index(terms_a, name="a"),
        columns=pd.MultiIndex.from_arrays(
            [["score"] * len(terms_b), terms_b], names=[None, "b"]
        ),
    )

    weights_a = {"age": [0.67, 0.4, 1.0, 1.0, 0.4]}
    weights_b = {"age": [1.0, 1.0, 1.0]}

    # compute pairwise best match weighted average
    test_data["scorer"].min_score_mask = None
    score_bmwa = test_data["scorer"].best_match_weighted_average(
        df, weights_a, weights_b
    )

    assert round(score_bmwa, 3) == 0.352

    # because both patients were described to have ID, but only patient a
    # has ataxia and so we mask good phenotype matches from being weighted
    # down by default we expect to get a better similarity score
    test_data["scorer"].min_score_mask = 0.05
    score_bmwa = test_data["scorer"].best_match_weighted_average(
        df, weights_a, weights_b
    )

    assert round(score_bmwa, 3) == 0.365


def test_age_weight(test_data):
    # Test age based weight distribution and best_match_weighted_average calculation

    terms_a = [
        "HP:0001251",
        "HP:0001263",
        "HP:0001290",
        "HP:0004322",
    ]  # ATAX, DD, HYP, SS
    terms_b = ["HP:0001263", "HP:0001249", "HP:0001290"]  # DD, ID, HYP

    test_data["hpo_network"] = annotate(
        test_data["hpo_network"],
        test_data["phenotype_to_diseases"],
        test_data["num_diseases_annotated"],
        test_data["alt2prim"],
        ages_distribution_file=test_data["ages_distribution_file"],
    )

    age_a = 9.0
    age_b = 4.0

    # calculate weights based on patients age
    weights_a = {"age": calculate_age_weights(terms_a, age_b, test_data["hpo_network"])}
    weights_b = {"age": calculate_age_weights(terms_b, age_a, test_data["hpo_network"])}

    # make pairwise scores matrix
    df = pd.DataFrame(
        [
            [4.22595743e-02, 3.92122308e-02, 3.04851573e-04],
            [1.07473687e-01, 5.05101655e-01, 3.78305515e-04],
            [3.69780479e-04, 3.78305515e-04, 4.64651944e-01],
            [4.17139800e-04, 4.12232546e-04, 3.67984322e-04],
        ],
        index=pd.Index(terms_a, name="a"),
        columns=pd.MultiIndex.from_arrays(
            [["score"] * len(terms_b), terms_b], names=[None, "b"]
        ),
    )
    # compute pairwise best match weighted average
    score_bmwa = test_data["scorer"].best_match_weighted_average(
        df, weights_a, weights_b
    )

    assert round(float(score_bmwa), 4) == 0.3742

    # set all weights to 1.0, result should be the same as BMA without weights
    weights_a = {"disease_frequency": [1.0] * len(terms_a)}
    weights_b = {"disease_frequency": [1.0] * len(terms_b)}
    score_bmwa = test_data["scorer"].best_match_weighted_average(
        df, weights_a, weights_b
    )

    assert round(float(score_bmwa), 4) == 0.2985

    # test term not in network
    terms_a = ["HP:Not_a_term"]
    weights_a = calculate_age_weights(terms_a, age_b, test_data["hpo_network"])
    assert weights_a == [1.0]

    # term in network no age
    terms_a = ["HP:0000001"]
    weights_a = calculate_age_weights(terms_a, age_b, test_data["hpo_network"])
    assert weights_a == [1.0]


def test_score_pairs_age(test_data):
    # Test reading in records files and calculating pairwise scores
    # read in records
    test_data["hpo_network"] = annotate(
        test_data["hpo_network"],
        test_data["phenotype_to_diseases"],
        test_data["num_diseases_annotated"],
        test_data["alt2prim"],
        ages_distribution_file=test_data["ages_distribution_file"],
    )

    records = parse_input(
        os.path.join(test_data["parent_dir"], "data/test.score-short.txt"),
        test_data["hpo_network"],
        test_data["alt2prim"],
    )

    # create instance the scorer class
    scorer = Scorer(
        test_data["hpo_network"], summarization_method="BMWA", min_score_mask=None
    )

    # select which patients to test in pairwise best_match_weighted_average
    input_records = [x for x in records if x["record_id"] in ["118200", "118210"]]

    results = scorer.score_records(
        input_records,
        input_records,
        [
            (0, 1),
        ],
    )
    assert len(results) == 1

    # the right answer
    answer = np.average(
        [0.017, 0.231, 0.325, 0.0, 0.042, 0.231, 0.325],
        weights=[0.481, 1.0, 1.0, 0.0446, 1.0, 1.0, 1.0],
    )

    assert pytest.approx(float(results[0][2]), 0.01) == answer

    # Test identical records for which one age exist and one doesn't
    input_records = [x for x in records if x["record_id"] in ["118210", "118211"]]

    results = scorer.score_records(
        input_records,
        input_records,
        [
            (0, 1),
        ],
    )
    assert len(results) == 1

    # the right answer
    answer = np.average([0.226, 0.231, 0.325], weights=[0.481, 1.0, 1.0])

    assert pytest.approx(float(results[0][2]), 0.1) == answer


def test_alpha_zero(test_data):
    """the root term should contain all diseases therefore the IC should be zero"""

    root_term_ic = test_data["hpo_network"].nodes["HP:0000118"]["ic"]
    assert root_term_ic == 0.0


def test_leaves_diff_branches_score_zero(test_data):
    """two leaves in different branches
    two leaves therefore beta is zero
    different branches therefore alpha is zero
    define I = (0.0 / (0.0 + 0.0)) as zero and not nan"""
    term_a = "HP:0001290"
    term_b = "HP:0011351"

    score_two_leaves_diff_branches = test_data["scorer"].score_hpo_pair_hrss(
        term_a, term_b
    )
    assert score_two_leaves_diff_branches == 0.0


def test_score_hrss_basic(test_data):
    test_data["scorer"].scoring_method = "HRSS"
    terms_a = ["HP:0001290", "HP:0000118"]
    terms_b = ["HP:0001290", "HP:0011351"]

    assert pytest.approx(0.162, 0.01) == test_data["scorer"].score_term_sets_basic(
        terms_a, terms_b
    )


def test_score_resnik_basic(test_data):
    test_data["scorer"].scoring_method = "Resnik"
    terms_a = ["HP:0001290", "HP:0000118"]
    terms_b = ["HP:0001290", "HP:0011351"]
    assert pytest.approx(1.283, 0.01) == test_data["scorer"].score_term_sets_basic(
        terms_a, terms_b
    )


def test_score_jaccard_basic(test_data):
    test_data["scorer"].scoring_method = "Jaccard"
    terms_a = ["HP:0001290", "HP:0000118"]
    terms_b = ["HP:0001290", "HP:0011351"]

    assert pytest.approx(0.333, 0.01) == test_data["scorer"].score_term_sets_basic(
        terms_a, terms_b
    )


def test_score_word2vec_basic(test_data):
    test_data["scorer"] = Scorer(test_data["hpo_network"], scoring_method="word2vec")
    terms_a = ["HP:0001290", "HP:0000118"]
    terms_b = ["HP:0001290", "HP:0011351"]

    assert pytest.approx(
        test_data["scorer"].score_term_sets_basic(terms_a, terms_b), 0.01
    ) == pytest.approx(0.156, 0.01)


def test_score_word2vec_out_of_vocab(test_data):
    test_data["scorer"] = Scorer(test_data["hpo_network"], scoring_method="word2vec")
    terms_a = ["HP:NOT_A_TERM", "HP:0000118"]
    terms_b = ["HP:0001290", "NOT_A_TERM"]

    assert pytest.approx(
        test_data["scorer"].score_term_sets_basic(terms_a, terms_b), 0.01
    ) == pytest.approx(0.063, 0.01)


def test_score_word2vec_empty(test_data):
    test_data["scorer"] = Scorer(test_data["hpo_network"], scoring_method="word2vec")
    terms_a = []
    terms_b = ["HP:0001290", "HP:0011351"]

    assert test_data["scorer"].score_term_sets_basic(terms_a, terms_b) == 0.0
