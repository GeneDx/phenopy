import os
import pytest
import pandas as pd
from phenopy.config import logger
from phenopy.weights import (
    get_truncated_normal,
    hpo_age_to_weight,
    make_age_distributions,
)


def test_age_weights(test_data):
    assert hpo_age_to_weight(test_data["hpo_network"], "HP:0001251", 9.0) == 1.0
    assert (
        pytest.approx(
            hpo_age_to_weight(test_data["hpo_network"], "HP:0001251", 5.0), 0.01
        )
        == 1.0
    )


def test_make_age_distributions(test_data):
    with pytest.raises(SystemExit) as se:
        make_age_distributions("notafilepath/notafile")

    assert se.type == SystemExit
    assert se.value.code == 1

    with pytest.raises(SystemExit) as se:
        make_age_distributions("notafilepath/notafile", logger=logger)

    assert se.type == SystemExit
    assert se.value.code == 1

    ages_truth = pd.DataFrame(
        [
            {
                "hpid": "HP:0001251",
                "age_dist": get_truncated_normal(6.0, 3.0, 0.0, 6.0),
            },
            {
                "hpid": "HP:0001263",
                "age_dist": get_truncated_normal(1.0, 1.0, 0.0, 1.0),
            },
            {
                "hpid": "HP:0001290",
                "age_dist": get_truncated_normal(1.0, 1.0, 0.0, 1.0),
            },
            {
                "hpid": "HP:0004322",
                "age_dist": get_truncated_normal(10.0, 3.0, 0.0, 10.0),
            },
            {
                "hpid": "HP:0001249",
                "age_dist": get_truncated_normal(6.0, 3.0, 0.0, 6.0),
            },
        ]
    ).set_index("hpid")

    phenotype_ages_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "data/phenotype_age.tsv"
    )
    df = make_age_distributions(phenotype_ages_file)
    assert set(ages_truth.index) == set(df.index)

    for hpid in ages_truth.index:
        assert pytest.approx(
            ages_truth.loc[hpid]["age_dist"].mean(), 0.001
        ) == pytest.approx(df.loc[hpid]["age_dist"].mean(), 0.001)


def test_get_truncated_normal(test_data):
    assert pytest.approx(get_truncated_normal(6.0, 1.0, 0.0, 6.0).mean(), 0.01) == 5.20
    assert (
        pytest.approx(get_truncated_normal(6.0, 1.0, 0.0, 6.0).cdf(3.0), 0.01) == 0.0027
    )
    assert (
        pytest.approx(get_truncated_normal(6.0, 1.0, 0.0, 6.0).cdf(12.0), 0.01) == 1.0
    )
