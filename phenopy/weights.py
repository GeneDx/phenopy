import logging

import networkx as nx
import pandas as pd
import sys

from functools import lru_cache
from typing import List
import numpy as np


@lru_cache(maxsize=1300000)
def hpo_age_to_weight(hpo_network: nx.MultiGraph, term_id: str, age: int) -> float:
    """
    calculate weight based on truncated normal distribution CDF
    """
    if term_id not in hpo_network.nodes or age is None:
        return 1.0
    elif "age_dist" in hpo_network.nodes[term_id]:
        return hpo_network.nodes[term_id]["age_dist"].cdf(float(age))
    else:
        return 1.0


def calculate_age_weights(
    terms: List, age: int, hpo_network: nx.MultiGraph
) -> List[float]:
    """
    Calculates an age-based weight vector given an iterable of terms.
    """
    weights = []
    for term_id in terms:
        weights.append(hpo_age_to_weight(hpo_network, term_id, age))

    return weights


def get_truncated_normal(
    mean: float, sd: float, lower: float, upper: float, instances: int = 1000000
) -> np.ndarray:
    """
    Simulates a truncated normal distribution
    """
    # Create the normal distribution
    distribution = np.random.normal(mean, sd, instances)

    # Truncate all values outside of the range
    distribution = np.array([i for i in distribution if lower <= i <= upper])

    return distribution


def get_empirical_cdf(value: float, distribution: np.ndarray) -> float:
    """
    Calculates the empirical cumulative distribution function for a given value within
    a given distribution.
    """
    # Sort the distribution
    data_sorted = np.sort(distribution)

    # Determine the CDF for the values within the distribution
    cdf = np.linspace(0, 1, len(distribution))

    # Establish as a dataframe
    df = pd.DataFrame(list(zip(data_sorted, cdf)), columns=["value", "cdf"])

    # Return the maximum CDF value for the given value
    return df[df["value"] <= value]["cdf"].max()


def make_age_distributions(
    phenotype_age_file: str, logger: logging.Logger = None
) -> pd.DataFrame:
    """
    Read in phenotype ages file and convert to pandas object with modeled distributions
    """

    try:
        df = pd.read_csv(phenotype_age_file, sep="\t", names=["hpid", "mean", "std"])

    except (FileNotFoundError, PermissionError) as e:

        if logger is not None:
            logger.critical(e)
        else:
            sys.stderr.write(str(e))
        exit(1)

    distributions = []
    for rec in df.to_dict("records"):

        try:
            # model truncated normal
            dist = get_truncated_normal(
                mean=rec["mean"], sd=rec["std"], lower=0, upper=rec["mean"]
            )
            distributions.append({"hpid": rec["hpid"], "age_dist": dist})

        except ValueError as e:
            if logger is not None:
                logger.critical(e)
            else:
                sys.stderr.write(str(e))
            exit(1)

    return pd.DataFrame.from_dict(distributions).set_index("hpid")
