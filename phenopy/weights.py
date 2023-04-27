import logging

import networkx as nx
import pandas as pd
import sys

from functools import lru_cache
from scipy.stats import truncnorm
from typing import List


@lru_cache(maxsize=1300000)
def hpo_age_to_weight(
        hpo_network: nx.MultiGraph, term_id: str, age: int) -> float:
    """
    calculate weight based on truncated normal distribution CDF
    """
    if term_id not in hpo_network.nodes or age is None:
        return 1.0
    elif 'age_dist' in hpo_network.nodes[term_id]:
        return hpo_network.nodes[term_id]['age_dist'].cdf(float(age))
    else:
        return 1.0


def calculate_age_weights(
        terms: List,
        age: int,
        hpo_network:
        nx.MultiGraph) -> List[float]:
    """
    Calculates an age-based weight vector given an iterable of terms.
    """
    weights = []
    for term_id in terms:
        weights.append(hpo_age_to_weight(hpo_network, term_id, age))

    return weights


def get_truncated_normal(mean: float = 0.0, sd: float = 1.0,
                         low: float = 0.0, upp: float = 10.0) -> truncnorm:
    """
    Model truncated normal given summary stats
    """
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)


def make_age_distributions(
        phenotype_age_file: str,
        logger: logging.Logger = None) -> pd.DataFrame:
    """
    Read in phenotype ages file and convert to pandas object with modeled distributions
    """

    try:
        df = pd.read_csv(phenotype_age_file, sep='\t', names=['hpid', 'mean', 'std'])

    except (FileNotFoundError, PermissionError) as e:

        if logger is not None:
            logger.critical(e)
        else:
            sys.stderr.write(str(e))
        exit(1)

    distributions = []
    for rec in df.to_dict('records'):

        try:
            # model truncated normal
            dist = get_truncated_normal(
                mean=rec['mean'],
                sd=rec['std'],
                low=0,
                upp=rec['mean']
            )
            distributions.append({'hpid': rec['hpid'], 'age_dist': dist})

        except ValueError as e:
            if logger is not None:
                logger.critical(e)
            else:
                sys.stderr.write(str(e))
            exit(1)

    return pd.DataFrame.from_dict(distributions).set_index('hpid')
