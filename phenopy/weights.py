import pandas as pd
import sys

from scipy.stats import truncnorm


def get_truncated_normal(mean=0.0, sd=1.0, low=0.0, upp=10.0):
    """
    Model truncated normal given summary stats
    :param mean: mean
    :param sd: standard deviation
    :param low: lower boundary
    :param upp: upper boundary
    :return: distribution
    """
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)


def age_to_weights(age_dist, age):
    """calculate weight based on truncated normal distribution CDF
    :param age_dist: distribution of ages for phenotype
    :param age: age of patient
    :return: cdf
    """

    if age is None:
        return 1.0
    else:
        return age_dist.cdf(age)


def make_age_distributions(phenotype_age_file, logger=None):
    """
    Read in phenotype ages file and convert to pandas object with modeled distributions
    :param phenotype_age_file: path to tab file containing hpid, mean phenotype age, standard deviation
    :param logger: to log
    :return: pandas df
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
            dist = get_truncated_normal(mean=rec['mean'], sd=rec['std'], low=0, upp=rec['mean'])
            distributions.append({'hpid': rec['hpid'], 'age_dist': dist})

        except ValueError as e:
            if logger is not None:
                logger.critical(e)
            else:
                sys.stderr.write(str(e))
            exit(1)

    return pd.DataFrame.from_dict(distributions).set_index('hpid')

