import pandas as pd
import sys

from scipy.stats import truncnorm


def get_truncated_normal(mean=0.0, sd=1.0, low=0.0, upp=10.0):
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)


def age_to_weights(age_dist, age):
    """calculate weight based on truncated normal distribution CDF"""
    if age is None:
        return 1.0
    elif age > age_dist.interval(0.99)[1]:
        return 1.0
    else:
        return age_dist.cdf(age)


def make_age_distributions(phenotype_age_file, logger=None):
    """
    Read in phenotype ages file and convert to pandas object with modeled distributions
    :param phenotype_age_file:
    :param logger:
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

    distros = []

    for rec in df.to_dict('records'):

        try:
            X = get_truncated_normal(mean=rec['mean'], sd=rec['std'], low=0, upp=rec['mean'])
            distros.append({'hpid': rec['hpid'], 'age_dist': X})
        except ValueError as e:
            if logger is not None:
                logger.critical(e)
            else:
                sys.stderr.write(str(e))
            exit(1)

    return pd.DataFrame.from_dict(distros).set_index('hpid')

