from scipy.stats import truncnorm


def get_truncated_normal(mean=0.0, sd=1.0, low=0.0, upp=10.0):
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)


def age_to_weights(age_dist, age):
    """calculate weight based on truncated normal distribution CDF"""
    if age > age_dist.interval(0.99)[1]:
        return 1.0
    else:
        return age_dist.cdf(age)
