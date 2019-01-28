import sys

import pandas as pd

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score


def clustering_assign(X, linkage, k, samples):
    """cluster and provide class labels for every sample in the record file
    :param X: The array of similarity scores
    :type X: array
    :param linakge: the linkage type to use {single, average, complete}
    :param k: The number of clusters to define
    :param samples: The array of sample ids for each record in X
    :type samples: array
    """
    clusterer = AgglomerativeClustering(n_clusters=k, linkage=linkage, affinity='precomputed')
    cluster_labels = clusterer.fit_predict(X)

    for sample, label in zip(samples, cluster_labels):
        try:
            sys.stdout.write('\t'.join([
                f'{sample}',
                f'{label}',
            ]))
            sys.stdout.write('\n')
        finally:
            sys.stdout.flush()


def clustering_grid_search(X, linkage, k):
    """cluster and provide metrics on a phenosim result file.
    """
    clusterer = AgglomerativeClustering(n_clusters=k, linkage=linkage, affinity='precomputed')
    cluster_labels = clusterer.fit_predict(X)
    silhouette_avg = silhouette_score(X, cluster_labels, metric='precomputed').round(4)
    cluster_sizes = pd.Series(cluster_labels).sort_index().value_counts().tolist()

    try:
        sys.stdout.write('\t'.join([
            f'{k}',
            f'{linkage}',
            f'{silhouette_avg}',
            f'{cluster_sizes}',
        ]))
        sys.stdout.write('\n')
    finally:
        sys.stdout.flush()
