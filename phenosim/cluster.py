import sys

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score


def clustering_assign(X, method, samples, k, **kwargs):
    """cluster and provide class labels for every sample in the record file
    :param X: The array of similarity scores
    :type X: array
    :param samples: The array of sample ids for each record in X
    :type samples: array
    :param k: The number of clusters to define
    :param kwargs: Keyword arguements to pass to scikit-learn or pyclustering
    """
    clusterer = AgglomerativeClustering(n_clusters=k, linkage=method, affinity='precomputed')
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


def clustering_grid_search(X, method, k):
    """cluster and provide metrics on a phenosim result file.
    """
    clusterer = AgglomerativeClustering(n_clusters=k, linkage=method, affinity='precomputed')
    cluster_labels = clusterer.fit_predict(X)
    silhouette_avg = silhouette_score(X, cluster_labels, metric='precomputed').round(4)

    try:
        sys.stdout.write('\t'.join([
            f'{k}',
            f'{method}',
            f'{silhouette_avg}',
        ]))
        sys.stdout.write('\n')
    finally:
        sys.stdout.flush()
