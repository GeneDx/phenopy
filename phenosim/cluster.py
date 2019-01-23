import sys

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score


def cluster_phenosim(X, method, k):
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
