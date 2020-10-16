
import matplotlib.pyplot as plt
import numpy as np


def plot_basic_dbscan(X, core_samples_mask, labels):
    plt.figure(figsize=(15, 12))
    unique_labels = set(labels)
    colors = [plt.cm.RdBu_r(each)
              for each in np.linspace(0, 1, len(unique_labels))]
    for k, col in zip(unique_labels, colors):
        k_label = f"Cluster #{k} "
        if k < 0:
            # Black used for noise.
            col = [0, 0, 0, 1]

        class_member_mask = (labels == k)

        xy = X[class_member_mask & core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                 markeredgecolor='k', markersize=20, label=k_label)

        xy = X[class_member_mask & ~core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                 markeredgecolor='k', markersize=10)

    plt.xlabel("UMAP1", fontsize=20)
    plt.ylabel("UMAP2", fontsize=20)
    return plt