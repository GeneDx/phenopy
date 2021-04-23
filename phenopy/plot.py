
import matplotlib.pyplot as plt
import numpy as np
from phenopy.cluster import extract_tfidf_features, hpo2name

from collections import Counter


def plot_basic_dbscan(X, core_samples_mask, labels):
    plt.figure(figsize=(15, 12))
    unique_labels = list(set(labels))
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


def plot_tfidf_features(result, labels, feature_to_hps, plot_noise=True):
    plt.figure(figsize=(15, 5))
    fig, axes = plt.subplots(1,1, figsize=(15, 5))
    original_df = result.copy()
    original_df['cluster_id'] = labels
    cut_at = 0.01
    cluster_info = []
    tfidf_features_df = extract_tfidf_features(result, labels)
    label_top_features = 2

    tfidf_features_df['index'] = tfidf_features_df.index
    clusters = sorted(tfidf_features_df['cluster_id'].unique())
    colors = [plt.cm.RdBu_r(each)
              for each in np.linspace(0, 1, len(clusters))]
    x_labels_pos = []
    if not plot_noise:
        clusters = [x for x in clusters if x >= 0]
    for k, col in zip(clusters, colors):
        if k < 0:
            # Black used for noise.
            col = [0, 0, 0, 1]

        df_sub = tfidf_features_df[tfidf_features_df['cluster_id']==k]
        df_sub = df_sub.reset_index(drop=True)
        xx = df_sub['index']
        yy = df_sub['norm_tfidf']
        n_size = (df_sub['n'].values / (np.ones((1, yy.shape[0])) * yy.shape[0])) * 100
        x_labels_pos.append((xx.max() - (xx.max() - xx.min())/2))
        plt.scatter(xx,yy, marker='o', color=col, s=n_size)
        if label_top_features > 0:
            sigs = df_sub.sort_values("norm_tfidf", ascending=False).head(label_top_features)
        else:
            sigs = tfidf_features_df[(tfidf_features_df['cluster_id']==k) & (tfidf_features_df['norm_tfidf'] > cut_at)]
        if not sigs.empty:
            all_terms = original_df.query(f"cluster_id == {k}")[['terms']].explode('terms')['terms'].to_list()
            for x, y, feat_id in sigs[['index','norm_tfidf','hp_features']].values:
                feat_id = int(feat_id)
                feat_terms = feature_to_hps[feat_id]
                rel_terms = [x for x in all_terms if x in feat_terms]
                if rel_terms:
                    mc_hpid = Counter(rel_terms).most_common(1)[0][0]
                    text = hpo2name(mc_hpid)
                else:
                    text = 'None'
                if not text:
                    text = "No text"
                print(f"{k} {feat_id} {mc_hpid} {text}, TF-IDF:{y}, index:{x}")
                cluster_info.append(dict(cluster=k,hpid=mc_hpid,name=text,tfidf=y))
                plt.scatter(x,y, marker='o', color=col, edgecolor='k', s=100)
                plt.text(x, y+0.01, text.split('(')[0], fontsize=20, color='grey', rotation=90)

    plt.xticks(x_labels_pos, labels=clusters)
    if label_top_features == 0:
        plt.axhline(y=cut_at,color="grey", linestyle='--')
    # make xaxis invisible
    #axes.xaxis.set_visible(False)
    # make spines (the box) invisible
    plt.setp(axes.spines.values(), visible=True)
    # remove ticks and labels for the left axis
    #axes.tick_params(left=False, labelleft=False)
    #remove background patch (only needed for non-white background)
    axes.patch.set_visible(False)
    plt.xlabel('Cluster Num.', size=20, rotation=0)
    plt.ylabel('TF-IDF (normalized)', size=20)
    plt.xticks(rotation=0, size=20)
    plt.yticks(rotation=0, size=10)
    return plt