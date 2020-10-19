import os
import umap

from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.feature_selection import SelectKBest, chi2
from sklearn.preprocessing import LabelBinarizer

from collections import Counter

import pandas as pd
import numpy as np


from phenopy.config import logger
from phenopy import generate_annotated_hpo_network

phenopy_data_directory = os.path.join(os.getenv('HOME'), '.phenopy/data')
obo_file = os.path.join(phenopy_data_directory, 'hp.obo')
disease_to_phenotype_file = os.path.join(phenopy_data_directory, 'phenotype.hpoa')
hpo_network, alt2prim, disease_records = generate_annotated_hpo_network(obo_file, \
                                                                        disease_to_phenotype_file,
                                                                        ages_distribution_file=None)


def hpo2name(term):
    if term in hpo_network.nodes():
        return hpo_network.nodes()[term]['name']
    else:
        return np.nan


def hpo2depth(term):
    if term in hpo_network.nodes():
        return hpo_network.nodes()[term]['depth']
    else:
        return np.nan


def hpo2ic(term):
    if term in hpo_network.nodes():
        return hpo_network.nodes()[term]['ic']
    else:
        return np.nan


#dipr_tfidf = pd.read_pickle("dipr_tfidf.pkl")


# def hpo2tfidf(term, hp_to_feature):
#     min_val = dipr_tfidf['norm_tfidf'].min()
#     feature = hp_to_feature[term]
#     max_tfidf = dipr_tfidf.query(f"hp_feature == {feature}")["norm_tfidf"].max()
#     if max_tfidf >= min_val:
#         return max_tfidf
#     else:
#         return min_val


def hp_to_weight(term, hp_to_feature, method="ic"):
    if method == "ic":
        return hpo2ic(term)
    elif method == "depth":
        return hpo2depth(term)
    elif method == "tfidf":
        raise NotImplementedError
        #TODO #hpo2tfidf(term, hp_to_feature)
    else:
        return 1


def process_kfile(kfile, k=1000):
    try:
        kdf = pd.read_csv(kfile, sep="\t")
    except FileNotFoundError:
        logger.critical(
            'Phenotype groups file not found')
        exit(1)
    if k == 1000:
        feature_to_hps = dict(zip(kdf['phenotype_group_k1000'],kdf['HPO_id']))
        hp_to_feature = dict(zip(kdf['HPO_id'],kdf['phenotype_group_k1000']))
    elif k == 1500:
        feature_to_hps = dict(zip(kdf['phenotype_group_k1500'], kdf['HPO_id']))
        hp_to_feature = dict(zip(kdf['HPO_id'], kdf['phenotype_group_k1500']))

    n_features = len(feature_to_hps.keys())

    return feature_to_hps, hp_to_feature, n_features


def prep_cluster_data(df, hp_to_feature, weights="default"):

    df['hp_features'] = df['terms'].apply(lambda hpos: [hp_to_feature[hp] for hp in hpos])
    df['hp_weights'] = df['terms'].apply(lambda hpos: [hp_to_weight(hp, hp_to_feature, method=weights) for hp in hpos])
    df['hp_feature_counts'] = df['hp_features'].apply(Counter)

    return df


def feature_array(cntr, weights=None, n_features=None, multiplier=1.0, superweighted_features=None):

    f_array = [0] * n_features
    if not weights:
        weights = np.ones(len(cntr))
    if not superweighted_features:
        superweighted_features = []
    i = 0
    for feature_index, count in cntr.items():
        if feature_index in superweighted_features:
            f_array[feature_index] = count * weights[i] * multiplier
        else:
            f_array[feature_index] = count * weights[i]
        i += 1
    return f_array


def prep_feature_array(res_df, n_features, superweighted_features=None, multiplier=1.0, scaler=MinMaxScaler()):
    if not superweighted_features:
        superweighted_features = []
    X = res_df.apply(lambda x:
                    feature_array(x['hp_feature_counts'],
                                        weights=x['hp_weights'],
                                        n_features=n_features,
                                        multiplier=multiplier,
                                        superweighted_features=superweighted_features),
                    axis=1).tolist()

    return scaler.fit_transform(np.array(X))


def compute_dx_yield(res_df):
    dx_yield = res_df.groupby('cluster_id')['is_diagnosed'].value_counts().unstack().reset_index()
    dx_yield.columns = ['cluster_id', 'False', 'True']
    dx_yield.fillna(0, inplace=True)
    dx_yield['dx_yield'] = (dx_yield['True'] / (dx_yield['True'] + dx_yield['False'])) * 100
    dx_yield['N'] = dx_yield['True'] + dx_yield['False']
    dx_yield = dx_yield.sort_values(by=['dx_yield'], ascending=False).reset_index()
    return dx_yield


def compute_ages(res_df):
    # Median age per cluster in months
    ages = res_df.groupby("cluster_id")['age'].median().reset_index()
    return ages


def compute_genes(res_df):
    genes = res_df.query("is_diagnosed == 1")
    genes['gene_name'] = genes['gene_name'].apply(lambda x: [x])
    genes = genes.groupby('cluster_id')['gene_name'].sum().reset_index()

    # Most common genes per cluster
    genes['gene'] = genes.apply(lambda x: str(Counter(x['gene_name']).most_common()), axis=1)
    return genes


def apply_umap(X_vects, n_neighbors=30, n_components=2, min_dist=0.01, metric='euclidean'):
    xa_sub_vectors_embedded = umap.UMAP(n_neighbors=n_neighbors,
                                        n_components=n_components,
                                        min_dist=min_dist,
                                        random_state=123,
                                        metric=metric,
                                        ).fit_transform(X_vects)
    #xa_sub_vectors_embedded_df = pd.DataFrame(xa_sub_vectors_embedded, columns=['umap1', 'umap2'])
    return xa_sub_vectors_embedded


def tfidf(term_i, cluster_j, df):
    # n clusters
    N = df['cluster_id'].nunique()
    # n times term in cluster_j
    try:
        tf_cl_j = df[(df['cluster_id'] == cluster_j) & (df['hp_features'] == term_i)]['n'].sum()
    except:
        print(term_i, cluster_j)
    # n of docs with i
    df_i = df[df['hp_features'] == term_i]['cluster_id'].nunique()
    # n of term_i in cluster_j
    w = tf_cl_j * np.log(N / df_i)
    return w

#X = xa_sub_vectors_embedded_df[['umap1', 'umap2']].values


def dbscan(X, eps=0.40, min_samples=10):
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    stats = dict(
        n_clusters=n_clusters_,
        n_noise=n_noise_,
        silhouette_score=metrics.silhouette_score(X, labels)
    )

    return labels, core_samples_mask, stats


def process_dbscan_results(res_df, X_vects):

    selector = SelectKBest(chi2, k='all')
    cluster_labels = LabelBinarizer().fit_transform(res_df['cluster_id'])
    selector.fit(X_vects, cluster_labels)

    # Add here Merge features first by hrss before feature selection

    p_values = selector.pvalues_
    scores = -np.log10(selector.pvalues_)

    top_k_feature_indices = selector.get_support()

    indx_features = np.where(top_k_feature_indices == True)

    df_pvals = pd.DataFrame(p_values, columns=['pval']).reset_index()
    df_pvals['-log(p_val)'] = scores
    df_pvals['pval'].fillna(p_values[~np.isnan(p_values)].max(), inplace=True)

    max_score = scores[(~np.isinf(scores)) & (~np.isnan(scores))].max()
    min_score = scores[(~np.isinf(scores)) & (~np.isnan(scores))].min()

    df_pvals['-log(p_val)'].replace(np.inf, max_score, inplace=True)
    df_pvals['-log(p_val)'].fillna(min_score, inplace=True)

    cl_feat = res_df.groupby('cluster_id')['hp_features'].sum().reset_index()
    cl_feat = cl_feat.explode('hp_features')
    cl_feat = pd.DataFrame(cl_feat.groupby(['cluster_id', 'hp_features'])['hp_features'].count()).rename(
        columns={'hp_features': 'n'}).reset_index()

    # df_sig_phenos = xa_cl_feat.query("pval < 0.05")
    cl_feat['tfidf'] = cl_feat.apply(lambda x: tfidf(x['hp_features'], x['cluster_id'], cl_feat.copy()), axis=1)

    return cl_feat