import os
from umap import UMAP
from sklearn.preprocessing import MinMaxScaler
from sklearn import metrics

import pandas as pd
import numpy as np
from phenopy.config import logger
from phenopy.build_hpo import generate_annotated_hpo_network

import matplotlib.pyplot as plt
from itertools import product
from phenopy.score import Scorer
from collections import Counter
from sklearn.cluster import DBSCAN
from phenopy.config import config
from sklearn.feature_selection import SelectKBest, chi2
from sklearn.preprocessing import LabelBinarizer
from multiprocessing import Pool
import hdbscan

phenopy_data_directory = os.path.join(os.getenv('HOME'), '.phenopy/data')
obo_file = os.path.join(phenopy_data_directory, 'hp.obo')
disease_to_phenotype_file = os.path.join(phenopy_data_directory, 'phenotype.hpoa')
hpo_network, alt2prim, disease_records = generate_annotated_hpo_network(obo_file,
                                                                        disease_to_phenotype_file,
                                                                        ages_distribution_file=None)


plt.rcParams["figure.figsize"] = (10, 10)


class Cluster:
    def __init__(self, data, scoring_method='HRSS', use_pdr=True, kfile=None, k=1000):

        terms_columns = {'hpo_terms', 'terms', 'txt2hpo'}
        terms_column = set(data.columns).intersection(terms_columns)
        self.terms_column = list(terms_column)[0]

        id_columns = {'record_id', 'id', 'accession_num', 'case_id'}
        id_column = set(data.columns).intersection(id_columns)
        self.id_column = list(id_column)[0]
        self.use_pdr = use_pdr
        self.k = k
        self.scoring_method = scoring_method
        self.alt2prim = alt2prim
        self.data = self.prep_data(data)
        self.scorer = Scorer(hpo_network, scoring_method=scoring_method)
        self.n_cpus = os.cpu_count()
        self.feature_to_hps = {}
        self.hp_to_feature = {}
        self.n_components = None
        self.n_features = None
        self.feature_array = None
        self.tfidf_features = None
        self.pairwise_map = None

        if kfile is None:
            self.kfile = config['phenotype_groups']['phenotype_groups_file']
        else:
            self.kfile = kfile

    def prep_data(self, df_input):
        df = df_input.copy()
        if df.empty:
            print("Empty df")
            return None
        if not self.id_column:
            print("No id column found")
            return None
        if not self.terms_column:
            print("No HPO terms found")
            return None
        if isinstance(df.iloc[0, :][self.terms_column], str):
            sep_by = None
            for sep_by in [',', '|', '\t', ';']:
                if df[self.terms_column].str.contains(",").unique()[0]:
                    break
            if not sep_by:
                print("Could not find HPO term sep")
                return None
            df[self.terms_column] = df[self.terms_column].str.split(sep_by)
        # convert alternate ids to primary
        df_ = df.explode(self.terms_column)
        df_ = df_.merge(
            pd.DataFrame(zip(self.alt2prim.keys(), self.alt2prim.values()), columns=[self.terms_column, 'alt']),
            how="left")
        df_['alt'] = df_["alt"].fillna(df_[self.terms_column])
        df_ = df_.drop(columns=[self.terms_column]).rename(columns={"alt": self.terms_column}).groupby(self.id_column)[
            self.terms_column].agg(list).reset_index()
        df = df.drop(columns=[self.terms_column]).merge(df_, on=self.id_column)
        return df

    def process_kfile(self):
        try:
            kdf = pd.read_table(self.kfile)
        except FileNotFoundError:
            logger.critical(
                'Phenotype groups file not found')
            exit(1)

        if self.k in [1000, 1500]:
            df_k = kdf.groupby(f"phenotype_group_k{self.k}")["HPO_id"].apply(list).reset_index()
            self.feature_to_hps = dict(zip(df_k[f"phenotype_group_k{self.k}"], df_k["HPO_id"]))
            self.hp_to_feature = dict(zip(kdf['HPO_id'], kdf[f"phenotype_group_k{self.k}"]))
        else:
            print(f"No features computed for k={self.k}")

        self.n_features = len(self.feature_to_hps.keys())

    def prep_cluster_data(self, weights="default"):
        self.data['hp_features'] = self.data[self.terms_column].apply(
            lambda hpos: [self.hp_to_feature[hp] for hp in hpos])
        self.data['hp_weights'] = self.data[self.terms_column].apply(
            lambda hpos: [self.hp_to_weight(hp, method=weights) for hp in hpos])
        self.data['hp_feature_counts'] = self.data['hp_features'].apply(Counter)

    def make_feature_array(self, cntr, weights=None, n_features=None, multiplier=1.0, superweighted_features=None):
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

    def prep_feature_array(self, superweighted_features=None, multiplier=1.0, scaler=MinMaxScaler()):
        if not superweighted_features:
            superweighted_features = []
        X = self.data.apply(lambda x:
                            self.make_feature_array(x['hp_feature_counts'],
                                                    weights=x['hp_weights'],
                                                    n_features=self.n_features,
                                                    multiplier=multiplier,
                                                    superweighted_features=superweighted_features),
                            axis=1).tolist()
        self.feature_array = scaler.fit_transform(np.array(X))

    def hp_to_weight(self, term, method="ic"):
        if method == "ic":
            return self.hpo2ic(term)
        elif method == "depth":
            return self.hpo2depth(term)
        elif method == "tfidf":
            raise NotImplementedError
            # TODO #hpo2tfidf(term, hp_to_feature)
        else:
            return 1

    @staticmethod
    def hpo2name(self, term):
        if term in hpo_network.nodes():
            return hpo_network.nodes()[term]['name']
        else:
            return np.nan

    @staticmethod
    def hpo2depth(self, term):
        if term in hpo_network.nodes():
            return hpo_network.nodes()[term]['depth']
        else:
            return np.nan

    @staticmethod
    def hpo2ic(self, term):
        if term in hpo_network.nodes():
            return hpo_network.nodes()[term]['ic']
        else:
            return np.nan

    def sim_bma(self, a, b):
        a = list(set(a))
        b = list(set(b))

        results = self.scorer.score_term_sets_basic(
            a, b
        )
        if results:
            return results
        else:
            return 0.0

    def bma_sim(self, p):
        return [p[0], p[1], self.sim_bma(self.pairwise_map[p[0]], self.pairwise_map[p[1]])]

    def compute_sim(self):
        df_input_cp = self.data.copy()
        self.pairwise_map = dict(zip(df_input_cp[self.id_column], df_input_cp[self.terms_column]))
        computed = []
        case_product = product(self.pairwise_map.keys(), self.pairwise_map.keys())
        pool = Pool(self.n_cpus)
        for res in pool.map(self.bma_sim, case_product):
            computed.append(res)
        pool.close()

        df_computed = pd.DataFrame(data=computed, columns=['case1', 'case2', 'sim'])
        df_computed = df_computed.sort_values(['case1', 'case2'])
        n_unique_names = df_computed['case2'].nunique()
        names = df_computed['case2'].unique()
        df_computed.drop_duplicates(inplace=True)
        values = np.array_split(df_computed['sim'].values, n_unique_names)
        values = np.array(values, dtype="float32")
        # values = MinMaxScaler().fit_transform(values)
        values[values > 1.0] = 1.0
        values = pd.DataFrame(values)
        values.columns = names
        values.index = names
        return values

    def umap(self, metric=None, n_neighbors=15, min_dist=0.0, n_components=2):
        self.n_components = n_components
        # is 1 - sim a valid distance metric
        self.dimention_labels = [f'umap{x}' for x in range(1, n_components + 1)]
        if self.use_pdr:
            self.process_kfile()
            self.prep_cluster_data()
            self.prep_feature_array()
            if metric is None:
                # check for semantic dimentionality reduction
                if len(self.feature_to_hps.keys()) != len(self.feature_to_hps.values()):
                    metric = 'braycurtis'
                else:
                    metric = 'euclidean'
            umap = UMAP(
                n_neighbors=n_neighbors,
                metric=metric,
                min_dist=min_dist,
                random_state=42,
                n_components=n_components
            )
            X = umap.fit_transform(self.feature_array)
            self.data[self.dimention_labels] = X[:, :]
        else:
            values = self.compute_sim()
            X = UMAP(metric="precomputed").fit_transform(1 - values.values)
            df_umap = pd.DataFrame(X, index=values.columns.tolist(), columns=self.dimention_labels)
            df_umap[self.id_column] = df_umap.index
            self.data = self.data.drop(self.dimention_labels, axis=1, errors='ignore')
            self.data = self.data.merge(df_umap, on=self.id_column)
            # for d in enumerate(self.dimention_labels):

    def hdbscan(self, min_samples=9, min_cluster_size=10):
        clusterer = hdbscan.HDBSCAN(metric='euclidean', min_samples=min_samples, min_cluster_size=min_cluster_size)
        clusterer.fit(self.data[[f"umap{x}" for x in range(1, self.n_components + 1)]].values)
        self.data['cluster_id'] = clusterer.labels_.tolist()

    def dbscan(self, eps=0.40, min_samples=10):
        db = DBSCAN(eps=eps, min_samples=min_samples).fit(self.data[['umap1', 'umap2']].values)
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_
        # Number of clusters in labels, ignoring noise if present.
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise = list(labels).count(-1)
        if n_clusters > 1:
            silhouette_score = metrics.silhouette_score(self.data[['umap1', 'umap2']].values, labels)
        else:
            silhouette_score = np.nan

        self.dbscan_stats = dict(
            n_clusters=n_clusters,
            n_noise=n_noise,
            silhouette_score=silhouette_score,
        )
        self.dbscan_core_samples_mask = core_samples_mask
        self.data['cluster_id'] = labels

    def precomputed_umap(self):
        self.umap(self.compute_sim(), metric='precomputed')
        self.hdbscan()

    def features_umap(self, metric, k=1000):
        self.process_kfile(k=k)
        self.prep_cluster_data()
        self.prep_feature_array()
        self.umap(self.feature_array, metric=metric)
        self.hdbscan()

    def visualize(self, title='', color_by=None, annotate_with=None, show_quads=False):
        if 'umap1' not in self.data.columns and 'umap2' not in self.data.columns:
            print("At least 2 dimentions expected")
            return None

        if isinstance(color_by, str):
            if color_by not in self.data.columns:
                print(f"No column named '{color_by}' in df")
                return None
        elif isinstance(color_by, (list, set)):
            for cb in color_by:
                if cb not in self.data.columns:
                    print(f"No column named '{cb}' in df")
                    return None

        fig, axs = plt.subplots()

        if color_by is not None:
            for name, grp in self.data.groupby(color_by):
                axs.scatter(x=grp["umap1"], y=grp["umap2"], label=name, s=120)
        else:
            axs.scatter(x=self.data["umap1"], y=self.data["umap2"], s=120)

        if annotate_with and annotate_with in self.data.columns:
            for x, y, label in self.data[["umap1", "umap2", annotate_with]].values:
                axs.text(x, y, label, fontsize=10)

        if show_quads:
            umap1_range = (self.data['umap1'].min(), self.data['umap1'].max())
            umap2_range = (self.data['umap2'].min(), self.data['umap2'].max())
            umap1_half = umap1_range[0] + (umap1_range[1] - umap1_range[0]) / 2
            umap2_half = umap2_range[0] + (umap2_range[1] - umap2_range[0]) / 2
            umap1_offset = (umap1_range[1] - umap1_half) / 50
            umap2_offset = (umap2_range[1] - umap2_half) / 50

            axs.text(umap1_range[0] + umap1_offset, umap2_range[1] - umap2_offset, 'A', fontsize=30)
            axs.text(umap1_range[1] - umap1_offset, umap2_range[1] - umap2_offset, 'B', fontsize=30)
            axs.text(umap1_range[0] + umap1_offset, umap2_range[0] + umap2_offset, 'C', fontsize=30)
            axs.text(umap1_range[1] - umap1_offset, umap2_range[0] + umap2_offset, 'D', fontsize=30)
            plt.axhline(umap2_half)
            plt.axvline(umap1_half)

        plt.title(f"{title}", size=20)
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        if color_by is not None:
            axs.legend(title=color_by, bbox_to_anchor=(1.05, 1), loc='upper left')

        return plt

    def extract_tfidf_features(self, drop_duplicates=True):
        input_df = self.data.copy()
        if drop_duplicates:
            input_df['hp_features'] = input_df['hp_features'].apply(lambda x: list(set(x)))
        tfidf_feat = input_df.groupby('cluster_id')['hp_features'].sum().reset_index()
        tfidf_feat = tfidf_feat.explode('hp_features')
        tfidf_feat = pd.DataFrame(
            tfidf_feat.groupby(['cluster_id', 'hp_features'])['hp_features']
                .count()
        ).rename(columns={'hp_features': 'n'}).reset_index()

        tfidf_feat = tfidf_feat.sort_values(by=['cluster_id'], ascending=False).reset_index(drop=True)
        tfidf_feat = tfidf_feat.merge(
            tfidf_feat
                .groupby("hp_features")["cluster_id"]
                .nunique()
                .reset_index()
                .rename(columns={"cluster_id": "tf_i"})
            , on="hp_features", how="left"
        )
        n_clusters = tfidf_feat['cluster_id'].nunique()
        tfidf_feat['tfidf_n'] = tfidf_feat \
            .apply(
            lambda x: x["n"] * np.log(n_clusters / x['tf_i']), axis=1)

        tfidf_feat = tfidf_feat.merge(
            tfidf_feat
                .groupby("cluster_id")['n']
                .max()
                .reset_index()
                .rename(
                columns={"n": "max_tf"}
            )
            , how='left')

        tfidf_feat['norm_tfidf'] = tfidf_feat \
            .apply(
            lambda x: 0.4 + (1 - 0.4) * (x['n'] / x['max_tf']) * np.log(n_clusters / x['tf_i'])
            , axis=1)
        tfidf_feat = tfidf_feat.sort_values("cluster_id").reset_index(drop=True)
        tfidf_feat['index'] = tfidf_feat.index
        self.tfidf_features = tfidf_feat

    def to_quad(self, umap1, umap2, umap1_range, umap2_range):
        umap2_half = umap2_range[0] + (umap2_range[1] - umap2_range[0]) / 2
        umap1_half = umap1_range[0] + (umap1_range[1] - umap1_range[0]) / 2
        if umap2 > umap2_half and umap1 < umap1_half:
            return 'A'
        if umap2 > umap2_half and umap1_half >= 0:
            return 'B'
        if umap2 <= umap2_half and umap1 >= umap1_half:
            return 'D'
        if umap2 <= umap2_half and umap1 <= umap1_half:
            return 'C'

    def label_quads(self):
        umap1_range = (self.data['umap1'].min(), self.data['umap1'].max())
        umap2_range = (self.data['umap2'].min(), self.data['umap2'].max())
        self.data['quad'] = self.data.apply(
            lambda x: self.to_quad(
                x['umap1'],
                x['umap2'],
                umap1_range,
                umap2_range), axis=1)

    def cluster_stats(self):

        selector = SelectKBest(chi2, k='all')
        cluster_labels = LabelBinarizer().fit_transform(self.data['cluster_id'])
        selector.fit(self.feature_array, cluster_labels)

        # Add here Merge features first by hrss before feature selection

        p_values = selector.pvalues_
        scores = -np.log10(selector.pvalues_)

        #
        # top_k_feature_indices = selector.get_support()
        #
        #
        # indx_features = np.where(top_k_feature_indices == True)

        df_pvals = pd.DataFrame(p_values, columns=['pval']).reset_index()
        df_pvals['-log(p_val)'] = scores
        df_pvals['pval'].fillna(p_values[~np.isnan(p_values)].max(), inplace=True)
        df_pvals.rename(columns={"index": "hp_features"}, inplace=True)

        max_score = scores[(~np.isinf(scores)) & (~np.isnan(scores))].max()
        min_score = scores[(~np.isinf(scores)) & (~np.isnan(scores))].min()

        df_pvals['-log(p_val)'].replace(np.inf, max_score, inplace=True)
        df_pvals['-log(p_val)'].fillna(min_score, inplace=True)

        sig_features = self.data.groupby('cluster_id')['hp_features'].sum().reset_index()
        sig_features = sig_features.explode('hp_features')
        sig_features = pd.DataFrame(sig_features \
                                    .groupby(['cluster_id', 'hp_features'])['hp_features'] \
                                    .count()
                                    ).rename(columns={'hp_features': 'n'}).reset_index()
        self.sig_features = sig_features \
            .merge(df_pvals, on='hp_features', how="left") \
            .sort_values("-log(p_val)", ascending=False)
        all_terms = self.data.explode(self.terms_column)[self.terms_column].unique()
        self.sig_features[self.terms_column] = self.sig_features['hp_features'].apply(lambda x: self.feature_to_hps[x])
        self.sig_features[self.terms_column] = self.sig_features[self.terms_column].apply(
            lambda x: set(x) & set(all_terms))
        self.sig_features['hpo_names'] = self.sig_features[self.terms_column].apply(
            lambda x: [self.hpo2name(t) for t in x])

    def plot_tfidf_features(self, label_top_features=2, plot_noise=True):
        plt.figure(figsize=(15, 5))
        fig, axes = plt.subplots(1, 1, figsize=(15, 5))
        original_df = self.data.copy()

        cut_at = 0.01
        cluster_info = []
        self.extract_tfidf_features()
        clusters = sorted(self.tfidf_features['cluster_id'].unique())

        colors = [plt.cm.RdBu_r(each)
                  for each in np.linspace(0, 1, len(clusters))]
        x_labels_pos = []
        if not plot_noise:
            clusters = [x for x in clusters if x >= 0]
        for k, col in zip(clusters, colors):
            if k < 0:
                # Black used for noise.
                col = [0, 0, 0, 1]

            df_sub = self.tfidf_features[self.tfidf_features['cluster_id'] == k].copy()
            df_sub = df_sub.reset_index(drop=True)
            xx = df_sub['index']
            yy = df_sub['norm_tfidf']
            n_size = (df_sub['n'].values / (np.ones((1, yy.shape[0])) * yy.shape[0])) * 100
            x_labels_pos.append((xx.max() - (xx.max() - xx.min()) / 2))
            plt.scatter(xx, yy, marker='o', color=col, s=n_size)
            if label_top_features > 0:
                sigs = df_sub.sort_values("norm_tfidf", ascending=False).head(label_top_features)
            else:
                sigs = self.tfidf_features[
                    (self.tfidf_features['cluster_id'] == k) & (self.tfidf_features['norm_tfidf'] > cut_at)
                    ]
            if not sigs.empty:
                all_terms = original_df.query(f"cluster_id == {k}")[[self.terms_column]] \
                    .explode(self.terms_column)[self.terms_column].to_list()
                for x, y, feat_id in sigs[['index', 'norm_tfidf', 'hp_features']].values:
                    feat_id = int(feat_id)
                    feat_terms = self.feature_to_hps[feat_id]

                    rel_terms = [x for x in all_terms if x in feat_terms]

                    if rel_terms:
                        mc_hpid = Counter(rel_terms).most_common(1)[0][0]
                        text = self.hpo2name(mc_hpid)
                    else:
                        mc_hpid = None
                        text = 'None'
                    if not text:
                        text = "No text"
                    cluster_info.append(dict(cluster=k, hpid=mc_hpid, name=text, tfidf=y))
                    plt.scatter(x, y, marker='o', color=col, edgecolor='k', s=100)
                    plt.text(x, y + 0.01, text.split('(')[0], fontsize=20, color='grey', rotation=90)

        plt.xticks(x_labels_pos, labels=clusters)
        if label_top_features == 0:
            plt.axhline(y=cut_at, color="grey", linestyle='--')

        plt.setp(axes.spines.values(), visible=True)

        axes.patch.set_visible(False)
        plt.xlabel('Cluster Num.', size=20, rotation=0)
        plt.ylabel('TF-IDF (normalized)', size=20)
        plt.xticks(rotation=0, size=20)
        plt.yticks(rotation=0, size=10)
        return plt


