import itertools
import networkx as nx
import numpy as np
import pandas as pd

from phenopy.weights import age_to_weights


class Scorer:
    def __init__(self, hpo_network, agg_score='BMA', min_score_mask=0.05):
        self.hpo_network = hpo_network
        self.scores_cache = {}
        self.agg_score = agg_score
        self.min_score_mask = min_score_mask

    def find_lca(self, term_a, term_b):
        """
        Determine the lowest common ancestor for a two terms

        :param term_a: HPO term A.
        :param term_b: HPO term B.
        :return: Least Common Ancestor for two terms, "HP:0000001"
        """
        # if either term is HP:0000001 return it
        if any(term == 'HP:0000001' for term in [term_a, term_b]):
            return 'HP:0000001'

        # if one of the terms is a child of the other return the parent
        if self.hpo_network.has_edge(term_a, term_b):
            return term_b
        if self.hpo_network.has_edge(term_b, term_a):
            return term_a

        # find common breadth-first-search predecessors
        parents = []
        for term in [term_a, term_b]:
            parents.append(
                {p[0] for p in nx.bfs_predecessors(self.hpo_network, term)})
        common_parents = parents[0].intersection(
            parents[1])
        # lca node
        # find the ancestor with the highest IC
        # break ties by choosing the node with the greatest depth
        return max(common_parents, key=lambda n: (self.hpo_network.node[n]['ic'], self.hpo_network.node[n]['depth']))

    def calculate_beta(self, term_a, term_b):
        """calculates the beta term in HRSS equation

        :param term_a: Any HPO term.
        :param term_b: Any HPO term.
        :return: `float` beta
        """
        # find information content for the most informative leaf for each term
        mil_ic = []
        for term in [term_a, term_b]:
            if self.hpo_network.in_edges(term):
                # children terms generator
                children = nx.ancestors(self.hpo_network, term)
                # append the max IC leaf (choose the one with the max depth)
                leaves = {p for p in children if self.hpo_network.out_degree(
                    p) >= 1 and self.hpo_network.in_degree(p) == 0}
                mil = max(leaves, key=lambda n: (self.hpo_network.node[n]['ic'], self.hpo_network.node[n]['depth']))
                mil_ic.append(self.hpo_network.node[mil]['ic'])
            # the node is a leaf
            else:
                mil_ic.append(self.hpo_network.node[term]['ic'])

        # calculate beta_ic
        beta_ic = ((mil_ic[0] - self.hpo_network.node[term_a]['ic'])
                   + (mil_ic[1] - self.hpo_network.node[term_b]['ic'])) / 2.0
        return beta_ic

    def calculate_gamma(self, term_a, term_b, term_lca):
        """
        Calculate gamma term for the HRSS algorithm.

        :param term_a: HPO term A.
        :param term_b: HPO term B.
        :param term_lca: Lowest common ancestor term.
        :return: `int` (term pair distance to lca)
        """
        # calculate gamma
        # "such that the value equals zero if the two terms are the same"
        if term_a == term_b:
            return 0

        # if one of the terms is a child of the other
        term_a_child = False
        term_b_child = False

        if self.hpo_network.has_edge(term_a, term_b):
            term_b_child = True
        if self.hpo_network.has_edge(term_b, term_a):
            term_a_child = True

        if term_a_child or term_b_child:
            return 1

        a_to_lca = nx.shortest_path_length(self.hpo_network, term_a, term_lca)
        b_to_lca = nx.shortest_path_length(self.hpo_network, term_b, term_lca)

        return a_to_lca + b_to_lca

    def score_hpo_pair_hrss(self, terms):
        """
        Scores the comparison of a pair of terms, using Hybrid Relative Specificity Similarity (HRSS) algorithm.

        :param terms: HPO terms A and B.
        :return: `float` (term pair comparison score)
        """
        term_a, term_b = terms

        if f'{term_a}-{term_b}' in self.scores_cache:
            return self.scores_cache[f'{term_a}-{term_b}']

        elif f'{term_b}-{term_a}' in self.scores_cache:
            return self.scores_cache[f'{term_b}-{term_a}']

        # calculate beta_ic
        beta_ic = self.calculate_beta(term_a, term_b)

        # find lowest common ancestors for the two terms
        lca_node = self.find_lca(term_a, term_b)

        # calculate alpha_ic
        alpha_ic = self.hpo_network.node[lca_node]['ic']

        # calculate gamma
        gamma = self.calculate_gamma(term_a, term_b, lca_node)

        # calculate the pairs score
        ic = (alpha_ic / float(alpha_ic + beta_ic))
        pair_score = (1.0 / float(1.0 + gamma)) * ic

        # cache this pair score
        self.scores_cache[f'{term_a}-{term_b}'] = pair_score

        return pair_score

    def score(self, record_a, record_b):
        """
        Scores the comparison of terms in list A to terms in list B.

        :param terms_a: List of HPO terms A.
        :param terms_b: List of HPO terms B.
        :return: `float` (comparison score)
        """

        # if either set is empty return 0.0
        terms_a = record_a['terms']
        terms_b = record_b['terms']
        if not terms_a or not terms_b:
            return 0.0

        # calculate weights for record_a and record_b
        weights_a = record_a['weights'] if record_a['weights'] is not None else []
        weights_b = record_b['weights'] if record_b['weights'] is not None else []

        term_pairs = itertools.product(terms_a, terms_b)
        df = pd.DataFrame(
            [(pair[0], pair[1], self.score_hpo_pair_hrss(pair))
             for pair in term_pairs],
            columns=['a', 'b', 'score']
        ).set_index(
            ['a', 'b']
        ).unstack()
        return

        # if self.agg_score == 'BMA':
        #     return self.best_match_average(df)
        # elif self.agg_score == 'maximum':
        #     return self.maximum(df)
        # elif self.agg_score == 'BMWA':
        #     # age weights scoring for scrore product
        #     if len(weights) == 2:
        #         return self.bmwa(df, weights_a=weights_a, weights_b=weights_b)
        #     # disease weights scoring for score
        #     elif len(weights) == 1:
        #         return self.bmwa(df, weights_a=np.ones(df.shape[0]), weights_b=weights_b)
        #     else:
        #         sys.stderr.write('weights cannot be an empty list or have more than two elements.')
        #         sys.exit(1)
        # else:
        #     return 0.0

    def score_records(self, a_records, b_records, record_pairs, thread_index=0, threads=1, use_weights=False):
        """
        Score list pair of records.
        :param a_records: Input records dictionary.
        :param b_records: Score against records. If not provided both pairs members are pulled from "a_records".
        :param record_pairs: List of record pairs to score.
        :param thread_index: Thread index for multiprocessing.
        :param threads: Total number of threads for multiprocessing.
        :param use_weights: Use disease weights.
        """
        results = []
        # iterate over record pairs starting, stopping, stepping taking multiprocessing threads in consideration
        for record_a, record_b in itertools.islice(record_pairs, thread_index, None, threads):
            # set weights if needed
            weights = None
            # if use_weights is True:
            #     weights = [self.get_disease_weights(b_records[record_b]['terms'], b_records[record_b]['record_id'])]

            score = self.score(
                a_records[record_a],
                b_records[record_b],
            )

            results.append((
                a_records[record_a]['record_id'],
                b_records[record_b]['record_id'],
                str(score),
            ))
        return results

    @staticmethod
    def best_match_average(df):
        """Returns the Best-Match average of a termlist to termlist similarity matrix."""
        max1 = df.max(axis=1).values
        max0 = df.max(axis=0).values
        return np.average(np.append(max1, max0)).round(4)

    @staticmethod
    def maximum(df):
        """Returns the maximum similarity value between to term lists"""
        return df.values.max().round(4)

    def bmwa(self, df, weights_a, weights_b):
        """Returns Best-Match Weighted Average of a termlist to termlist similarity matrix."""
        max1 = df.max(axis=1).values
        max0 = df.max(axis=0).values

        scores = np.append(max1, max0)
        weights = np.array(np.append(weights_a, weights_b))

        # mask good matches from weighting
        # mask threshold based on >75% of pairwise scores of all hpo terms
        # TODO: expose min_score cutoff value to be set in config
        if self.min_score_mask is not None:
            masked_weights = np.where(scores > self.min_score_mask, 1.0, weights)
            weights = masked_weights

        # if weights add up to zero, calculate unweighted average
        if np.sum(weights) == 0.0:
            weights = np.ones(len(weights))

        return np.average(scores, weights=weights).round(4)

    def calculate_age_weights(self, terms, age):
        """
        Calculates an age-based weight vector given an iterable of terms.
        :param terms: iterable of hpo terms
        :param age: numeric age of patient
        :return: list of weights in same order as terms
        """
        weights = []
        for node_id in terms:
            if node_id not in self.hpo_network.node:
                weights.append(1.0)
            elif self.hpo_network.node[node_id]['weights']['age_exists']:
                weights.append(age_to_weights(self.hpo_network.node[node_id]['weights']['age_dist'], age))
            else:
                weights.append(1.0)
        return weights

    def get_disease_weights(self, terms, disease_id):
        """Lookup the disease weights for the input terms to the disease_id"""
        return [self.hpo_network.node[node_id]['weights']['disease_frequency'][disease_id] for node_id in terms]
