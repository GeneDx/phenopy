import itertools
import networkx as nx
import numpy as np
import pandas as pd

from functools import lru_cache
from phenopy.weights import calculate_age_weights


class Scorer:
    def __init__(self, hpo_network, summarization_method='BMWA', min_score_mask=0.05):
        self.hpo_network = hpo_network
        if summarization_method not in ['BMA', 'BMWA', 'maximum']:
            raise ValueError('Unsupported summarization method, please choose from BMA, BMWA, or maximum.')
        self.summarization_method = summarization_method
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
        for i, term in enumerate([term_a, term_b]):
            parents.append(
                {p[0] for p in nx.bfs_predecessors(self.hpo_network, term)})
            parents[i].add(term)
        common_parents = parents[0].intersection(
            parents[1])
        # lca node
        # find the ancestor with the highest IC
        # break ties by choosing the node with the greatest depth
        return max(common_parents, key=lambda n: (self.hpo_network.nodes[n]['ic'], self.hpo_network.nodes[n]['depth']))

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
                mil = max(leaves, key=lambda n: (self.hpo_network.nodes[n]['ic'], self.hpo_network.nodes[n]['depth']))
                mil_ic.append(self.hpo_network.nodes[mil]['ic'])
            # the node is a leaf
            else:
                mil_ic.append(self.hpo_network.nodes[term]['ic'])

        # calculate beta_ic
        beta_ic = ((mil_ic[0] - self.hpo_network.nodes[term_a]['ic'])
                   + (mil_ic[1] - self.hpo_network.nodes[term_b]['ic'])) / 2.0
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

    @lru_cache(maxsize=72000000)
    def score_hpo_pair_hrss(self, term_a, term_b):
        """
        Scores the comparison of a pair of terms, using Hybrid Relative Specificity Similarity (HRSS) algorithm.

        :param term_a: HPO term A
        :param term_b: HPO term B
        :return: `float` (term pair comparison score)
        """

        # calculate beta_ic
        beta_ic = self.calculate_beta(term_a, term_b)

        # find lowest common ancestors for the two terms
        lca_node = self.find_lca(term_a, term_b)

        # calculate alpha_ic
        alpha_ic = self.hpo_network.nodes[lca_node]['ic']

        if (alpha_ic == 0.0) and (beta_ic == 0.0):
            return 0.0

        gamma = self.calculate_gamma(term_a, term_b, lca_node)
        I = (alpha_ic / (alpha_ic + beta_ic))
        D = (1.0 / (1.0 + gamma))
        return I * D


    def score(self, record_a, record_b):
        """
        Scores the comparison of terms in list A to terms in list B.

        :param record_a: record A.
        :param record_b: record B.
        :return: `float` (comparison score)
        """
        if self.summarization_method not in ['BMA', 'BMWA', 'maximum']:
            raise ValueError('Unsupported summarization method, please choose from BMA, BMWA, or maximum.')

        # if either set is empty return 0.0
        terms_a = record_a['terms']
        terms_b = record_b['terms']
        if not terms_a or not terms_b:
            return 0.0

        # calculate weights for record_a and record_b
        weights_a = record_a['weights'].copy() if record_a['weights'] is not None else []
        weights_b = record_b['weights'].copy() if record_b['weights'] is not None else []

        # set weights
        # if we have age of record_a use it to set age weights for record_b
        if 'age' in record_a:
            weights_b['age'] = calculate_age_weights(record_b['terms'], record_a['age'], self.hpo_network)

        # if we have age of record_b use it to set age weights for record_a
        if 'age' in record_b:
            weights_a['age'] = calculate_age_weights(record_a['terms'], record_b['age'], self.hpo_network)

        term_pairs = itertools.product(terms_a, terms_b)
        df = pd.DataFrame(
            [(pair[0], pair[1], self.score_hpo_pair_hrss(pair[0], pair[1]))
             for pair in term_pairs],
            columns=['a', 'b', 'score']
        ).set_index(
            ['a', 'b']
        ).unstack()

        if self.summarization_method == 'maximum':
            return self.maximum(df)
        elif self.summarization_method == 'BMWA' and any([weights_a, weights_b]):
            return self.best_match_weighted_average(df, weights_a=weights_a, weights_b=weights_b)
        else:
            return self.best_match_average(df)

    def score_records(self, a_records, b_records, record_pairs, thread_index=0, threads=1):
        """
        Score list pair of records.
        :param a_records: Input records dictionary.
        :param b_records: Score against records. If not provided both pairs members are pulled from "a_records".
        :param record_pairs: List of record pairs to score.
        :param thread_index: Thread index for multiprocessing.
        :param threads: Total number of threads for multiprocessing.
        """
        results = []
        # iterate over record pairs starting, stopping, stepping taking multiprocessing threads in consideration
        for record_a, record_b in itertools.islice(record_pairs, thread_index, None, threads):

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
        return np.average(np.append(max1, max0))

    @staticmethod
    def maximum(df):
        """Returns the maximum similarity value between to term lists"""
        return df.values.max()

    def best_match_weighted_average(self, df, weights_a, weights_b):
        """Returns Best-Match Weighted Average of a termlist to termlist similarity matrix."""
        max_a = df.max(axis=1).values
        max_b = df.max(axis=0).values
        scores = np.append(max_a, max_b)

        weights_matrix = {}
        for w in weights_a:
            # init weight list if necessary
            if w not in weights_matrix:
                weights_matrix[w] = []

            # extend weight with the values of a
            weights_matrix[w].extend(weights_a[w])

            # for columns not in b, fill in with 1s for each b row
            if w not in weights_b:
                weights_matrix[w].extend([1 for _ in range(max_b.shape[0])])

        for w in weights_b:
            # for columns not in a, fill in with 1s for each a row
            if w not in weights_matrix:
                weights_matrix[w] = [1 for _ in range(max_a.shape[0])]

            # extend weight with the values of b
            weights_matrix[w].extend(weights_b[w])

        weights_df = pd.DataFrame.from_dict(weights_matrix)
        weights = weights_df.min(axis=1)

        # mask good matches from weighting
        # mask threshold based on >75% of pairwise scores of all hpo terms
        # TODO: expose min_score cutoff value to be set in config
        if self.min_score_mask is not None:
            masked_weights = np.where(scores > self.min_score_mask, 1.0, weights)
            weights = masked_weights

        # if weights add up to zero, calculate unweighted average
        if np.sum(weights) == 0.0:
            weights = np.ones(len(weights))

        return np.average(scores, weights=weights)
