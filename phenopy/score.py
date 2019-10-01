import itertools
import sys

import networkx as nx
import numpy as np
import pandas as pd

from itertools import product
from phenopy.weights import age_to_weights

class Scorer:
    def __init__(self, hpo_network, agg_score='BMA'):
        self.hpo_network = hpo_network

        self.scores_cache = {}

        self.alt2prim = {}
        self.generate_alternate_ids()
        self.agg_score = agg_score

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
        return max(common_parents, key=lambda n: self.hpo_network.node[n]['ic'])

    def generate_alternate_ids(self):
        """Create a key, value store of alternate terms to canonical terms."""
        for n in self.hpo_network.nodes(data=True):
            n = n[0]
            try:
                for alt in self.hpo_network.node[n]['alt_id']:
                    self.alt2prim[alt] = n
            except KeyError:
                # no alternate HPO ids for this term
                continue

    def convert_alternate_ids(self, termlist):
        """return a list of terms with list of alternate HPO ids converted to canonical ones."""
        return [self.alt2prim[t] if t in self.alt2prim else t for t in termlist]

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
                # append the max IC leaf
                mil_ic.append(max({self.hpo_network.node[p]['ic'] for p in children if self.hpo_network.out_degree(
                    p) >= 1 and self.hpo_network.in_degree(p) == 0}))
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

        if (term_a_child or term_b_child):
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

    def score(self, terms_a, terms_b, weights=[]):
        """
        Scores the comparison of terms in list A to terms in list B.

        :param terms_a: List of HPO terms A.
        :param terms_b: List of HPO terms B.
        :param agg_score: The aggregation method to use for summarizing the similarity matrix between two term sets
            Must be one of {'BMA', 'BMWA'}
        :param weights: List (length=2) of weights, one for each dimension of score matrix.
        :return: `float` (comparison score)
        """
        # convert alternate HPO ids to canonical ones
        terms_a = set(self.convert_alternate_ids(terms_a))
        terms_b = set(self.convert_alternate_ids(terms_b))

        # filter out hpo terms not in the network and unique them
        terms_a = list(filter(lambda x: x in self.hpo_network.node, terms_a))
        terms_b = list(filter(lambda x: x in self.hpo_network.node, terms_b))

        # if either set is empty return 0.0
        if not terms_a or not terms_b:
            return 0.0

        term_pairs = itertools.product(terms_a, terms_b)
        df = pd.DataFrame(
            [(pair[0], pair[1], self.score_hpo_pair_hrss(pair))
             for pair in term_pairs],
            columns=['a', 'b', 'score']
        ).set_index(
            ['a', 'b']
        ).unstack()

        if self.agg_score == 'BMA':
            return self.best_match_average(df)
        elif self.agg_score == 'maximum':
            return self.maximum(df)
        elif self.agg_score == 'BMWA' and len(weights) == 2:
            return self.bmwa(df, weights_a=weights[0], weights_b=weights[1])
        else:
            return 0.0

    def score_pairs(self, records, lock, thread=0, number_threads=1, stdout=True):
        """
        Score list pair of records.

        :param records: Records dictionary.
        :param record_pairs: List of record pairs to score.
        :param weight_method: (age) List of methods by which to adjust score.
        :param thread: Thread index for multiprocessing.
        :param number_threads: Total number of threads for multiprocessing.
        :param stdout:(True,False) write results to standard out
        :return: `list` of `tuples`
        """
        # iterate over record pairs starting, stopping, stepping taking multiprocessing threads in consideration
        results = []

        record_pairs = product([x['sample'] for x in records], repeat=2)

        record_terms = {x['sample']: x['terms'] for x in records}

        for record_a, record_b in itertools.islice(record_pairs, thread, None, number_threads):

            if self.agg_score == 'BMWA':

                record_age = {x['sample']: x['age'] for x in records}

                weights_a = self.calculate_age_weights(record_terms[record_a], record_age[record_b])
                weights_b = self.calculate_age_weights(record_terms[record_b], record_age[record_a])

                score = self.score(record_terms[record_a],
                               record_terms[record_b], weights=[weights_a, weights_b])
            else:

                score = self.score(record_terms[record_a],
                                   record_terms[record_b])

            if stdout:
                try:
                    lock.acquire()
                    sys.stdout.write('\t'.join([record_a, record_b, str(score)]))
                    sys.stdout.write('\n')
                finally:
                    sys.stdout.flush()
                    lock.release()
            else:
                results.append((record_a, record_b, str(score)))

        return results

    def best_match_average(self, df):
        """Returns the Best-Match average of a termlist to termlist similarity matrix."""

        max1 = df.max(axis=1).values
        max0 = df.max(axis=0).values
        return np.average(np.concatenate((max1, max0)))

    def maximum(self, df):
        """Returns the maximum similarity value between to term lists"""
        return df.values.max().round(4)

    def bmwa(self, df, weights_a, weights_b):
        """Returns Best-Match Weighted Average of a termlist to termlist similarity matrix."""

        max1 = df.max(axis=1).values
        max0 = df.max(axis=0).values
        return round(np.average(np.concatenate((max1, max0)), weights=np.concatenate((weights_a, weights_b))), 4)

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
