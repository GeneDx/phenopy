import itertools
import sys

import networkx as nx
import pandas as pd


class Scorer:
    def __init__(self, hpo_network):
        self.hpo_network = hpo_network

        self.scores_cache = {}

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
        return max(common_parents, key=lambda n: self.hpo_network.node[n]['depth'])

    def calculate_beta(self, term_a, term_b):
        """calculates the gamma term in HRSS equation

        :param term_a: Any HPO term.
        :param term_b: Any HPO term.
        :return: `float` gamma
        """
        # find information content for the most informative leaf for each term
        mil_ic = []
        for term in [term_a, term_b]:
            if self.hpo_network.in_edges(term):
                # children terms generator
                children = nx.ancestors(self.hpo_network, term)
                if children:
                    # append the max IC leaf
                    mil_ic.append(max({self.hpo_network.node[p]['ic'] for p in children if self.hpo_network.out_degree(
                        p) >= 1 and self.hpo_network.in_degree(p) == 0}))
                # node is a leaf
                else:
                    mil_ic.append(self.hpo_network.node[term]['ic'])
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

        # calculat beta_ic
        beta_ic = self.calculate_beta(term_a, term_b)

        # find lowest common ancestors for the two terms
        lca_node = self.find_lca(term_a, term_b)

        # calculate alpha_ic
        alpha_ic = self.hpo_network.node[lca_node]['ic']

        # calculate gamma
        gamma = self.calculate_gamma(term_a, term_b, lca_node)

        # calculate the pairs score
        pair_score = (1.0 / float(1.0 + gamma)) * \
            (alpha_ic / float(alpha_ic + beta_ic))

        # cache this pair score
        self.scores_cache[f'{term_a}-{term_b}'] = pair_score

        return pair_score

    def score(self, terms_a, terms_b):
        """
        Scores the comparison of terms in list A to terms in list B.

        :param terms_a: List of HPO terms A.
        :param terms_b: List of HPO terms B.
        :return: `float` (comparison score)
        """
        # filter out hpo terms not in the network and unique them
        terms_a = set(filter(lambda x: x in self.hpo_network.node, terms_a))
        terms_b = set(filter(lambda x: x in self.hpo_network.node, terms_b))

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

        return round(((df.max(axis=1).sum() + df.max(axis=0).sum()) / (len(df.index) + len(df.columns))), 4)

    def score_pairs(self, records, record_pairs, lock, thread=0, number_threads=1):
        """
        Score list pair of records.

        :param records: Records dictionary.
        :param record_pairs: List of record pairs to score.
        :param thread: Thread index for multiprocessing.
        :param number_threads: Total number of threads for multiprocessing.
        :return: `list` of `tuples`
        """
        # iterate over record pairs starting, stopping, stepping taking multiprocessing threads in consideration
        for record_a, record_b in itertools.islice(record_pairs, thread, None, number_threads):
            score = self.score(records[record_a], records[record_b])
            lock.acquire()
            try:
                sys.stdout.write('\t'.join([
                    f'{record_a}-{record_b}',
                    str(score),
                ]))
                sys.stdout.write('\n')
            finally:
                sys.stdout.flush()
                lock.release()
