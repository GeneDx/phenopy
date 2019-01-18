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
        # find common breadth-first-search predecessors
        try:
            bfs_predecessors = []
            for term in [term_a, term_b]:
                bfs_predecessors.append({p[0] for p in nx.bfs_predecessors(self.hpo_network, term)})
            common_bfs_predecessors = bfs_predecessors[0].intersection(bfs_predecessors[1])
            # lca node
            return max(common_bfs_predecessors, key=lambda n: self.hpo_network.node[n]['depth'])
        except ValueError:
            raise ValueError(term_a, term_b)

    def calculate_gamma(self, term_a, term_b, term_lca):
        """
        Calculate gamma term for the HRSS algorithm.

        :param term_a: HPO term A.
        :param term_b: HPO term B.
        :param term_lca: Lowest common ancestor term.
        :return: `float` (term pair comparison score)
        """
        # calculate gamma
        if term_a == term_b:
            return 0
        elif term_b in nx.ancestors(self.hpo_network, term_a):
            if nx.shortest_path_length(self.hpo_network, term_b, term_a) == 1:
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

        # find information content for the most informative leaf for each term
        mil_ic = []
        for term in [term_a, term_b]:
            if self.hpo_network.in_edges(term):
                mil_ic.append(max(
                    {self.hpo_network.node[p]['ic'] for p in self.hpo_network.predecessors(term) if 'ic' in self.hpo_network.node[p]}
                ))
            else:
                mil_ic.append(self.hpo_network.node[term]['ic'])

        # calculate beta_ic?
        beta_ic = ((mil_ic[0] - self.hpo_network.node[term_a]['ic']) + (mil_ic[1] - self.hpo_network.node[term_b]['ic'])) / 2.0

        # find lowest common ancestors for the two terms
        lca_node = self.find_lca(term_a, term_b)

        # calculate alpha_ic
        alpha_ic = self.hpo_network.node[lca_node]['ic']

        # calculate gamma
        gamma = self.calculate_gamma(term_a, term_b, lca_node)

        # calculate the pairs score
        pair_score = (1.0 / float(1.0 + gamma)) * (alpha_ic / float(alpha_ic + beta_ic))

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
            [(pair[0], pair[1], self.score_hpo_pair_hrss(pair)) for pair in term_pairs],
            columns=['a', 'b', 'score']
        ).set_index(
            ['a', 'b']
        ).unstack()

        return round((df.max(axis=1).mean() + df.max(axis=0).mean()) / 2.0, 4)

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
