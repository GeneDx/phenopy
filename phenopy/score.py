import itertools
import gensim
import networkx as nx
import numpy as np
import pandas as pd

from functools import lru_cache
from multiprocessing import Pool
from phenopy.weights import calculate_age_weights
from phenopy.config import config
from typing import (
    Tuple,
    Dict,
    List,
    Set,
)


class Scorer:
    def __init__(
        self,
        hpo_network: nx.MultiDiGraph,
        summarization_method: str = "BMWA",
        min_score_mask: float = 0.05,
        scoring_method: str = "HRSS",
    ) -> None:

        # Establish hpo_network
        self.hpo_network = hpo_network

        # Establish summarization method
        if summarization_method not in ["BMA", "BMWA", "maximum"]:
            raise ValueError(
                "Unsupported summarization method, please choose from "
                "BMA, BMWA, or maximum."
            )
        self.summarization_method = summarization_method

        # Assign min_score_mask
        self.min_score_mask = min_score_mask

        # Assign scoring method
        if scoring_method not in ["HRSS", "Resnik", "Jaccard", "word2vec"]:
            raise ValueError(
                "Unsupported semantic similarity scoring method, please "
                "choose from HRSS, Resnik, Jaccard, or word2vec."
            )
        self.scoring_method = scoring_method

        # Load the word vectors if using word2vec
        if scoring_method == "word2vec":
            try:
                self.word_vectors = gensim.models.KeyedVectors.load(
                    config.get("models", "phenopy.wv.model")
                )
            except FileNotFoundError:
                raise ValueError(
                    "Please make sure that a word2vec model is in "
                    "your project data directory."
                )

    def find_lca(self, term_a: str, term_b: str) -> str:
        """
        Determine the lowest common ancestor for two HPO terms
        """

        # if either term is HP:0000001 return it
        if any(term == "HP:0000001" for term in [term_a, term_b]):
            return "HP:0000001"

        # if one of the terms is a child of the other return the parent
        if self.hpo_network.has_edge(term_a, term_b):
            return term_b
        if self.hpo_network.has_edge(term_b, term_a):
            return term_a

        # find common breadth-first-search predecessors
        parents = []
        for i, term in enumerate([term_a, term_b]):
            parents.append({p[0] for p in nx.bfs_predecessors(self.hpo_network, term)})
            parents[i].add(term)

        # Find the intersection between the two sets of parents
        common_parents = parents[0].intersection(parents[1])

        # lca node - find the ancestor with the highest IC
        # break ties by choosing the node with the greatest depth
        return max(
            common_parents,
            key=lambda n: (
                self.hpo_network.nodes[n]["ic"],
                self.hpo_network.nodes[n]["depth"],
            ),
        )

    def calculate_beta(self, term_a: str, term_b: str) -> float:
        """
        calculates the beta term in HRSS equation
        """
        # find information content for the most informative leaf for each term
        mil_ic = []

        # For each term, if it has children, find the most informative leaf
        for term in [term_a, term_b]:
            if self.hpo_network.in_edges(term):

                # children terms generator
                children = nx.ancestors(self.hpo_network, term)

                # Establish the leaf nodes
                leaves = {
                    p
                    for p in children
                    if self.hpo_network.out_degree(p) >= 1
                    and self.hpo_network.in_degree(p) == 0
                }

                # append the max IC leaf (choose the one with the max depth)
                mil = max(
                    leaves,
                    key=lambda n: (
                        self.hpo_network.nodes[n]["ic"],
                        self.hpo_network.nodes[n]["depth"],
                    ),
                )
                mil_ic.append(self.hpo_network.nodes[mil]["ic"])

            # the node is a leaf
            else:
                mil_ic.append(self.hpo_network.nodes[term]["ic"])

        # calculate beta_ic
        beta_ic = (
            (mil_ic[0] - self.hpo_network.nodes[term_a]["ic"])
            + (mil_ic[1] - self.hpo_network.nodes[term_b]["ic"])
        ) / 2.0

        return beta_ic

    def calculate_gamma(self, term_a: str, term_b: str, term_lca: str) -> int:
        """
        Calculate gamma term for the HRSS algorithm.
        """
        # calculate gamma
        # "such that the value equals zero if the two terms are the same"
        if term_a == term_b:
            return 0

        # If one of the terms is a child of the other return 1
        term_a_child = self.hpo_network.has_edge(term_a, term_b)
        term_b_child = self.hpo_network.has_edge(term_b, term_a)
        if term_a_child or term_b_child:
            return 1

        # Otherwise calculate the shortest-path length to the LCA
        a_to_lca = nx.shortest_path_length(self.hpo_network, term_a, term_lca)
        b_to_lca = nx.shortest_path_length(self.hpo_network, term_b, term_lca)

        return a_to_lca + b_to_lca

    @lru_cache(maxsize=72000000)
    def score_hpo_pair_hrss(self, term_a: str, term_b: str) -> float:
        """
        Scores the comparison of a pair of terms, using Hybrid Relative Specificity
        Similarity (HRSS) algorithm.
        """

        # calculate beta_ic
        beta_ic = self.calculate_beta(term_a, term_b)

        # find lowest common ancestors for the two terms
        lca_node = self.find_lca(term_a, term_b)

        # calculate alpha_ic
        alpha_ic = self.hpo_network.nodes[lca_node]["ic"]
        if self.scoring_method == "Resnik":
            return alpha_ic

        # Return 0 if alpha_ic and beta_ic are both 0
        if (alpha_ic == 0.0) and (beta_ic == 0.0):
            return 0.0

        # calculate gamma
        gamma = self.calculate_gamma(term_a, term_b, lca_node)

        # Assign the I and D variables in the HRSS equation
        i_variable = alpha_ic / (alpha_ic + beta_ic)
        d_variable = 1.0 / (1.0 + gamma)

        return i_variable * d_variable

    def score(self, record_a: Dict, record_b: Dict) -> Tuple[str, str, float]:
        """
        Scores the comparison of terms listed in record A to terms listed in record B.
        """
        if self.summarization_method not in ["BMA", "BMWA", "maximum"]:
            raise ValueError(
                "Unsupported summarization method, please choose from "
                "BMA, BMWA, or maximum."
            )

        # if either set is empty return 0.0
        terms_a = record_a["terms"]
        terms_b = record_b["terms"]
        if not terms_a or not terms_b:
            return record_a["record_id"], record_b["record_id"], 0.0

        # If specified, calculate the Jaccard similarity
        if self.scoring_method == "Jaccard":
            intersection = len(list(set(terms_a).intersection(terms_b)))
            union = (len(terms_a) + len(terms_b)) - intersection
            comparison_score = float(intersection) / union
            return record_a["record_id"], record_b["record_id"], comparison_score

        # If specified, calculate the word2vec similarity
        elif self.scoring_method == "word2vec":

            # Ensure that all HPO terms are in the vocab
            in_vocab_terms_a = [
                x for x in terms_a if x in self.word_vectors.key_to_index
            ]
            in_vocab_terms_b = [
                x for x in terms_b if x in self.word_vectors.key_to_index
            ]

            # If both records have terms in the vocab (both are non-empty lists)
            if in_vocab_terms_a and in_vocab_terms_b:
                return self.word_vectors.n_similarity(
                    in_vocab_terms_a, in_vocab_terms_b
                )

            # One record or the other has no terms in the word2vec vocab
            else:
                return record_a["record_id"], record_b["record_id"], 0.0

        # calculate weights for record_a and record_b
        if record_a["weights"] is not None:
            weights_a = record_a["weights"].copy()
        else:
            weights_a = []
        if record_b["weights"] is not None:
            weights_b = record_b["weights"].copy()
        else:
            weights_b = []

        # set weights
        # if we have age of record_a use it to set age weights for record_b
        if "age" in record_a:
            weights_b["age"] = calculate_age_weights(
                record_b["terms"], record_a["age"], self.hpo_network
            )

        # if we have age of record_b use it to set age weights for record_a
        if "age" in record_b:
            weights_a["age"] = calculate_age_weights(
                record_a["terms"], record_b["age"], self.hpo_network
            )

        # Creates a dataframe that houses the HRSS for each term pair
        df = self.get_term_pair_dataframe(terms_a, terms_b)

        # Return maximum if specified
        if self.summarization_method == "maximum":
            return record_a["record_id"], record_b["record_id"], self.maximum(df)

        # Retrun BMWA if specified
        elif self.summarization_method == "BMWA" and any([weights_a, weights_b]):
            score_output = self.best_match_weighted_average(
                df, weights_a=weights_a, weights_b=weights_b
            )

            return record_a["record_id"], record_b["record_id"], score_output

        # Otherwise return the best-match-average
        else:
            score_output = self.best_match_average(df)
            return record_a["record_id"], record_b["record_id"], score_output

    def score_term_sets_basic(self, terms_a: str, terms_b: str) -> float:
        """
        Calculate the semantic similarity of two lists of terms.
        This is intended to be used as a library function. It is not used by the CLI.
        """
        # Instantiate the two lists of HPO identifiers
        terms_a = set(terms_a)
        terms_b = set(terms_b)

        # Calculate the Jaccard similarity if specified
        if self.scoring_method == "Jaccard":
            intersection = len(list(set(terms_a).intersection(terms_b)))
            union = (len(terms_a) + len(terms_b)) - intersection
            return float(intersection) / union

        # Calculate the word vector similarity if word2vec is specified
        elif self.scoring_method == "word2vec":

            # Instantiate a list to house all HPO terms that are within the vocab
            in_vocab_terms_a = [
                x for x in terms_a if x in self.word_vectors.key_to_index
            ]
            in_vocab_terms_b = [
                x for x in terms_b if x in self.word_vectors.key_to_index
            ]

            # If both lists exist (both are non-empty lists) return their similarity
            if in_vocab_terms_a and in_vocab_terms_b:
                return self.word_vectors.n_similarity(
                    in_vocab_terms_a, in_vocab_terms_b
                )

            # Otherwise return 0.0
            else:
                return 0.0

        # Creates a dataframe that houses the HRSS for each term pair
        df = self.get_term_pair_dataframe(terms_a, terms_b)

        # If set to maximum, return the maximum, otherwise best-match-average
        if self.summarization_method == "maximum":
            return self.maximum(df)
        else:
            return self.best_match_average(df)

    def score_records(
        self, a_records: Dict, b_records: Dict, record_pairs: List, threads: int = 1
    ) -> List:
        """
        Scores a pair of records based on the specified number of threads
        """
        with Pool(processes=threads) as p:
            results = p.starmap(
                self.score,
                [
                    (
                        a_records[record_a],  # a records
                        b_records[record_b],  # b records
                    )
                    for (record_a, record_b) in record_pairs
                ],
            )

        return results

    @staticmethod
    def best_match_average(df: pd.DataFrame) -> float:
        """
        Returns the Best-Match average of a termlist to termlist similarity matrix.
        """
        # Determine the max values of the rows and columns
        max_column_values = df.max(axis=1).values
        max_row_values = df.max(axis=0).values
        return np.average(np.append(max_column_values, max_row_values))

    @staticmethod
    def maximum(dataframe: pd.DataFrame) -> float:
        """Returns the maximum similarity value between to term lists"""
        return dataframe.values.max()

    def best_match_weighted_average(
        self, df: pd.DataFrame, weights_a: Dict, weights_b: Dict
    ) -> float:
        """
        Returns Best-Match Weighted Average of a termlist to termlist similarity matrix.
        """
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
            # for columns not in a fill in with 1s for each a row
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

    def get_term_pair_dataframe(self, terms_a: Set, terms_b: Set) -> pd.DataFrame:
        """
        Creates a dataframes of pairwise HRSS scores between them
        """
        # Create the list of term pairs
        # e.g., ['a', 'b']['c', 'd'] -> [('a', 'c'), ('a', 'd'), ('b', 'c), ('b','d')]
        term_pairs = itertools.product(terms_a, terms_b)

        # Apply the HRSS score to each pair within the dataframe
        dataframe = (
            pd.DataFrame(
                [
                    (pair[0], pair[1], self.score_hpo_pair_hrss(pair[0], pair[1]))
                    for pair in term_pairs
                ],
                columns=["a", "b", "score"],
            )
            .set_index(["a", "b"])
            .unstack()
        )

        return dataframe
