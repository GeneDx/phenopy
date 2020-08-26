import argparse
import os
import pandas as pd

from itertools import chain
from pickle import load


def make_rank_dataframe(pairwise_sim_matrix, mimdf, ps2mimids):
    relevant_ranks_results = []
    for psid, ps_mim_ids in ps2mimids.items():
        # Grab the index of the "relevant" mim ids
        # Helps identify index in pairwise distance matrix
        ps_mim_idxs = mimdf[mimdf["omim_id"].isin(ps_mim_ids)].index.tolist()
        for query_mim_idx in ps_mim_idxs:
            ranks = return_relevant_ranks(
                pairwise_sim_matrix, query_mim_idx, ps_mim_idxs
            )
            query_mim = mimdf.iloc[query_mim_idx]["omim_id"]
            relevant_ranks_results.append([psid, query_mim, ranks])

    rankdf = pd.DataFrame(
        relevant_ranks_results, columns=["psid", "query_mim_id", "relevant_ranks"]
    )
    # for calculation of recall -- TP / total relevant
    rankdf["total_relevant"] = rankdf.apply(
        lambda row: len(row["relevant_ranks"]), axis=1
    )

    return rankdf


def return_relevant_ranks(pairwise_sim, query_idx, other_mim_indices):
    """Given a pairwise similarity matrix, compute the rank of the similarity between
    a query mim and another mim disease from the same PS.
    """
    other_idxs = other_mim_indices.copy()
    other_idxs.remove(query_idx)
    array = pairwise_sim[query_idx].copy()
    order = array.argsort()
    ranks = order.argsort()
    ranks = max(ranks) - ranks
    return sorted(ranks[other_idxs])


def phenoseries_ranks(
    phenoseries_mim_hpo_file, pairwise_scores_file, ps2mim_file, outdir=None
):

    if outdir is None:
        outdir = os.getcwd()

    # open the dictionary that contains list of mims for each ps
    with open(ps2mim_file, "rb") as handle:
        ps2mimids = load(handle)

    mimdf = pd.read_csv(phenoseries_mim_hpo_file, sep="\t", names=['omim_id', 'age', 'hpo_ids'])

    pairwise_scores = pd.read_csv(
        pairwise_scores_file, header=None, names=["mimid1", "mimid2", "score"], sep="\t"
    )
    pairwise_scores = pairwise_scores.set_index(["mimid1", "mimid2"]).unstack()
    pairwise_scores = (
        pairwise_scores["score"]
        .reset_index(drop=True)
        .fillna(pairwise_scores.T.droplevel(0).reset_index(drop=True))
        .set_index(pairwise_scores.index, drop=True)
    )
    # reindex with the mimdf index
    pairwise_scores = pairwise_scores.reindex(mimdf["omim_id"].astype(int).tolist())
    pairwise_scores = pairwise_scores[mimdf["omim_id"].astype(int).tolist()]

    rankdf = make_rank_dataframe(pairwise_scores.values, mimdf, ps2mimids)
    rankdf.to_csv(
        os.path.join(outdir, "phenoseries.ranks.dataframe.txt"), sep="\t", index=False
    )
    ranks = list(chain.from_iterable(rankdf["relevant_ranks"].tolist()))
    return rankdf, ranks


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--phenoseries-mim-hpo-file",
        "-p",
        default="phenoseries.mim_hpo_ids.txt",
        help="Path to the phenoseries to omim id file output by build.py",
    )
    parser.add_argument(
        "--pairwise-scores-file",
        "-s",
        default="phenoseries.hrss.bma.pairwise.txt",
        help="Path to the pairwise scores array output by score.py (or phenopy command line)",
    )
    parser.add_argument(
        "--ps2mim_file",
        "-m",
        default="phenoseries.ps_mim_ids.pkl",
        help="Path to the ps to omim id dictionary.",
    )
    parser.add_argument(
        "--outdir", "-o", default=os.getcwd(), help="Path where to store the results."
    )
    args = parser.parse_args()

    phenoseries_mim_hpo_file = args.phenoseries_mim_hpo_file
    pairwise_scores_file = args.pairwise_scores_file
    ps2mim_file = args.ps2mim_file
    outdir = args.outdir

    phenoseries_ranks(
        phenoseries_mim_hpo_file, pairwise_scores_file, ps2mim_file, outdir
    )
