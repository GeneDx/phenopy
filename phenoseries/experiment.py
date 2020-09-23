import argparse
import os
import numpy as np
import pandas as pd
import requests
import sys

from ast import literal_eval
from phenopy import generate_annotated_hpo_network
from phenopy.config import config, logger
from phenopy.score import Scorer
from phenopy.util import remove_parents, half_product

from txt2hpo.extract import Extractor

OMIM_API_URL = "https://api.omim.org/api/"
OMIM_DOWNLOADS_URL = "https://data.omim.org/downloads/"


def request_mimid_info(mimid):
    """
    request mimid description from OMIM
    """
    access = "entry?"
    api_key = os.getenv("OMIM_API_KEY")
    if api_key is None:
        api_key = config.get("omim", "omim_api_key")
    payload = {
        "mimNumber": mimid,
        "include": "text",
        "format": "json",
        "apiKey": api_key,
    }

    r = requests.get(OMIM_API_URL + access, params=payload)
    if r.status_code == 200:
        return r
    else:
        logger.critical("Please set the omim_api_key in your phenopy.ini config file")


def convert_and_filter_hpoids(terms, hpo, alt2prim):
    """Given a list of HPO ids, first try to convert synonyms to primary ids,
    then filter if terms are not in the ontology"""
    terms = [alt2prim[term] if term in alt2prim else term for term in terms]
    terms = list(filter(lambda term: term in hpo.nodes, terms))
    terms = remove_parents(terms, hpo)
    return terms


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
    other_idxs = [idx-1 for idx in other_idxs]
    mim_sims = pairwise_sim[query_idx].copy()
    mim_sims_noself = np.delete(mim_sims, [query_idx])
    order = mim_sims_noself.argsort()
    ranks = order.argsort()
    ranks = max(ranks) - ranks
    # convert the ranks to 1-based
    ranks = np.array([r+1 for r in ranks])
    return sorted(ranks[other_idxs])


def run_phenoseries_experiment(outdir=None, phenotypic_series_filepath=None,
    min_hpos=2, min_entities=4, phenoseries_fraction=1.0,
    scoring_method="HRSS", threads=1, omim_phenotypes_file=None, pairwise_mim_scores_file=None):
    
    if outdir is None:
        outdir = os.getcwd
    
    # load HPO network
    # data directory
    phenopy_data_directory = os.path.join(os.getenv("HOME"), ".phenopy/data")

    # files used in building the annotated HPO network
    obo_file = os.path.join(phenopy_data_directory, "hp.obo")
    disease_to_phenotype_file = os.path.join(phenopy_data_directory, "phenotype.hpoa")

    hpo_network, alt2prim, disease_records = generate_annotated_hpo_network(
        obo_file, disease_to_phenotype_file, ages_distribution_file=None
    )

    # read the phenotypic series file as a DataFrame
    psdf = pd.read_csv(
        phenotypic_series_filepath,
        sep="\t",
        comment="#",
        names=["PS", "MIM", "Phenotype"],
    )
    # null phenotypes are actually null MIM id fields, so just drop these
    psdf = psdf.dropna().sample(frac=phenoseries_fraction, random_state=42)
    psdf.reset_index(inplace=True, drop=True)

    # create a dictionary for phenotypic series to list of omim ids mapping
    ps2mimids = {}
    for ps, mim_ids in psdf.groupby(["PS"])["MIM"]:
        # more than two mims in a ps
        if len(mim_ids) >= 2:
            ps2mimids[ps] = list(set([int(mid) for mid in mim_ids.tolist()]))

    # invert the ps2mimid dictionary for easy lookup of which ps a mim belongs to
    mim2psids = {}
    for mim_id, ps in psdf.groupby(["MIM"])["PS"]:
        mim2psids[int(mim_id)] = ps.tolist()

    fields_to_use = [
        "text",
        "description",
        "otherFeatures",
        "biochemicalFeatures",
        "diagnosis",
        "clinicalFeatures",
    ]

    if omim_phenotypes_file == "":
        logger.info("Scraping OMIM Diseases text")
        mim_texts = {}
        for mim_id in mim2psids:
            mim_response = request_mimid_info(mim_id)
            try:
                mim_info = mim_response.json()
            except AttributeError:
                break
            mim_text = mim_info["omim"]["entryList"][0]["entry"]["textSectionList"]

            all_mim_text = ""
            for text_section in mim_text:
                section_name = text_section["textSection"]["textSectionName"]
                if section_name in fields_to_use:
                    # unique_section_names.add(section_name)
                    all_mim_text += " " + text_section["textSection"]["textSectionContent"]

            mim_texts[mim_id] = all_mim_text
        # instantiate txt2hpo's Exctractor class to perform named entity recognition
        extractor = Extractor(remove_negated=True, max_neighbors=3, correct_spelling=False)

        # loop over the MIM ids and extract hpo ids from each MIM's text fields
        mim_hpos = {}
        for mim_id in mim2psids:
            mim_hpos[mim_id] = extractor.hpo(mim_texts[mim_id]).hpids

        mimdf = pd.DataFrame()
        mimdf["omim_id"] = list(mim2psids.keys())
        mimdf["hpo_terms"] = mimdf["omim_id"].apply(lambda mim_id: mim_hpos[mim_id])
        mimdf.to_csv(os.path.join(outdir, "omim_phenotypes.txt"), index=False, sep='\t')

    else:
        logger.info("You passed an OMIM disease to phenotype file")
        try:
            mimdf = pd.read_csv(omim_phenotypes_file, sep="\t")
            mimdf["omim_id"] = mimdf["omim_id"].astype(int)
            mimdf["hpo_terms"] = mimdf["hpo_terms"].apply(literal_eval)
            mim_hpos = dict(zip(mimdf["omim_id"], mimdf["hpo_terms"]))
        except FileNotFoundError:
            sys.exit("Please provide a valid file path")

    # do we need this?
    # mim_hpos = {mim_id: hpos for mim_id, hpos in mim_hpos.items()}

    # clean up HPO ids in lists
    for mim_id, hpo_ids in mim_hpos.items():
        mim_hpos[mim_id] = convert_and_filter_hpoids(hpo_ids, hpo_network, alt2prim)

    # remove entities (mims) that have less than min_hpos
    mims_to_remove = []
    for mim_id, hpo_ids in mim_hpos.copy().items():
        if len(hpo_ids) <= min_hpos:
            mims_to_remove.append(mim_id)

    # Now remove the entities (mim ids) with less than min_hpos
    experiment_ps2mimids = {}
    # remove these mims from ps
    for ps, mimids in ps2mimids.copy().items():
        experiment_ps2mimids[ps] = []
        for ps_mim_id in mimids:
            if ps_mim_id not in mims_to_remove:
                experiment_ps2mimids[ps].append(ps_mim_id)

    # After removing entities, make sure the series has min number of entities
    # get lists of mims and their PS
    remove_these_ps = []
    for ps, mimids in experiment_ps2mimids.items():
        if len(mimids) < min_entities:
            remove_these_ps.append(ps)

    for psid in remove_these_ps:
        del experiment_ps2mimids[psid]

    # Create a unique list of entity ids, for scoring later
    experiment_omims = set()
    for psid, mim_ids in experiment_ps2mimids.items():
        for mim in mim_ids:
            experiment_omims.add(mim)
    experiment_omims = list(experiment_omims)

    # make a DataFrame for entity ids
    mimdf = pd.DataFrame()
    mimdf["omim_id"] = experiment_omims
    mimdf["hpo_terms"] = mimdf["omim_id"].apply(lambda mim_id: mim_hpos[mim_id])

    if pairwise_mim_scores_file == "":
        scorer = Scorer(hpo_network, scoring_method=scoring_method)
        records = [
            {
                "record_id": mim_id,
                "terms": convert_and_filter_hpoids(hpo_terms, hpo_network, alt2prim),
                "weights": {},
            }
            for mim_id, hpo_terms in dict(zip(mimdf["omim_id"], mimdf["hpo_terms"])).items()
        ]

        results = scorer.score_records(
            records, records, half_product(len(records), len(records)), threads=threads
        )

        pairwise_scores = pd.DataFrame(
            results, columns=["mimid1", "mimid2", "phenopy-score"]
        )
        # convert to square form
        pairwise_scores = pairwise_scores.set_index(["mimid1", "mimid2"]).unstack()
        # This pandas method chain fills in the missing scores of the square matrix with the values from the transpose of df.
        pairwise_scores = (
            pairwise_scores["phenopy-score"]
            .reset_index(drop=True)
            .fillna(pairwise_scores.T.droplevel(0).reset_index(drop=True))
            .set_index(pairwise_scores.index, drop=True)
        )
        # reindex with the mimdf index
        pairwise_scores = pairwise_scores.reindex(mimdf["omim_id"].tolist())
        pairwise_scores = pairwise_scores[mimdf["omim_id"].tolist()]
        pd.DataFrame(pairwise_scores).to_csv(os.path.join(outdir, 'phenoseries.psim_matrix.txt'),
                                          sep='\t')
    else:
        pairwise_scores = pd.read_csv(pairwise_mim_scores_file, sep='\t')

    ranksdf = make_rank_dataframe(pairwise_scores.astype(float).values, mimdf, experiment_ps2mimids)
    ranksdf.to_csv(os.path.join(outdir, "phenoseries.rankdf.txt"), sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--outdir", "-o", default=os.getcwd(), help="Path where to store the results."
    )
    parser.add_argument(
        "--phenotypic-series-filepath",
        "-p",
        help="path to the omim text file defining phenotypic series to omim id",
    )
    parser.add_argument(
        "--min-hpos",
        "-n",
        default=4,
        type=int,
        help="The minimum number of hpo ids per entity (mim id, for example) to be considered for the experiment",
    )
    parser.add_argument(
        "--min-entities",
        "-m",
        default=2,
        type=int,
        help="The minimum number of entities (mim id, for example) per series to be considered for the experiment",
    )
    parser.add_argument(
        "--phenoseries-fraction",
        "-f",
        default=1.0,
        help="The fraction of phenoseries to use",
        type=float,
    )
    parser.add_argument(
        "--scoring-method",
        "-s",
        default="HRSS",
        help="The scoring method to use",
        type=str,
    )
    parser.add_argument(
        "--threads", "-t", default=4, help="The number of threads to use", type=int,
    )
    parser.add_argument(
        "--omim-phenotypes-file",
        "-a",
        default="",
        help="The full path to a pre-generated omim id to list of phenotypes file",
        type=str,
    )
    parser.add_argument(
        "--pairwise-mim-scores-file",
        "-b",
        default="",
        help="The full path to a pre-generated file with all the pairwise scores for each omim id in the experiment.",
        type=str,
    )

    args = parser.parse_args()

    outdir = args.outdir
    phenotypic_series_filepath = args.phenotypic_series_filepath
    min_hpos = args.min_hpos
    min_entities = args.min_entities
    phenoseries_fraction = args.phenoseries_fraction
    scoring_method = args.scoring_method
    threads = args.threads
    omim_phenotypes_file = args.omim_phenotypes_file
    pairwise_mim_scores_file = args.pairwise_mim_scores_file

    run_phenoseries_experiment(
        outdir = outdir,
        phenotypic_series_filepath = phenotypic_series_filepath,
        min_hpos = min_hpos,
        min_entities = min_entities,
        phenoseries_fraction = phenoseries_fraction,
        scoring_method = scoring_method,
        threads = threads,
        omim_phenotypes_file = omim_phenotypes_file,
        pairwise_mim_scores_file = pairwise_mim_scores_file,
        )
