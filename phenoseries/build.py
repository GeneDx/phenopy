import argparse
import os
from pickle import dump
import pandas as pd
import requests

from phenopy import generate_annotated_hpo_network
from phenopy.config import config, logger
from phenopy.util import standardize_phenotypes, remove_parents

from txt2hpo.extract import Extractor

OMIM_API_URL = "https://api.omim.org/api/"
OMIM_DOWNLOADS_URL = "https://data.omim.org/downloads/"


def request_mimid_info(mimid):
    """
    request mimid description from OMIM
    """
    access = "entry?"

    payload = {
        "mimNumber": mimid,
        "include": "text",
        "format": "json",
        "apiKey": config.get("omim", "omim_api_key"),
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


def build_phenoseries_experiment(
    phenotypic_series_filepath=None, outdir=None, min_hpos=4, min_entities=2
):
    """Create a DataFrame that has PS and list of MIM ids."""

    # load HPO network
    # data directory
    phenopy_data_directory = os.path.join(os.getenv("HOME"), ".phenopy/data")

    # files used in building the annotated HPO network
    obo_file = os.path.join(phenopy_data_directory, "hp.obo")
    disease_to_phenotype_file = os.path.join(phenopy_data_directory, "phenotype.hpoa")

    hpo_network, alt2prim, disease_records = generate_annotated_hpo_network(
        obo_file, disease_to_phenotype_file, ages_distribution_file=None
    )

    # read the phenotype.hpoa as a DataFrame
    disease_to_phenotype_file = os.path.join(phenopy_data_directory, "phenotype.hpoa")
    omimdf = pd.read_csv(disease_to_phenotype_file, sep="\t", skiprows=4)

    # process OMIM file
    # only omim annotations
    omimdf = omimdf.loc[omimdf["DatabaseID"].str.contains("OMIM")]
    # strip so only integer remains
    omimdf["mimid"] = omimdf["DatabaseID"].str.strip("OMIM:")

    # dictionary of mim id to list of hpo ids
    mim_hpos_hpoa = {}
    for mim_id, ser in omimdf.groupby("mimid"):
        mim_hpos_hpoa[mim_id] = list(set(ser["HPO_ID"].tolist()))

    # read the phenotypic series file as a DataFrame
    psdf = pd.read_csv(
        phenotypic_series_filepath,
        sep="\t",
        comment="#",
        names=["PS", "MIM", "Phenotype"],
    )
    # null phenotypes are actually null MIM id fields, so just drop these
    psdf = psdf.dropna()
    psdf.reset_index(inplace=True, drop=True)

    # create a dictionary for phenotypic series to list of omim ids mapping
    ps2mimids = {}
    for ps, mim_ids in psdf.groupby(["PS"])["MIM"]:
        # more than two mims in a ps
        if len(mim_ids) >= 2:
            ps2mimids[ps] = list(set(mim_ids.tolist()))

    # invert the ps2mimid dictionary for easy lookup of which ps a mim belongs to
    mim2psids = {}
    for mim_id, ps in psdf.groupby(["MIM"])["PS"]:
        mim2psids[mim_id] = ps.tolist()

    fields_to_use = [
        "text",
        "description",
        "otherFeatures",
        "biochemicalFeatures",
        "diagnosis",
        "clinicalFeatures",
    ]

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

    # do we need this?
    mim_hpos = {str(mim_id): hpos for mim_id, hpos in mim_hpos.items()}

    # clean up HPO ids in lists
    for mim_id, hpo_ids in mim_hpos.items():
        mim_hpos[mim_id] = standardize_phenotypes(hpo_ids, hpo_network, alt2prim)

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
    mimdf['age'] = '.'
    mimdf["hpo_string"] = mimdf["hpo_terms"].apply("|".join)

    # dump the mimdf DataFrame to csv
    mimdf[["omim_id", "age", "hpo_string"]].to_csv(
        os.path.join(outdir, "phenoseries.mim_hpo_ids.txt"),
        sep="\t",
        index=False,
        header=False,
    )

    # dump the dictionary that contains list of mims for each ps
    with open(os.path.join(outdir, "phenoseries.ps_mim_ids.pkl"), "wb") as handle:
        dump(experiment_ps2mimids, handle)


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
        help="The minimum number of hpo ids per entity (mim id, for example) to be considered for the experiment",
    )
    parser.add_argument(
        "--min-entities",
        "-m",
        default=2,
        help="The minimum number of entities (mim id, for example) per series to be considered for the experiment",
    )
    args = parser.parse_args()

    outdir = args.outdir
    phenotypic_series_filepath = args.phenotypic_series_filepath
    min_hpos = args.min_hpos
    min_entities = args.min_entities

    build_phenoseries_experiment(
        phenotypic_series_filepath=phenotypic_series_filepath,
        outdir=outdir,
        min_hpos=min_hpos,
        min_entities=min_entities,
    )
