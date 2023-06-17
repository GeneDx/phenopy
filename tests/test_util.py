import os
import pytest

from phenopy.util import parse, read_records_file, encode_phenotypes, parse_input


def test_parse(test_data):
    string = "age=13;sex=Male"
    assert parse(string, what="sex") == "Male"
    assert parse(string, what="age") == 13.0

    string = "age=13.64;sex=male"
    assert parse(string, what="sex") == "Male"
    assert parse(string, what="age") == 13.6

    string = "age=12.9;sex=female"
    assert parse(string, what="sex") == "Female"
    assert parse(string, what="age") == 12.9

    string = "sex=Female"
    assert parse(string, what="sex") == "Female"

    string = "sex=FEMALE"
    assert parse(string, what="sex") == "Female"

    string = "sex=F"
    assert parse(string, what="sex") == "Female"

    string = "age=1"
    assert parse(string, what="age") == 1.0

    string = "."
    assert not parse(string, what="age")

    string = ". "
    assert not parse(string, what="age")

    string = " . "
    assert not parse(string, what="age")

    string = "13?"
    assert not parse(string, what="age")

    string = "sex=NA"
    assert not parse(string, what="sex")

    string = "sex=Unknown"
    assert not parse(string, what="sex")


def test_encode_phenotypes_file(test_data):
    input_file = os.path.join(test_data["parent_dir"], "data/test.score-short.txt")
    records = parse_input(input_file, test_data["hpo_network"], test_data["alt2prim"])
    encoded_phenotypes = encode_phenotypes(
        [record["terms"] for record in records],
        test_data["phenotype_groups"],
        test_data["hpo_network"],
        test_data["alt2prim"],
    )
    assert sum(encoded_phenotypes[0]) == 4


def test_encode_1d_phenotypes(test_data):
    phenotypes = ["HP:0012759", "HP:0003011", "HP:0011442"]
    encoded_phenotypes = encode_phenotypes(
        phenotypes,
        test_data["phenotype_groups"],
        test_data["hpo_network"],
        test_data["alt2prim"],
        k=1000,
    )
    assert sum(encoded_phenotypes) == 3


def test_encode_2d_phenotypes(test_data):
    phenotypes = [
        ["HP:0012759", "HP:0003011", "HP:0011442"],
        ["HP:0012759", "HP:0003011"],
    ]
    encoded_phenotypes = encode_phenotypes(
        phenotypes,
        test_data["phenotype_groups"],
        test_data["hpo_network"],
        test_data["alt2prim"],
        k=1000,
    )
    assert sum(encoded_phenotypes[1]) == 2


def test_read_records_file(test_data):
    with pytest.raises(SystemExit) as se:
        read_records_file("notafilepath/notafile")

    assert se.type == SystemExit
    assert se.value.code == 1

    records_truth = [
        {
            "sample": "118200",
            "age": 9.0,
            "gender": "Female",
            "terms": "HP:0001263|HP:0001251|HP:0001290|HP:0004322".split("|"),
        },
        {
            "sample": "118210",
            "age": 4.0,
            "gender": None,
            "terms": "HP:0001249|HP:0001263|HP:0001290".split("|"),
        },
        {
            "sample": "118211",
            "age": None,
            "gender": None,
            "terms": "HP:0001249|HP:0001263|HP:0001290".split("|"),
        },
    ]
    records_path = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "data/test.score-short.txt"
    )
    records = read_records_file(records_path, no_parents=False)
    assert records == records_truth
