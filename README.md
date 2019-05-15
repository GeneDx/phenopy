[![Build Status](http://gaithersburg-ds:8111/guestAuth/app/rest/builds/buildType:(id:Phenosim_Build)/statusIcon)](http://gaithersburg-ds:8111/viewType.html?buildTypeId=Phenosim_Build)
[![Coverage](http://gaithersburg-ds:8111/guestAuth/app/rest/builds/buildType:(id:Phenosim_Build)/artifacts/content/htmlcov/coverage.svg)](http://gaithersburg-ds:8111/guestAuth/repository/download/Phenosim_Build/.lastFinished/htmlcov/index.html)

# Phenosim
Phenotype similarity scoring using Python. Phenosim is a lightweight but highly optimized command line tool to robustly and efficiently perform semantic similarity scoring on generic entities with phenotype annotations.

## Installation
Clone from the bitbucket repo, and install using pip:
```bash
git clone ssh://git@bitbucket.bioreference.com:7999/gv/phenosim.git
cd phenosim
pip install .
```

## Usage
Phenosim is primarily used as a command line tool. An entity, as described here, is presented as a sample, gene, or disease, but could be any concept that warrants annotation of phenotype terms.   
There are several different use cases for phenosim:

1. Score similarity of an entity defined by the HPO terms from an input file against all the genes in `.phenosim/data/phenotypes_to_genes.txt`. We provide a test input file in the repo.
```bash
phenosim score tests/data/test.score.txt

query-A2M       0.0
query-AAGAB     0.0371
query-A4GALT    0.0025
```

2. Score similarity of an entity defined by the HPO terms from an input file against a custom list of entities with HPO annotations, referred to as the `--records-file`.
```bash
phenosim score tests/data/test.score.txt --records-file tests/data/test.score-product.txt

query-001	0.0809
query-002	0.0813
query-003	0.0781
```

3. Score pairwise similarity of entities defined in the `--records-file`.
```bash
phenosim score-product tests/data/test.score-product.txt

001-001	0.8444
001-002	0.1097
001-003	0.1211
```

## Parameters
Anny of the scoring Phenosim uses [`multiprocessing`](https://docs.python.org/3.4/library/multiprocessing.html?highlight=process) to parallelize computations.  

**To use 4 threads add to any scoring function:**  
`--threads 4`

**To score HPO term lists with parent terms pruned from the list:**
`--no_parents=True`

## Generate no_parents phenotypes_to_genes.txt file
From `phenosim/util.py` use `export_pheno2genes_with_no_parents` to output a version of the `phenotypes_to_genes.txt` which has parent terms pruned from the annotations.


## Under the hood
Phenosim will create a `.phenosim` directory in your home directory. This has data files from the [Human Phenotype Ontology](https://hpo.jax.org/app/) and
```bash
./phenosim/
  data/
    hp.obo
    hpo_network.pickle
    phenotype_to_genes.txt
  log.txt
```

## Hybrid Relative Specificity Similarity (HRSS)
The underlying algorithm which determines the semantic similarity for any two HPO terms is based on an implementation of HRSS, [published here](https://www.ncbi.nlm.nih.gov/pubmed/23741529).
