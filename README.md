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

4. Perform a grid search to determine the ideal number of clusters via hierarchical clustering. The input of this function is the output of `score-product` and display performance metrics for up to `k` clusters where `--max-clusters k`. The output columns are `[k, linkage type, average silhouette, geometric mean of the cluster size, cluster sizes]`
```bash
phenosim cluster-grid-search score-product.results --max-clusters 3

2	complete	-0.1642	49.3558	[58, 42]
3	complete	-0.2312	26.4764	[58, 32, 10]
2	average	0.0686	9.9498	[99, 1]
3	average	0.0237	4.6104	[98, 1, 1]
2	single	-0.0787	9.9498 [99, 1]
3	single	-0.1363	4.6104 [98, 1, 1]
```

5. Assign cluster membership to entities on entities in the results file from `score-product`. Use `--linkage` and one of `{single, average, complete}` to select the type of linkage algorithm. Use `--k` and `<int>` to specify the number of clusters.

```bash
phenosim cluster-assign score-product.results --linkage complete --k 4

001	3
002	2
003	3
004	0
005	0
```

## Parameters
Anny of the scoring Phenosim uses [`multiprocessing`](https://docs.python.org/3.4/library/multiprocessing.html?highlight=process) to parallelize computations.  
To use 4 threads add to any scoring function:  
`--threads 4`

## Under the hood
Phenosim will create a `.phenosim` directory in your home directory. This has data files from the [Human Phenotype Ontology]() and
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
