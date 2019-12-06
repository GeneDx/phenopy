[![python-version](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![github-actions](https://github.com/GeneDx/phenopy/workflows/Python%20package/badge.svg)](https://github.com/GeneDx/phenopy/actions)
[![codecov](https://codecov.io/gh/GeneDx/phenopy/branch/develop/graph/badge.svg)](https://codecov.io/gh/GeneDx/phenopy)
[![DOI](https://zenodo.org/badge/207335538.svg)](https://zenodo.org/badge/latestdoi/207335538)

# phenopy
`phenopy` is a Python package to perform phenotype similarity scoring by semantic similarity. `phenopy` is a
lightweight but highly optimized command line tool and library to efficiently perform semantic similarity scoring on
generic entities with phenotype annotations from the [Human Phenotype Ontology (HPO)](https://hpo.jax.org/app/).

![Phenotype Similarity Clustering](https://raw.githubusercontent.com/GeneDx/phenopy/develop/notebooks/output/cluster_three_diseases.png)

## Installation
Install using pip:
```bash
pip install phenopy
```

Install from GitHub:
```bash
git clone https://github.com/GeneDx/phenopy.git
cd phenopy
python setup.py install
```

## Command Line Usage
### Initial setup
phenopy is designed to run with minimal setup from the user, to run phenopy with default parameters (recommended), skip ahead
to the [Commands overview](#Commands-overview).  

This section provides details about where phenopy stores data resources and config files. The following occurs when
you run phenopy for the first time.
 1. phenopy creates a `.phenopy/` directory in your home folder and downloads external resources from HPO into the
  `$HOME/.phenopy/data/` directory.
 2. phenopy creates a `$HOME/.phenopy/phenopy.ini` config file where users can set variables for phenopy to use
 at runtime.

### Commands overview
`phenopy` is primarily used as a command line tool. An entity, as described here, is presented as a sample, gene, or
disease, but could be any concept that warrants annotation of phenotype terms. 

Use `phenopy score` to perform semantic similarity scoring in various formats. Write the results of any command to file 
using `--output-file=/path/to/output_file.txt`

1. Score similarity of entities defined by the HPO terms from an input file against all the OMIM diseases in
    `.phenopy/data/phenotype.hpoa`. We provide a test input file in the repo. The default summarization method is to
     use `--summarization-method=BMWA` which weighs each diseases' phenotypes by the frequency of a phenotype seen in 
     each particular disease.
    ```bash
    phenopy score tests/data/test.score.txt  
    ```
    Output:
    ```
    #query	entity_id	score
    118200  210100  0.0
    118200  615779  0.0
    118200  613266  0.0052
    ...
    ```
    
2. Score similarity of entities defined by the HPO terms from an input file against all the OMIM diseases in
    `.phenopy/data/phenotype.hpoa`, to use the non-weighted summarization method use `--summarization-method=BMA` which
    uses a traditional *best-match average* summarization of semantic similarity scores when comparing terms from record *a* 
    with terms from record *b*.
    ```bash
    phenopy score tests/data/test.score.txt --summarization-method=BMWA
    ```
    Output:
    ```
    #query	entity_id	score
    118200  210100  0.0
    118200  615779  0.0
    118200  613266  0.0052
    ...
    ```

3. Score similarity of an entities defined by the HPO terms from an input file against a custom list of entities with HPO annotations, referred to as the `--records-file`. Both files are in the same format.
    ```bash
    phenopy score tests/data/test.score-short.txt --records-file tests/data/test.score-long.txt
    ```
    Output:
    ```
    #query  entity_id       score
    118200  118200  0.0169
    118200  300905  0.0156
    118200  601098  0.0171
    ...
    ```

4. Score pairwise similarity of entities defined by the HPO terms from an input file using `--self`.

    ```bash
    phenopy score tests/data/test.score-long.txt --threads 4 --self
    ```
    Output:
    ```
    #query  entity_id       score
    118200  118200  0.2284
    118200  118210  0.1302
    118200  118211  0.1302
    118210  118210  0.2048
    118210  118211  0.2048
    118211  118211  0.2048
    ```
5. Score age-adjusted pairwise similarity of entities defined in the input file, 
    using phenotype mean age and standard deviation defined in the `--ages_distribution_file`,
    select best-match weighted average as the scoring summarization method `--summarization-method BMWA`.

    ```bash
    phenopy score tests/data/test.score-short.txt --ages_distribution_file tests/data/phenotype_age.tsv --summarization-method BMWA --threads 4 --self
    ```
    Output:
    ```
    #query  entity_id       score
    118200  210100  0.0
    118200  177650  0.0127
    118200  241520  0.0
    ...
    ```
    
    The phenotype age file contains hpo-id, mean, sd as tab separated text as follows
    
    |  |  | |
    |------------|------|-----|
    | HP:0001251 | 6.0  | 3.0 |
    | HP:0001263 | 1.0  | 1.0 |
    | HP:0001290 | 1.0  | 1.0 |
    | HP:0004322 | 10.0 | 3.0 |
    | HP:0001249 | 6.0  | 3.0 |

  If no phenotype ages file is provided `--summarization-method=BMWA` can be selected to use default, open access literature-derived phenotype ages (~ 1,400 age weighted phenotypes).  
   ```bash
    phenopy score tests/data/test.score-short.txt  --summarization-method BMWA --threads 4
   ```


## Parameters
For a full list of command arguments use `phenopy [subcommand] --help`:
```bash
phenopy score --help
```
Output:
```
    --output_file=OUTPUT_FILE
        File path where to store the results. [default: - (stdout)]
    --records_file=RECORDS_FILE
        An entity-to-phenotype annotation file in the same format as "input_file". This file, if provided, is used to score entries in the "input_file" against entries here. [default: None]
    --annotations_file=ANNOTATIONS_FILE
        An entity-to-phenotype annotation file in the same format as "input_file". This file, if provided, is used to add information content to the network. [default: None]
    --ages_distribution_file=AGES_DISTRIBUTION_FILE
        Phenotypes age summary stats file containing phenotype HPO id, mean_age, and std. [default: None]
    --self=SELF
        Score entries in the "input_file" against itself.
    --summarization_method=SUMMARIZATION_METHOD
        The method used to summarize the HRSS matrix. Supported Values are best match average (BMA), best match weighted average (BMWA), and maximum (maximum). [default: BMWA]
    --threads=THREADS
        Number of parallel processes to use. [default: 1]
```
## Library Usage
The `phenopy` library can be used as a `Python` module, allowing more control for advanced users.   

```python
import os
from phenopy import generate_annotated_hpo_network
from phenopy.score import Scorer

# data directory
phenopy_data_directory = os.path.join(os.getenv('HOME'), '.phenopy/data')

# files used in building the annotated HPO network
obo_file = os.path.join(phenopy_data_directory, 'hp.obo')
disease_to_phenotype_file = os.path.join(phenopy_data_directory, 'phenotype.hpoa')

# if you have a custom ages_distribution_file, you can set it here.
ages_distribution_file = os.path.join(phenopy_data_directory, 'xa_age_stats_oct052019.tsv')

hpo_network, alt2prim, disease_records = \
    generate_annotated_hpo_network(obo_file,
                                   disease_to_phenotype_file,
                                   ages_distribution_file=ages_distribution_file
                                   )

scorer = Scorer(hpo_network)

terms_a = ['HP:0001882', 'HP:0011839']
terms_b = ['HP:0001263', 'HP:0000252']

print(scorer.score(terms_a, terms_b))
```
Output:
```
0.0005
```

The library can be used to prune parent phenotypes from the `phenotype.hpoa` and store pruned annotations as a file.
```python
from phenopy.util import export_phenotype_hpoa_with_no_parents
# saves a new file of phenotype disease annotations with parent HPO terms removed from phenotype lists.
disease_to_phenotype_no_parents_file = os.path.join(phenopy_data_directory, 'phenotype.noparents.hpoa') 
export_phenotype_hpoa_with_no_parents(disease_to_phenotype_file, disease_to_phenotype_no_parents_file, hpo_network)
```

### Config
While we recommend using the default settings for most users, the config file *can be* modified: `$HOME/.phenopy/phenopy.ini`.

To run phenopy with a different version of `hp.obo`, set the path of `obo_file` in `$HOME/.phenopy/phenopy.ini`.

## Contributing
We welcome contributions from the community. Please follow these steps to setup a local development environment.  
```bash
pipenv install --dev
```

To run tests locally:
```bash
pipenv shell
coverage run --source=. -m unittest discover --start-directory tests/
coverage report -m
```  

## References
The underlying algorithm which determines the semantic similarity for any two HPO terms is based on an implementation of HRSS, [published here](https://www.ncbi.nlm.nih.gov/pubmed/23741529).


## Citing Phenopy
Please use the following Bibtex to cite this software.
```
@software{arvai_phenopy_2019,
    title = {Phenopy},
    rights = {Attribution-NonCommercial-ShareAlike 4.0 International},
    url = {https://github.com/GeneDx/phenopy},
    abstract = {Phenopy is a Python package to perform phenotype similarity scoring by semantic similarity. 
        Phenopy is a lightweight but highly optimized command line tool and library to efficiently perform semantic 
        similarity scoring on generic entities with phenotype annotations from the Human Phenotype Ontology (HPO).},
    version = {0.3.0},
    author = {Arvai, Kevin and Borroto, Carlos and Gainullin, Vladimir and Retterer, Kyle},
    date = {2019-11-05},
    year = {2019},
    doi = {10.5281/zenodo.3529569}
}
```
