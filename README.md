# Scale-free network analysis

This repo contains the code for analyzing scale-free patterns in networks as described in this paper(https://www.nature.com/articles/s41467-019-08746-5). The data sets used in the paper are in the *degreesequences* directory.

## Citation information

If you use this code in your research, please cite this repo and/or the paper (linked above):

 Anna D. Broido & Aaron Clauset, Scale-free networks are rare, Nature Communications **10**, 1017 (2019).

## Usage

There are two ways to use this repo:

1) extract degree sequences from network data in the form of gml files and then analyze these for scale-free patterns

2) given degree sequences in the appropriate format (examples in the degree sequences folder), sort them into scale-free categories. Without gml information, this version of the pipeline treats each degree sequences as belonging to a unique network.  

### Dependencies
The code is written in Python2 and will not work in Python3. Additionally, the following packages must also be installed:
* NumPy
* Pandas
* SciPy
* mpmath
* python-igraph

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## A simple usage example
### With GML files
```python
import sys
import pandas as pd
sys.path.append('../code/')
import sfanalysis as sf

# location of gml files to analyze
gml_dir = 'gmls/'
# location to write degree sequences
deg_dir = 'degseqs/'
# make catalog of gmls and write degree sequence files
# each row of deg_df is a degree sequence file
deg_df = sf.write_degree_sequences(gml_dir, deg_dir)
# analyze all degree sequences (this will take a while for many or large data sets)
analysis_df = sf.analyze_degree_sequences(deg_dir, deg_df)
# categorize networks (by unique gml file) into scale-free categories
hyps_df = sf.categorize_networks(analysis_df)
```

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```



## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
