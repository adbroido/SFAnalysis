# Scale-free network analysis

This repo contains the code for analyzing scale-free patterns in networks as described in this paper(https://www.nature.com/articles/s41467-019-08746-5). The data sets used in the paper are in the *degreesequences* directory.

## Citation information

If you use this code in your research, please cite this repo and/or the paper linked above:

 Anna D. Broido & Aaron Clauset, Scale-free networks are rare, Nature Communications **10**, 1017 (2019).

## Usage

This repo assumes you have network data in the form of gml files. If you already have degree sequences, you may have to slightly alter the code and pipeline described below.

### Dependencies
The code is written in Python2. It will not work in Python3. Additionally, the following packages must also be installed:
* NumPy
* Pandas
* SciPy
* igraph

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
```python
import sfanalysis as sf
import pandas as pd

# location of gml files to analyze
gml_dir = '/filepath/'
# location to write degree sequences
deg_dir = '/filepath/'
# make catalog of gmls and write degree sequence files
# each row of deg_df is a degree sequence file
deg_df = sf.write_degree_sequences(gml_dir, deg_dir)
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
