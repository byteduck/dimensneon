# dimenSNEon - CSE 185
This python library that aims to act as a drop-in replacement for scanpy's t-SNE simulation functionality.

It is mainly intended for use in 10x scRNA-seq workflows, although it will theoretically work with any compatible AnnData structure.

## Installation instructions

Installation requires the `scanpy` library to be installed (and `leidenalg` as well for tests). You can install these using the following command (However, `setup.py` should install all required dependencies for you):

```shell
pip install scanpy leidenalg
```

You can then install the library by running:

```shell
python setup.py install
```

## Basic usage

dimenSNEon is a drop-in replacement for scanpy's t-SNE simulation functionality. You can import the library in your script by using:

```python
import dimensneon.tsne as dtsne
```

And then, once you have run PCA on your data (as you would before using scanpy's `sc.tl.tsne` function), you can calculate t-SNE embeddings for your data by using:

```python
# data is our AnnData object
dtsne.tsne(data)
```

Then, you can plot the results using scanpy's `sc.pl.tsne` as you normally would.

## Advanced usage

dimenSNEon allows you to tweak aspects of the simulation by supplying optional arguments after the AnnData object. These arguments are:

- `perplexity` (`float`, default = `30.0`): The perplexity value to target for the t-SNE simulation.
- `iterations` (`int`, default = `1000`): The number of iterations of the simulation to run.

## Credits

This repository was created by Aaron Sonin.

The mathematics and mechanics behind t-SNE are from the paper "Visualizing Data using t-SNE" by Laurens van der Maaten and Geoffrey Hinton. Citation: 
```
van der Maaten, Laurens & Hinton, Geoffrey. (2008). Viualizing data using t-SNE. Journal of Machine Learning Research. 9. 2579-2605. 
```

Additional inspiration was drawn from sklearn and scanpy's t-SNE functionality.

## Sample Data

The sample data downloaded by the tests in this repository are as follows:

- [10x scRNA sample data v3.0.0](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3?)
- [10x scRNA data from Functional, metabolic and transcriptional maturation of human pancreatic islets derived from stem cells](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167880)