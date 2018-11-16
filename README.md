2018koyeast
===========

Analysis on chemical genomics data across 4 yeast strains.

Introduction
------------

This repository contains the reproducible pipeline that starts from the
chemical genomics data to the final summary plots. Part of the data is generated
prior to this repository (e.g. the s-scores and the gene disruption score for
the yeast natural isolates). The pipeline is based on a [snakemake](https://snakemake.readthedocs.io/),
with the exception of the [jupyter notebooks](https://jupyter.org/), which have to be run one by one.
Most of the steps involve a script present in the `src` directory, which are mostly based on python 3
and the following libraries:

* `numpy`
* `scipy`
* `pandas`
* `scikit-learn`
* `statsmodels`
* `biopython`
* `limix`
* `xarray`
* `DendroPy`

If jupyter notebooks are used, the `matplotlib` and `seaborn` python libraries are also required.

Usage
-----

Provided that all input files are present the whole pipeline can be run by typing: `snakemake -p all`.
You can add `--cores XX` to parallelize some steps. Please make sure you have at least 12Gb of RAM available.

The required input files are:

* 
*
*


