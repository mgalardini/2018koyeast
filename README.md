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

Copyright
---------

Copyright 2018 EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
