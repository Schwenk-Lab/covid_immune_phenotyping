[![DOI](https://zenodo.org/badge/979315666.svg)](https://doi.org/10.5281/zenodo.15356855)

# Immune phenotyping of self-sampling citizens

Analysis and visualisation code accompanying *Immune phenotyping of self-sampling citizens*. Analysis was done using R 4.3.2.

`_targets.R`  
The analysis pipeline is written using the `targets` R package. This file defines the pipeline steps. Data is read from a `../data/` folder and results are saved in `../results/` containing subfolders for different parts of the anlaysis. Data that will be published is already preprocessed, which is why data wrangling steps in the beginning are commented out in the pipeline. Some steps in the downstream analysis may need to be skipped or the input data substituted if they use non-processed or non-adjusted data.

`target_params.R`  
Some parameters used in the pipeline.

`functions_....R`  
Scripts containing functions used for pipeline steps, divided into scripts for e.g. handling data wrangling, serology data, proteomics data.

`misc_functions.R`  
Various smaller functions used throughout the pipeline.

`rmd/..._report.Rmd`  
Scripts generating reports for showing or saving certain visualisations.

The remaining scripts contain single larger functions such as `comp_groups.R` for group comparison plots and `tidy_classification.R` for classification using Lasso. `ben_mixmod.jl` contains Julia code for normalisation using a mixed model written by [Ben Murrell](https://github.com/murrellb) ([Murrell group](https://github.com/MurrellGroup), Karolinska Institutet). This script is run in using a conda environment called `julia_env` containing Julia 1.7.2. 
