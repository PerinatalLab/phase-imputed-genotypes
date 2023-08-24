# Snakemake workflow: 

[![Snakemake](https://img.shields.io/badge/snakemake-≥{{cookiecutter.min_snakemake_version}}-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/{{cookiecutter.repo_name}}.svg?branch=master)](https://travis-ci.org/snakemake-workflows/{{cookiecutter.repo_name}})

This Snakemake workflow is to be used for phasing imputed genotype data from the Norwegian Mother, Father and Child Cohort study, which includes parent-offspring duos and trios, as well as non-related individuals.

## Authors

* Pol Sole-Navais (@psnavais) - pol.sole.navais@gu.se

## Dependencies  

* Global environment to run the Snakemake pipeline (not all required):  
+ python=3.8.10  
+ numpy=1.23.3  
+ pandas=1.5.0  
+ scipy=1.9.1  
+ bedtools=2.30.0  
+ bcftools=1.15.1  
+ mamba=0.25.0  
* Other dependencies can be found under the `workflow/envs/` folder (see below).  

## Structure  

All scripts needed to run this pipeline are contained in the `workflow` folder.  

* `Snakefile`: main Snakefile, calls other snakefiles and lists files to be collected  
* `envs`: contains the .yml environment files to run this pipeline  
* `rules`: contains the rules for each of the steps in the pipeline in different Snakefiles (*.smk*)  
* `scripts`: contains R and python scripts used by the pipeline  
* `report`: not used
* `schemas`: not used

## Usage  

Call snakemake from the main folder - it will automatically detect that the main `Snakefile` is in the `workflow` folder. Beware that the pipeline is computationally intensive (both memory and threads used) because it is run on imputed genetic data.  


