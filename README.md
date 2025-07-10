# eye-expression
Gene expression analysis of jellyfish eyes

This repository contains data and scripts associated with Picciani et al. 2025.

Brief description of directories:

For each species (Sarsia, Aurelia, Tripedalia) you can find the following:

	+/rawdata: raw sequencing data from whole-body jellyfish (Sarsia tubulosa, EMBL) and from tissue-specific samples (rhopalia or tentacle bulbs, manubrium, tentacles)

	+/scripts: custom scripts for data analysis

	+/jobs: slurm jobs for data analysis ran on Farnam

	+/results: results from each analysis step

The Expression_analyses folder contains the R code and primary data files for differential expression, GO enrichment, and cross-species comparisons.
