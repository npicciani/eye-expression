# Gene tree-based annotations

This is a snakemake workflow to generate gene trees with functionally annotated sequences that is divided into a few main steps:

1. Download of public data (proteins from target species) for building gene trees.
2. Functionally annotate each protein set using eggNOG-mapper v.2.1.10.
3. Generate gene trees using orthofinder v.2.5.4.
4. Append gene names from functional annotations to sequence headers in gene trees.

## How to get started

This workflow uses conda environments for every step. That means that you don't need to manually install any program. 
The programs necessary for running this workflow are defined in separate yaml files, which will be used by snakemake
in the very beginning to create all conda environments that will be used in snakemake rules.

1. The first thing to do here is to install the snakemake_base environment which contains snakemake and other basic programs. 
In McCleary, you can create a conda environment named `snakemake` by running this from any folder:

```
ml miniconda
conda env create -n snakemake --file snakemake_base.yml
```

2. You can git clone this repo:

```
git clone https://github.com/npicciani/nanomia.git
```

3. In `config/download_targets.tsv`, include species name and sources for taxa you would like to include in your gene trees.
4. Run snakemake. I have set up a cluster profile named slurm so my command looks like:
```
conda activate snakemake
snakemake --profile slurm
```

## Output files
The main output file is `results/annotations/gene_trees.master.annotated.txt`. In this file, each line corresponds to a gene tree in parenthetical notation and sequence headers with gene names.

You can also check the master annotation file (`results/annotations/annotations.master.txt`), which concatenates eggnog annotations for all sequences from all species. It includes additional information for each sequence, such as KEGG pathways, Gene Ontology labels, EC numbers, BiGG reactions, CAZy terms, etc. - check out the eggnog publication [here](https://academic.oup.com/mbe/article/38/12/5825/6379734).
