# load sample information from sample.tsv file (necessary??)
import pandas as pd
import os

samples = (
    pd.read_csv(config["samples"], delim_whitespace=True, dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
    )

units = (
    pd.read_csv(config["units"], delim_whitespace=True, dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))

def get_reads_R1(wildcards):
    if config["mergeReads"]["activate"]:
        sample=wildcards.sample
        return "results/merged/{sample}_1.fq.gz"
    if config["sra"]["activate"]:
        accession=wildcards.accession
        return "resources/rawdata/{accession}_1.fastq"
    sample_units = units.loc[wildcards.sample]
    return sample_units["fq1"]

def get_reads_R2(wildcards):
    if config["mergeReads"]["activate"]:
        sample=wildcards.sample
        return "results/merged/{sample}_2.fq.gz"
    if config["sra"]["activate"]:
        accession=wildcards.accession
        return "resources/rawdata/{accession}_2.fastq"
    sample_units = units.loc[wildcards.sample]
    return sample_units["fq2"]

def get_sam(wildcards):
    resultsPath="results/star/mapping"
    accession=wildcards.accession
    samFilename="Aligned.out.sam"
    return resultsPath + "/" + accession + "/" + samFilename

def get_fastqs(wildcards):
    fq= "fq{}".format(wildcards.read[-1]) #-1 is the very last string character; gives you fq1 and fq2
    return units.loc[wildcards.sample, fq].tolist()

def get_sra(wildcards):
    accession=wildcards.accession
    return accession
