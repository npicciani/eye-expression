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

# def is_paired_end(sample):
#     sample_units = units.loc[sample]
#     fq2_null = sample_units["fq2"].isnull()
#     sra_null = sample_units["sra"].isnull()
#     paired = ~fq2_null | ~sra_null
#     all_paired = paired.all()
#     all_single = (~paired).all()
#     assert (
#         all_single or all_paired
#     ), "invalid units for sample {}, must be all paired end or all single end".format(
#         sample
#     )
#     return all_paired

def get_reads_R1(wildcards):
    # if is_paired_end(wildcards.sample):
        # unit = units.loc[wildcards.sample]
    if config["mergeReads"]["activate"]:
        sample=wildcards.sample
        return "results/merged/{sample}_1.fq.gz"
    sample_units = units.loc[wildcards.sample]
    return sample_units["fq1"]

def get_reads_R2(wildcards):
    # if is_paired_end(wildcards.sample):
        # unit = units.loc[wildcards.sample]
    if config["mergeReads"]["activate"]:
        sample=wildcards.sample
        return "results/merged/{sample}_2.fq.gz"
    sample_units = units.loc[wildcards.sample]
    return sample_units["fq2"]

def get_sam(wildcards):
    resultsPath="results/star/mapping"
    sampleID=(wildcards.sample)
    samFilename="Aligned.out.sam"
    samFilepath=os.path.join(resultsPath, sampleID, samFilename)
    return samFilepath

def get_fastqs(wildcards):
    fq= "fq{}".format(wildcards.read[-1]) #-1 is the very last string character; gives you fq1 and fq2
    return units.loc[wildcards.sample, fq].tolist()
