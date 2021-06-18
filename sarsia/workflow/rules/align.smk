import pandas as pd

samples = (
    pd.read_csv(config["samples"], delim_whitespace=True, dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
    )
def get_reads_input_R1(wildcards):
    return config["samples"][wildcards.sample]

def get_reads_input_R2(wildcards):
    return config["samples"][wildcards.sample]


rule star_index:
    input:
        transcriptome=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"])        
    output:
        "results/star/index"
    threads: 12
    log:
        "logs/{input.transcriptome}_star_index.log"
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} "
        "--genomeDir 'results/star/index' "
        "--genomeFastaFiles 'results/reference/{input.transcriptome}' "
        "--genomeChrBinNbits 12 "
        "--genomeSAindexNbases 11"

rule star_align:
    input:
        index="results/star/index",
        R1=get_reads_input_R1,
        R2=get_reads_input_R2        
    output:
        "results/start/mapping/{sample}.aligned.out.sam"
    threads: 25
    log:
        "logs/{sample}_star_align.log"
    shell:
        "STAR --runThreadN {threads} "
        "--genomeDir {input.index} "
        "--readFilesCommand gunzip -c "
        "--outFileNamePrefix results/star/mapping/{sample}.aligned "
        "--readFilesIn {input.R1} {input.R2}"
