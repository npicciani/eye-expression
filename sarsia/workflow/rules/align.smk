rule star_index:
    input:
        transcriptomePath=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"])       
    output:
        directory("results/star/index")
    threads: 8
    log:
        logfile=expand("logs/{transcriptome}_longestORFperGene.fasta_star_index.log", transcriptome=config["reference"])
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} "
        "--genomeDir 'results/star/index' "
        "--genomeFastaFiles '{input.transcriptomePath}' "
        "--genomeChrBinNbits 12 "
        "--genomeSAindexNbases 11 "
        "> {log.logfile}"

rule star_pe_multi:
    input:
        fq1=get_reads_R1,
        fq2=get_reads_R2,
        index="results/star/index" #include index as input so that snakemake checks that it exists before executing this rule
    output:
        "results/star/mapping/{sample}/Aligned.out.sam"
    log:
       "logs/{sample}_star_align.log"
    params:
        index=lambda wc, input: input.index
    threads: 20
    wrapper:
        "v0.75.0/bio/star/align"