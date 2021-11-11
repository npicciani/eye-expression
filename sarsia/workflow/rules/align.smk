rule star_index:
    input:
        transcriptomePath=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["fileStem"])
    output:
        directory("results/star/index")
    threads: 8
    log:
        logfile=expand("logs/star_index/{transcriptome}_longestORFperGene.fasta_star_index.log", transcriptome=config["reference"]["fileStem"])
    conda:
        "../../workflow/envs/star.yaml" #star v2.7.9
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
       "logs/star_pe_multi/{sample}_star_align.log"
    params:
        index=lambda wc, input: input.index
    threads: 20
    wrapper:
        "v0.75.0/bio/star/align"
