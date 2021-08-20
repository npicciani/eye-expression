rule star_index:
    input:
        reference=config["reference"]
    output:
        directory("results/star/index")
    threads: 8
    log:
        logfile="logs/star_index/star_index.log"
    shell:
        "STAR --runMode genomeGenerate --runThreadN {threads} "
        "--genomeDir 'results/star/index' "
        "--genomeFastaFiles '{input.reference}' "
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
        index=lambda wc, input: input.index,
        extra="--clip3pNbases 30 30 --clip5pNbases 30 30"
    threads: 20
    wrapper:
        "v0.75.0/bio/star/align"
