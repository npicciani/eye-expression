# rule star_index:
#     input:
#         transcriptomePath=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"])        
#     output:
#         directory("results/star/index/")
#     threads: 8
#     log:
#         logfile=expand("logs/{transcriptome}_longestORFperGene.fasta_star_index.log", transcriptome=config["reference"])
#     shell:
#         "STAR --runMode genomeGenerate --runThreadN {threads} "
#         "--genomeDir 'results/star/index' "
#         "--genomeFastaFiles '{input.transcriptomePath}' "
#         "--genomeChrBinNbits 12 "
#         "--genomeSAindexNbases 11 "
#         "> {log.logfile}"

rule star_align:
    input:
        index="results/star/index",
        R1= lambda wildcards: "resources/rawdata/{sample}_{read}_1.fq.gz" for sample in config["samples"][wildcards.sample],
        R2="resources/rawdata/{sample}_{read}_2.fq.gz"
    output:
       "results/star/mapping/{sample}_{read}.aligned.out.sam"
    threads: 8
    # log:
    #    expand("logs/{sample}_star_align.log", sample=samples.loc[:,"sample_name"])
    shell:
        "STAR --runThreadN {threads} "
        "--genomeDir {input.index} "
        "--readFilesCommand gunzip -c "
        "--outFileNamePrefix results/star/mapping/{sample}.aligned "
        "--readFilesIn {input.R1} {input.R2}"
