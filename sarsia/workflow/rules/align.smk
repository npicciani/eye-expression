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

# rule star_align:
#     input:
#         index="results/star/index",
#         R1= "resources/rawdata/{sample}_1.fq.gz", sample=samples.loc[:,"sample_name"],
#         R2="resources/rawdata/{sample}_2.fq.gz", sample=samples.loc[:,"sample_name"]
#     output:
#         expand("results/star/mapping/{sample}.aligned.out.sam", sample=samples.loc[:,"sample_name"])
#     threads: 8
#     # log:
#     #    expand("logs/{sample}_star_align.log", sample=samples.loc[:,"sample_name"])
#     shell:
#         "STAR --runThreadN {threads} "
#         "--genomeDir {input.index} "
#         "--readFilesCommand gunzip -c "
#         "--outFileNamePrefix results/star/mapping/{wildcards.sample}.aligned "
#         "--readFilesIn {input.R1} {input.R2}"

rule star_pe_multi:
    input:
        fq1=get_reads_R1,
        fq2=get_reads_R2
    output:
        "results/star/mapping/{sample}/Aligned.out.sam"
    log:
       "logs/{sample}_star_align.log"
    params:
        index="results/star/index"
    threads: 8
    wrapper:
        "v0.75.0/bio/star/align"