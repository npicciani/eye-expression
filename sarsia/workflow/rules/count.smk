rule count_features:
    input:
        gtf=expand("results/reference/{transcriptome}_longestORFperGene.fasta.eggnog.gtf", transcriptome=config["reference"]["fileStem"]),
        reference=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["fileStem"]),
        mappedReads=expand("results/star/mapping/{sample}/Aligned.out.sam", sample=samples.loc[:,"sample_name"])
    output:
        "results/counts/counts_matrix.txt"
    threads: 20
    log:
        "logs/count_features/count_features.log"
    shell:
        "featureCounts -a {input.gtf} "
        "-o {output} -t exon -g gene_id "
        "-G {input.reference} "
        "-T {threads} {input.mappedReads}"
