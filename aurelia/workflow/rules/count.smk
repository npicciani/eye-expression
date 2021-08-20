rule count_features:
    input:
        gtf=config["reference_gtf"],
        reference=config["reference"],
        mappedReads=expand("results/star/mapping/{sample}/Aligned.out.sam", sample=samples.loc[:,"sample_name"])
    output:
        "results/counts/counts_matrix.txt"
    threads: 20
    log:
        "logs/count_features.log"
    shell:
        "featureCounts -a {input.gtf} "
        "-o {output} -t exon -g gene_id "
        "-G {input.reference} "
        "-T {threads} {input.mappedReads}"
