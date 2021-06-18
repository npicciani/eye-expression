rule count_feature:
    input:
        gtf="results/reference/{transcriptome}_longestORFperGene.fasta.eggnog.gtf",
        reference="results/reference/{transcriptome}_longestORFperGene.fasta",
        mappedReads="results/star/mapping/{sample}.aligned.out.sam"
    output:
        "results/counts/count_matrix.txt"
    threads: 20
    log: 
        "logs/count_feature.log"
    shell:
        "featureCounts -a {input.gtf} "
        "-o {output} -t exon -g gene_id "
        "-G {input.reference} "
        "-T {threads} {input.mappedReads}"