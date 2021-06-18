rule busco_scores:
    input:
        "results/reference/{transcriptome}_longestORFperGene.fasta"
    output:
        "results/busco/{transcriptome}.metazoa"
    threads: 20
    log:
        "logs/{transcriptome}_busco_scores.log"
    shell:
        "busco -i {input} -o {output} -l metazoa -m tran -c {threads}"

rule mapping_stats:
    input:
        "results/star/mapping/{sample}.aligned.out.sam"
    output:
        "results/star/mapping_stats/{sample}.samstats.tab"
    shell:
        "samtools stats {input} > {output}"
