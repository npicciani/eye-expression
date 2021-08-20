rule mapping_stats:
    input:
        get_sam
    threads: 8
    conda:
        "../envs/samtools.yaml" #samtools=1.12
    output:
        "results/star/mapping/{sample}/Aligned.out.sam.stats"
    shell:
        "samtools stats --threads {threads} {input} > {output}"
