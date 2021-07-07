rule merge_fastqs:
    input:
        get_fastqs
    output:
        "results/merged/{sample}{read}.fq.gz"
    log:
        "logs/merge_fastqs/{sample}{read}.log"
    wildcard_constraints:
        read="_1|_2" #could use other patterns depending on fastq filenames
    shell:
        "cat {input} > {output} 2> {log}"