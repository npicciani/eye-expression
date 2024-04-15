rule get_SRA_fastqs:
    output:
        "resources/rawdata/{accession}_1.fastq",
        "resources/rawdata/{accession}_2.fastq"
    params:
        get_sra
    conda:
        "../../workflow/envs/sra-tools.yaml"
    shell:
        "fasterq-dump --split-files --outdir resources/rawdata {params}"

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