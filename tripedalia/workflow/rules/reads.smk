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

rule SRA_fastqs:
    input:
        get_sra
    output:
        "resources/rawdata/{accession}_1.fastq",
        "resources/rawdata/{accession}_2.fastq"
    conda:
        "../../workflow/envs/sra-tools.yaml"
 #   params:
 #       accessions=expand("{accs}", accs=units.loc[:,"sra"]) #dump all files at once
    shell:
#        "fastq-dump --split-files --outdir resources/rawdata {params.accessions}"
        "fastq-dump --split-files --outdir resources/rawdata {input}"
