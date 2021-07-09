rule busco_scores:
    input:
        expand("results/reference/{transcriptome}.fixed_longestORFperGene.fasta", transcriptome=config["reference"]["filename"])
    output:
        path=directory("results/busco")
    threads: 20
    conda: 
        "../envs/busco.yaml" #busco=5.1.3
    log:
        expand("logs/busco_scores/{transcriptome}.fixed_longestORFperGene.fasta_busco_scores.log", transcriptome=config["reference"]["filename"])
    params:
        download_path="results/busco/busco_downloads",
        mode="transcriptome",
        lineage="metazoa",
        stem=expand("{transcriptome}.fixed_longestORFperGene.fasta", transcriptome=config["reference"]["filename"])
    shell:
        "busco -i {input} -o {params.stem} --force --out_path {output.path} -l {params.lineage} -m {params.mode} --download_path {params.download_path} -c {threads}"

rule mapping_stats:
    input:
        get_sam
    threads: 8
    conda: 
        "../envs/samtools.yaml" #samtools=1.12
    output:
        "results/star/mapping/{accession}/Aligned.out.sam.stats"
    shell:
        "samtools stats --threads {threads} {input} > {output}"