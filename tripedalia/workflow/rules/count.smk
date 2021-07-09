rule count_features:
    input:
        gtf=expand("results/reference/{transcriptome}.fixed_longestORFperGene.fasta.eggnog.gtf", transcriptome=config["reference"]["filename"]),
        reference=expand("results/reference/{transcriptome}.fixed_longestORFperGene.fasta", transcriptome=config["reference"]["filename"]),
        mappedReads=expand("results/star/mapping/{accession}/Aligned.out.sam", accession=units.loc[:,"sra"])
        
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