rule deseq2_init:
    input:
        counts="results/counts/counts_matrix.txt"
    output:
        "results/deseq2/dds.rds",
        "results/deseq2/counts_normalized.csv"
    params:
        samples=config["samples"],
        model=config["diffexp"]["model"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log"
    script:
        "../scripts/deseq2-init.R"
