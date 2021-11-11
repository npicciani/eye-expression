rule get_opsins:
    input:
        script="workflow/scripts/pia.py",
        peptides=expand("results/reference/{transcriptome}_longestORFperGene.pep", transcriptome=config["reference"]["fileStem"])
    output:
        expand("results/opsins/{transcriptome}_longestORFperGene.opsins.fasta", transcriptome=config["reference"]["fileStem"])
    conda:
        "../envs/pia.yaml"
    threads: 8
    params:
        outdir="$PWD/results/opsins",
        bait="resources/pia/baits1217.fasta",
        alignment="resources/pia/all_0512_gt1_rs50_s65.aln",
        tree="resources/pia/rep17.treefile"
    shell:
        "python {input.script} {input.peptides} {params.outdir} {params.bait} {params.alignment} {params.tree}"
