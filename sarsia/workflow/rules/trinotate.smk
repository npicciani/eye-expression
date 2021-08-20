rule trinotate:
    """
    Functionally annotate reference assembly using Trinotate.
    """
    input:
        fastafile=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["filename"]),
        peptidefile=expand("results/reference/{transcriptome}_longestORFperGene.pep.transdecoder_dir/longest_orfs.pep", transcriptome=config["reference"]["filename"]),
        transdecoder_peptidefile=expand("results/reference/{transcriptome}.transdecoder_dir/longest_orfs.pep", transcriptome=config["reference"]["filename"]),
        geneTranscript_map=expand("results/reference/{transcriptome}_longestORFperGene.fasta.geneID_to_transcript.txt", transcriptome=config["reference"]["filename"]),
        script="resources/trinotate.sh"
    output:
        outdir="results/trinotate"
    conda:
        "../envs/trinotate.yaml"
    shell:
        "bash {input.script} {input.fastafile} {input.peptidefile} {input.transdecoder_peptidefile} {input.geneTranscript_map} {output.outdir}"
