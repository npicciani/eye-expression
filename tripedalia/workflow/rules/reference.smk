rule fix_TSA:
    """
    Edit the headers from TSA nucleotide fasta file.
    From: >GHAQ01000008.1 TSA: Tripedalia cystophora tri_comp79_c0_seq1, transcribed RNA sequence
    To: >comp79_c0_seq1
    """
    input:
        transcriptomePath=config["reference"]["path"],
        script="workflow/scripts/fastaClean.py"
    output:
        fixedTranscriptome=expand("results/reference/{transcriptome}.fixed", transcriptome=config["reference"]["filename"])
    shell:
        """
        mkdir results/reference
        python {input.script} {input.transcriptomePath} > {output.fixedTranscriptome}
        """

rule generate_longest_ORFs:
    """
    Infer open reading frames in reference transcriptome.
    """
    input:
        transcriptomePath=expand("results/reference/{transcriptome}.fixed", transcriptome=config["reference"]["filename"])
    output:
        expand("results/reference/{transcriptome}.fixed.transdecoder_dir/longest_orfs.pep", transcriptome=config["reference"]["filename"])
    params:
        reference=config["reference"]["filename"]
    conda:
        "../../workflow/envs/transdecoder.yaml" #transdecoder v5.5.0
    shell:
        "TransDecoder.LongOrfs -t {input.transcriptomePath} --output_dir results/reference/{params.reference}.fixed.transdecoder_dir"

rule keep_longest_ORF_per_gene:
    """
    Keep longest open reading frames per gene (as per gene identifier type defined in the config yaml).
    Details on functions and arguments in the python script itself.
    """
    input:
        longestORFs=expand("results/reference/{transcriptome}.fixed.transdecoder_dir/longest_orfs.pep", transcriptome=config["reference"]["filename"]),
        transcriptomePath=expand("results/reference/{transcriptome}.fixed", transcriptome=config["reference"]["filename"]),
        script="workflow/scripts/keepLongestORFperGene.py"
    output:
        peptides=expand("results/reference/{transcriptome}.fixed_longestORFperGene.pep", transcriptome=config["reference"]["filename"]),
        nucleotides=expand("results/reference/{transcriptome}.fixed_longestORFperGene.fasta", transcriptome=config["reference"]["filename"])
    params:
        geneID_type=config["geneIDType"]
    shell:
        "python {input.script} -p {input.longestORFs} -t "
        "{input.transcriptomePath} -identifier {params.geneID_type} -o results/reference"

rule make_GTF:
    input:
        nucleotides=expand("results/reference/{transcriptome}.fixed_longestORFperGene.fasta", transcriptome=config["reference"]["filename"]),
        peptides=expand("results/reference/{transcriptome}.fixed_longestORFperGene.pep", transcriptome=config["reference"]["filename"]),
        script="workflow/scripts/makeGTF_emapper.py"
    output:
        expand("results/reference/{transcriptome}.fixed_longestORFperGene.fasta.eggnog.gtf", transcriptome=config["reference"]["filename"]),
        geneTranscript_map = expand("results/reference/{transcriptome}.fixed_longestORFperGene.fasta.geneID_to_transcript.txt", transcriptome=config["reference"]["filename"])
    threads: 15
    conda:
        "../envs/emapper.yaml"
    params:
        time="5:00:00",
        mem="50GB",
        geneID_type=config["geneIDType"],
        emapper="$CONDA_PREFIX/bin/emapper.py",
        python="$CONDA_PREFIX/bin/python",
        gonames_file="resources/go_terms_2019.txt",
        data_folder="$CONDA_PREFIX/lib/python3.7/site-packages/data",
        download_eggnog_databases="$CONDA_PREFIX/bin/download_eggnog_data.py"
    shell:
        """
        mkdir {params.data_folder}
        python {params.download_eggnog_databases} -y
        python {input.script} {input.nucleotides} {input.peptides} {params.geneID_type} results/reference {params.python} {params.emapper} {params.gonames_file}
        """
