def get_geneID_type(wildcards):
    return config["geneIDType"][wildcards.geneID_type]

rule generate_longest_ORFs:
    """
    Infer open reading frames in reference transcriptome.
    """
    input:
        expand("results/reference/{transcriptome}", transcriptome=config["reference"])
    output:
        expand("results/reference/{transcriptome}.transdecoder_dir/longest_orfs.pep", transcriptome=config["reference"])
    shell:
        "TransDecoder.LongOrfs -t {input} --output_dir results/reference/{wildcards.transcriptome}.transdecoder_dir"

rule keep_longest_ORF_per_gene:
    """
    Keep longest open reading frames per gene (as per gene identifier type defined in the config yaml). 
    Details on functions and arguments in the python script itself.
    """
    input:
        longestORFs=expand("results/reference/{transcriptome}.transdecoder_dir/longest_orfs.pep", transcriptome=config["reference"]),
        transcriptomeFasta=expand("results/reference/{transcriptome}", transcriptome=config["reference"]),
        script="scripts/keepLongestORFperGene.py"
    output:
        peptides=expand("results/reference/{transcriptome}_longestORFperGene.pep", transcriptome=config["reference"]),
        nucleotides=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"])
    shell:
        "python {input.script} -p {input.longestORFs} -t "
        "{input.transcriptomeFasta} -identifier {geneID_type} -o reference"

rule make_GTF:
    input:
        nucleotides=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]),
        peptides=expand("results/reference/{transcriptome}_longestORFperGene.pep", transcriptome=config["reference"]),
        script="scripts/makeGTF_emapper.py"
    output:
        expand("results/reference/{transcriptome}_longestORFperGene.fasta.eggnog.gtf", transcriptome=config["reference"])
    threads: 15
    params:
        time="70:00:00",
        mem="100GB"
    shell:
        "python {input.script} {input.nucleotides} {input.peptides} {gene_ID_type} reference"