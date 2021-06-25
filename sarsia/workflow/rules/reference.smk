rule generate_longest_ORFs:
    """
    Infer open reading frames in reference transcriptome.
    """
    input:
        transcriptomePath=expand("resources/{transcriptome}", transcriptome=config["reference"]),
    output:
        expand("results/reference/{transcriptome}.transdecoder_dir/longest_orfs.pep", transcriptome=config["reference"])
    params:
        referenceStem=expand("{transcriptome}", transcriptome=config["reference"])
    shell:
        "TransDecoder.LongOrfs -t {input.transcriptomePath} --output_dir results/reference/{params.referenceStem}.transdecoder_dir"

rule keep_longest_ORF_per_gene:
    """
    Keep longest open reading frames per gene (as per gene identifier type defined in the config yaml). 
    Details on functions and arguments in the python script itself.
    """
    input:
        longestORFs=expand("results/reference/{transcriptome}.transdecoder_dir/longest_orfs.pep", transcriptome=config["reference"]),
        transcriptomePath=expand("resources/{transcriptome}", transcriptome=config["reference"]),
        script="workflow/scripts/keepLongestORFperGene.py"
    output:
        peptides=expand("results/reference/{transcriptome}_longestORFperGene.pep", transcriptome=config["reference"]),
        nucleotides=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"])
    params:
        geneID_type=config["geneIDType"]
    shell:
        "python {input.script} -p {input.longestORFs} -t "
        "{input.transcriptomePath} -identifier {params.geneID_type} -o results/reference"

rule make_GTF:
    input:
        nucleotides=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]),
        peptides=expand("results/reference/{transcriptome}_longestORFperGene.pep", transcriptome=config["reference"]),
        script="workflow/scripts/makeGTF_emapper.py"
    output:
        expand("results/reference/{transcriptome}_longestORFperGene.fasta.eggnog.gtf", transcriptome=config["reference"])
    threads: 15
    params:
        time="70:00:00",
        mem="100GB",
        geneID_type=config["geneIDType"]
    shell:
        "python {input.script} {input.nucleotides} {input.peptides} {params.geneID_type} results/reference"