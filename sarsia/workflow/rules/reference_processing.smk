rule generate_longest_ORFs:
    """
    Infer open reading frames in reference transcriptome.
    """
    input:
        "reference/{transcriptome}"
    output:
        "reference/{transcriptome}.transdecoder_dir/longest_orfs.pep"
    shell:
        "TransDecoder.LongOrfs -t {input} --output_dir reference/{wildcards.transcriptome}.transdecoder_dir"

rule keep_longest_ORF_per_gene:
    """
    Keep longest open reading frames per gene (as per gene identifier). 
    Details on functions and arguments in the python script itself.
    """
    input:
        longestORFs="reference/{transcriptome}.transdecoder_dir/longest_orfs.pep",
        transcriptomeFasta="reference/{transcriptome}",
        script="scripts/keepLongestORFperGene.py"
    output:
        peptides="reference/{transcriptome}_longestORFperGene.pep",
        nucleotides="reference/{transcriptome}_longestORFperGene.fasta"
    shell:
        "python {input.script} -p {input.longestORFs} -t "
        "{input.transcriptomeFasta} -identifier {GENE_ID_TYPE} -o reference"

rule make_GTF:
    input:
        nucleotides="reference/{transcriptome}_longestORFperGene.fasta",
        peptides="reference/{transcriptome}_longestORFperGene.pep",
        script="scripts/makeGTF_emapper.py"
    output:
        "reference/{transcriptome}_longestORFperGene.fasta.eggnog.gtf"
    threads: 15
    params:
        time="70:00:00",
        mem="100GB"
    shell:
        "python {input.script} {input.nucleotides} {input.peptides} {GENE_ID_TYPE} reference"