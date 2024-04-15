rule generate_longest_ORFs:
    """
    Infer open reading frames in reference transcriptome.
    """
    input:
        transcriptomePath=expand("{transcriptome_path}", transcriptome_path=config["reference"]["path"]),
    output:
        expand("results/reference/{transcriptome}.transdecoder_dir/longest_orfs.pep", transcriptome=config["reference"]["filestem"])
    params:
        referenceFilename=expand("{transcriptome}", transcriptome=config["reference"]["filestem"])
    conda:
        "../../workflow/envs/transdecoder.yaml" #transdecoder v5.5.0
    shell:
        "TransDecoder.LongOrfs -t {input.transcriptomePath} --output_dir results/reference/{params.referenceFilename}.transdecoder_dir"

rule keep_longest_ORF_per_gene:
    """
    Keep longest open reading frames per gene (as per gene identifier type defined in the config yaml).
    Details on functions and arguments in the python script itself.
    """
    input:
        longestORFs=expand("results/reference/{transcriptome}.transdecoder_dir/longest_orfs.pep", transcriptome=config["reference"]["filestem"]),
        transcriptomePath=expand("{transcriptome_path}", transcriptome_path=config["reference"]["path"]),
        script="workflow/scripts/keepLongestORFperGene.py"
    output:
        peptides=expand("results/reference/{transcriptome}_longestORFperGene.pep", transcriptome=config["reference"]["filestem"]),
        nucleotides=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["filestem"])
    params:
        geneID_type=config["geneIDType"]
    shell:
        "python {input.script} -p {input.longestORFs} -t "
        "{input.transcriptomePath} -identifier {params.geneID_type} -o results/reference"

rule make_GTF:
    input:
        nucleotides=expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["filestem"]),
        peptides=expand("results/reference/{transcriptome}_longestORFperGene.pep", transcriptome=config["reference"]["filestem"]),
        script="workflow/scripts/makeGTF_emapper.py"
    output:
        expand("results/reference/{transcriptome}_longestORFperGene.fasta.eggnog.gtf", transcriptome=config["reference"]["filestem"]),
        expand("results/reference/{transcriptome}_longestORFperGene.fasta.geneID_to_transcript.txt", transcriptome=config["reference"]["filestem"])
    threads: 15
    conda:
    	"../../workflow/envs/emapper.yaml"
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