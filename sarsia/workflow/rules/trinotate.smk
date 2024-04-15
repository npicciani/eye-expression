rule trinotate_dbs:
    """
    Download databases for Trinotate v3.2.2
    """
    output:
        "results/trinotate/Pfam-A.hmm",
        "results/trinotate/uniprot_sprot.pep"
    conda:
        "../envs/trinotate.yaml"
    params:
        sqlite = "$CONDA_PREFIX/bin/Build_Trinotate_Boilerplate_SQLite_db.pl"
    shell:
        """
        # Build databases
        mkdir results/trinotate
        perl {params.sqlite} results/trinotate/Trinotate
        mv uniprot_sprot.* results/trinotate
        mv Pfam-A.hmm.gz results/trinotate

        cd results/trinotate
        makeblastdb -in uniprot_sprot.pep -dbtype prot
        gunzip Pfam-A.hmm.gz
        hmmpress Pfam-A.hmm
        """
rule trinotate_hits:
    """
    Run blastx, hmmer, signalp and rnammer then load results on sql database
    """
    input:
        "results/trinotate/uniprot_sprot.pep",
        "results/trinotate/Pfam-A.hmm",
        fastafile = expand("results/reference/{transcriptome}_longestORFperGene.fasta", transcriptome=config["reference"]["filestem"]),
        peptidefile = expand("results/reference/{transcriptome}_longestORFperGene.pep", transcriptome=config["reference"]["filestem"]),
        transdecoder_peptidefile = expand("results/reference/{transcriptome}.transdecoder_dir/longest_orfs.pep", transcriptome=config["reference"]["filestem"]),
        geneTranscript_map = expand("results/reference/{transcriptome}_longestORFperGene.fasta.geneID_to_transcript.txt", transcriptome=config["reference"]["filestem"]),
    output:
        "results/trinotate/blastx.outfmt6",
        "results/trinotate/blastp.outfmt6",
        "results/trinotate/pfam.log",
        "results/trinotate/signalp.out",
        "results/trinotate/tmhmm.out",
        "results/trinotate/trinotate_annotation_report.xls"
    threads: 12
    conda:
        "../envs/trinotate.yaml"
    params:
        signalp = expand("{signalpPath}", signalpPath=config["signalp"]),
        tmhmm = expand("{tmhmmPath}", tmhmmPath=config["tmhmm"]),
        rnammer = expand("{rnammerPath}", rnammerPath=config["rnammer"]),
        rnammerUtil = "$CONDA_PREFIX/bin/RnammerTranscriptome.pl"
    shell:
        """
        cd results/trinotate
        # Blast against uniprot
        blastx -query ../../{input.fastafile} -db uniprot_sprot.pep -num_threads 12 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx.outfmt6
        blastp -query ../../{input.peptidefile} -db uniprot_sprot.pep -num_threads 12 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp.outfmt6

        # Run HMMER to identify protein domains
        hmmscan --cpu {threads} --domtblout TrinotatePFAM.out Pfam-A.hmm ../../{input.peptidefile} > pfam.log

        # Run signalP to predict signal peptides
        {params.signalp} -f short -n signalp.out ../../{input.peptidefile}

        # Run tmHMM to predict transmembrane regions
        {params.tmhmm} --short < ../../{input.peptidefile} > tmhmm.out

        # Run RNAMMER to identify rRNA genes
        perl {params.rnammerUtil} --transcriptome ../../{input.fastafile} --path_to_rnammer {params.rnammer}

        # Load results into SQL database
        Trinotate Trinotate.sqlite init --gene_trans_map ../../{input.geneTranscript_map} \
                                --transcript_fasta ../../{input.fastafile} \
                                --transdecoder_pep ../../{input.transdecoder_peptidefile}
        Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6 #load protein hits
        Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6 #load transcript hits
        Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
        Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out
        Trinotate Trinotate.sqlite LOAD_signalp signalp.out
        Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls
        """

rule trinotate_extract:
    """
    Extract GO terms assigned to each transcript
    """
    input:
        "results/trinotate/trinotate_annotation_report.xls"
    output:
        "results/trinotate/go_annotations.txt"
    conda:
        "../envs/trinotate.yaml"
    params:
        extract = "$CONDA_PREFIX/bin/extract_GO_assignments_from_Trinotate_xls.pl",
        conda_folder = "$CONDA_PREFIX/lib/site_perl/5.26.2/x86_64-linux-thread-multi"
    shell:
        """
        mkdir {params.conda_folder}/obo
        wget https://github.com/Trinotate/Trinotate/blob/master/PerlLib/obo/go-basic.obo.gz -P {params.conda_folder}/obo #fix conda installation bug

        cd results/trinotate
        # Extract GO annotations
        perl {params.extract} \
                --Trinotate_xls trinotate_annotation_report.xls \
                -G --include_ancestral_terms \
                > go_annotations.txt
        """
