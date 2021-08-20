#!/bin/bash

# #load environment trinotate2
# #programs
# blastx=/home/nnp9/project/conda_envs/trinotate2/bin/blastx
# blastp=/home/nnp9/project/conda_envs/trinotate2/bin/blastp
# hmmscan=/home/nnp9/project/conda_envs/trinotate2/bin/hmmscan
# signalp=/home/nnp9/project/conda_envs/trinotate2/bin/signalp-4.1/signalp
# tmhmm=/home/nnp9/project/conda_envs/trinotate2/bin/tmhmm-2.0c/bin/tmhmm
# rnammerUtil=/home/nnp9/project/conda_envs/trinotate2/bin/Trinotate-Trinotate-v3.2.2/util/rnammer_support/RnammerTranscriptome.pl
# rnammer=/home/nnp9/project/conda_envs/trinotate2/bin/rnammer-1.2/rnammer
# trinotate=/home/nnp9/project/conda_envs/trinotate2/bin/Trinotate-Trinotate-v3.2.2/Trinotate
# extract=/home/nnp9/project/conda_envs/trinotate2/bin/Trinotate-Trinotate-v3.2.2/util/extract_GO_assignments_from_Trinotate_xls.pl

## files and folders
fastafile = $1
peptidefile = $2
transdecoder_peptidefile = $3
geneTranscript_map = $4
outdir = $5

## databases
uniprotdb=/gpfs/ysm/project/dunn/nnp9/conda_envs/trinotate2/bin/Trinotate-Trinotate-v3.2.2/uniprot_sprot.pep
pfamdb=/gpfs/ysm/project/dunn/nnp9/conda_envs/trinotate2/bin/Trinotate-Trinotate-v3.2.2/Pfam-A.hmm

$blastx -query $fastafile -db $uniprotdb -num_threads 12 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > $outdir/blastx.outfmt6
$blastp -query $peptidefile -db $uniprotdb -num_threads 12 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > $outdir/blastp.outfmt6

Running HMMER to identify protein domains
$hmmscan --cpu 12 --domtblout $outdir/TrinotatePFAM.out $pfamdb $peptidefile > $outdir/pfam.log

Running signalP to predict signal peptides
$signalp -f short -n $outdir/signalp.out $peptidefile

Running tmHMM to predict transmembrane regions
$tmhmm --short < $peptidefile > $outdir/tmhmm.out

cd $outdir

#perl $rnammerUtil --transcriptome $fastafile --path_to_rnammer $rnammer

#concatenating all the information
$trinotate Trinotate.sqlite init --gene_trans_map $geneTranscript_map \
--transcript_fasta $fastafile --transdecoder_pep $transdecoder_peptidefile

$trinotate Trinotate.sqlite LOAD_swissprot_blastp $outdir/blastp.outfmt6 #load protein hits
$trinotate Trinotate.sqlite LOAD_swissprot_blastx $outdir/blastx.outfmt6 #load transcript hits
$trinotate Trinotate.sqlite LOAD_pfam $outdir/TrinotatePFAM.out
$trinotate Trinotate.sqlite LOAD_tmhmm $outdir/tmhmm.out
$trinotate Trinotate.sqlite LOAD_signalp $outdir/signalp.out

$trinotate Trinotate.sqlite report > $outdir/trinotate_annotation_report.xls

cd $outdir

perl $extract --Trinotate_xls $outdir/trinotate_annotation_report.xls \
                         -G --include_ancestral_terms \
                         > go_annotations.txt
