#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created by Natasha Picciani on Oct-14-2019 using Python 3.7
# Version 2
# Contact natasha.picciani@gmail.com with questions/comments
# Last modified: Jun-7-2021
# Changes in v2: added nucleotide output file matching selected transcripts that generate longest ORF per gene
# Usual Command example: python keepLongestORFperGene.py -p longest_orfs.pep -t transcriptome.fasta -o ./test -identifier type_1
# Adapted for snakemake workflow

"""
Read a FASTA file with open reading frames (ORFs), and keep the longest
ORF for each gene. Read the original nucleotide fasta file used to generate ORFs
and select the nucleotide sequence corresponding to the longest ORFs.
Return two fasta files: one containing the aminoacid sequences of longest ORFs
per gene (.pep) and another containing their corresponding nucleotides (.fasta).
"""

import re
import argparse
from pathlib import PurePosixPath
from Bio import SeqIO


def getID(header, type, whichID):

    """
    Search the header of a fasta file and return the gene identifier of a sequence."

    Arguments:
    header -- the header string (the record.id value of a sequence when using SeqIO package).
    type -- the type of gene identifier; type_1 = "compXX_cXX_seqXX", type_2="TRINITY_DNXXXXX_cX_gX_iX",
                    type_3="SegXX.XX.XXXX".
    whichID -- the kind of identifier to be retrieved; geneID = "gene", transcriptID = "transcript"
    """

    if type == "type_1":
        searchStr = "((comp\d+_c\d+)_seq\d+)"
    if type == "type_2":
        searchStr = "((TRINITY_DN\d+_c\d+_g\d+)_i\d+)"
    if type == "type_3":
        searchStr = "((Seg\d+)\..+)"
    if type == "type_4":
        searchStr = "(t\d+aep.+(DN.+g.+)[_it]\d+)"
    parts = header.split("\t")
    p = re.match(searchStr, parts[0])
    geneID = p.group(2)
    transcriptID = p.group(1)

    if whichID == "gene":
        return geneID
    if whichID == "transcript":
        return transcriptID


def keepLongest(ORFs, transcriptomeFile, identifier_type, output_directory):

    """
    Read a fasta file with ORFs (open reading frames), a fasta file with original transcriptome,
    select the longest ORF per each gene from ORFs file and select their nucleotide sequence
    from transcriptome. Return two FASTA files: one with aminoacid sequences and
    another with nucleotides.

    Arguments:
    ORFs -- FASTA file with ORF sequences
    transcriptomeFile -- FASTA file with nucleotide sequences
    identifier_type -- the type of gene identifier; type_1 = "compXX_cXX_seqXX", type_2="TRINITY_DNXXXXX_cX_gX_iX",
                    type_3="SegXX.XX.XXXX".
    output_directory -- path to directory for placing output files
    """

    seqs = {}
    currentGeneID = ""
    longestORFperGenePep = (
        output_directory
        + "/"
        + PurePosixPath(transcriptomeFile).name
        + "_longestORFperGene.pep"
    )  # uncomment for normal script
    longestORFperGeneFasta = (
        output_directory
        + "/"
        + PurePosixPath(transcriptomeFile).name
        + "_longestORFperGene.fasta"
    )  # uncomment for normal script
    # 	longestORFperGenePep = snakemake.output[0] #snakemake output file included in rule "keep_longest_ORF_per_gene"; comment for normal script
    # 	longestORFperGeneFasta = snakemake.output[1] #snakemake output file included in rule "keep_longest_ORF_per_gene"; comment for normal script
    # 	unique = 0 #uncomment for testing
    # 	duplicates = 0 #uncomment for testing

    with open(longestORFperGenePep, "w") as outfile1:

        for record in SeqIO.parse(ORFs, "fasta"):  # find gene names
            geneID = getID(record.id, identifier_type, whichID="gene")
            if currentGeneID == "":  # start the loop with an empty GeneID string
                currentGeneID = geneID
                seqs[record.description] = [
                    geneID,
                    record.seq,
                ]  # add first sequence to dictionary "seqs", which has geneID as a list value
            else:
                if (
                    geneID != currentGeneID
                ):  # if the previous geneID is different from that of the current record
                    # 					unique += 1 #uncomment for testing
                    seqs[record.description] = [
                        geneID,
                        record.seq,
                    ]  # add unique sequences to dictionary "seqs"
                    currentGeneID = geneID  # re-set current geneID
                else:  # if the geneID is the same as that of the previous record
                    # 					duplicates += 1 #uncomment for testing
                    for (
                        key,
                        value,
                    ) in (
                        seqs.items()
                    ):  # search for the previous record in the dictionary "seqs"
                        if (
                            value[0] == currentGeneID
                        ):  # and save its sequence length and name
                            dupeLength = len(value[1])
                            dupeID = key
                    if dupeLength < len(
                        record.seq
                    ):  # if length of previous record is smaller that of than current record
                        del seqs[dupeID]  # get rid of previous record
                        seqs[record.description] = [
                            geneID,
                            record.seq,
                        ]  # add current record to dictionary
                        currentGeneID = geneID

        for name, seq in seqs.items():  # print sequences to an output file
            name = name.split(" ")
            name = name[0]
            print(">" + name, file=outfile1)
            print(seq[1], file=outfile1)

    with open(longestORFperGeneFasta, "w") as outfile2:
        querylist = []
        for record in SeqIO.parse(longestORFperGenePep, "fasta"):
            transcript = getID(record.id, identifier_type, whichID="transcript")
            querylist.append(transcript)

        input_seq_iterator = SeqIO.parse(transcriptomeFile, "fasta")
        short_seq_iterator = (
            record
            for record in input_seq_iterator
            for item in querylist
            if item.strip("\n") == record.id
        )
        SeqIO.write(short_seq_iterator, outfile2, "fasta")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Read FASTA files with open reading frames(ORFs) and nucleotides, and keep the longest ORF per gene"
    )

    parser.add_argument(
        "-p",
        metavar="FASTA_file",
        type=str,
        required=True,
        help="FASTA file with ORFs, peptides",
    )
    parser.add_argument(
        "-t",
        metavar="FASTA_file",
        type=str,
        required=True,
        help="FASTA file with nucleotides, transcriptome",
    )
    parser.add_argument(
        "-o",
        metavar="output_directory",
        type=str,
        required=True,
        default=".",
        help="directory for placing output file",
    )
    parser.add_argument(
        "-identifier",
        metavar="identifier_type",
        type=str,
        choices=["type_1", "type_2", "type_3", "type_4"],
        required=False,
        default="type_1",
        help="type of gene identifier in the ORF FASTA file",
    )
    parser.add_argument("-version", action="version", version="1.0")

    args = parser.parse_args()

    # Uncomment this following lines for normal script
    outdir = args.o
    ORFs = args.p
    transcriptomeFile = args.t
    idType = args.identifier

    # Snakemake input arguments; comment for normal script
    # ORFs = snakemake.input[0]
    # transcriptomeFile = snakemake.input[1]
    # idType = snakemake.input[2]
    # outdir = snakemake.input[3]

    # Select longest ORFs per gene
    keepLongest(ORFs, transcriptomeFile, idType, output_directory=outdir)
