# -*- coding: utf-8 -*-
# Written by Natasha Picciani
# Feb 22, 2023

import sys
import ahocorasick

annotations_master = sys.argv[1]
gene_trees_master = sys.argv[2]
outdir = sys.argv[3]
gene_trees_appended = f"{outdir}/gene_trees.master.annotated.txt"

def make_GO_dictionary(emmaper_annotation_file):
    """Make dictionary with Protein IDS and corresponding GOs pulled
    from functional annotation file generated with Emapper 2.1.10
    """
    GO_dict={}
    with open(emmaper_annotation_file, "r") as infile:
        for line in infile:
            line = line.strip("\n")
            columns = line.split("\t")
            proteinID = columns[0]
            goterms_field = columns[9]  # GOs
            GO_dict[proteinID] = goterms_field
    return GO_dict

def make_gene_names_dictionary(emmaper_annotation_file):
    """Make dictionary with Protein IDS and corresponding Gene Names pulled
    from functional annotation file generated with Emapper 2.1.10
    """
    gene_names_dict={}
    with open(emmaper_annotation_file, "r") as infile:
        for line in infile:
            line = line.strip("\n")
            columns = line.split("\t")
            proteinID = columns[0]
            gene_name = columns[8]
            proteinID_gene_name = f"{proteinID}_{gene_name}"
            gene_names_dict[proteinID] = proteinID_gene_name
    return gene_names_dict

def make_automaton(input_dictionary):
    """Build an Aho-Corasick automaton from a dictionary file and return
    it.
    Source: https://codereview.stackexchange.com/questions/170987/optimize-search-and-replace-in-one-file-based-on-dictionary-in-another-file
    """
    automaton = ahocorasick.Automaton()
    for key,value in input_dictionary.items():
        automaton.add_word(key, (key, value))
    automaton.make_automaton()
    return automaton

def apply_automaton(automaton, input_filename, output_filename):
    """Apply an Aho-Corasick automaton to an input file, replacing the
    first occurrence of a key in each line with the corresponding
    value, and writing the result to the output file.
    Source: https://codereview.stackexchange.com/questions/170987/optimize-search-and-replace-in-one-file-based-on-dictionary-in-another-file
    """
    with open(input_filename) as infile, open(output_filename, 'w') as outfile:
        for line in infile:
            for end, (key, value) in automaton.iter(line):
                line= line.replace(key, value)
            outfile.write(line)


gene_names_dict=make_gene_names_dictionary(annotations_master)
gene_names_automaton=make_automaton(gene_names_dict)
apply_automaton(gene_names_automaton,gene_trees_master,gene_trees_appended)