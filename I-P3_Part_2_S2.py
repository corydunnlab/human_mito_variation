#!/usr/bin/env python
# coding: utf-8

# Funding for this project was received from the European Research Council and the Sigrid Jusélius Foundation.
# Author: Cory Dunn, with some code obtained from (https://github.com/corydunnlab/mammal_mitoprot_evolve, https://www.biorxiv.org/content/10.1101/2021.03.10.434614v1)
# Institution: University of Helsinki
# Author Email: cory.dunn@helsinki.fi
# License: GPLv3
# Version: A

# Load dependencies

import pandas as pd
import numpy as np
import os
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
import fileinput
from Bio import AlignIO
import numpy.polynomial.polynomial as poly

# Features set by user

genbank_file_to_use = 'mito_synonymous_mammals_AND_A_punctatus.gb'
file_prefix = 'S2_' + 'MADPROPS_APR_21_2021_FINAL_'
chosen_genes_set = {'ND1','ND2','COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CYTB'} # vertebrates
chosen_genes = list(chosen_genes_set)
chosen_accession = 'NC_012920.1'
reference_for_ungap = 'NC_012920_1_Homo_sapiens'
besttree = 'mito_synonymous_best_tree_rooted_Anolis_punctatus.nwk'

selected_accession = reference_for_ungap
ancestral_tree = file_prefix + 'coding_DNA_sequence_concatenate_ASR.raxml.ancestralTree'
codon_table = 2 ## vertebrate mitochondrial table
verbose_flag = 'y'

# Ungap coding DNA sequence concatenate

os.system('python ungap_on_reference.py -i ' + file_prefix + 'coding_DNA_sequence_concatenate_MAFFT_FFT_NS_2.fasta -o ' + file_prefix + 'coding_DNA_sequence_concatenate_MAFFT_FFT_NS_2_' + reference_for_ungap + '_ungap.fasta -r ' + reference_for_ungap) # ungap on reference

# Ancestral prediction on ungapped coding DNA sequence concatenate

os.system("raxml-ng --ancestral --msa " + file_prefix + 'coding_DNA_sequence_concatenate_MAFFT_FFT_NS_2_' + reference_for_ungap + '_ungap.fasta --tree ' + besttree + " --model GTR --prefix " + file_prefix + "coding_DNA_sequence_concatenate_ASR") # ancestral prediction by raxml-ng

# Run MAFFT with FFT_NS_2 on the output FASTA CDS files, ungap, then do ancestral prediction using the provided tree

for gene in chosen_genes:
    os.system("mafft --thread 4 --retree 2 --inputorder " + file_prefix + gene + ".fasta > " + file_prefix + gene + "_MAFFT_FFT_NS_2.fasta") # alignment
    os.system('python ungap_on_reference.py -i ' + file_prefix + gene + '_MAFFT_FFT_NS_2.fasta -o ' + file_prefix + gene + '_MAFFT_' + reference_for_ungap + '_ungap.fasta -r ' + reference_for_ungap) # ungap on reference
    os.system("raxml-ng --ancestral --msa " + file_prefix + gene + '_MAFFT_' + reference_for_ungap + '_ungap.fasta --tree ' + besttree + " --model GTR --prefix " + file_prefix + gene + "_ASR") # ancestral prediction by raxml-ng
    for line in fileinput.input(file_prefix + gene + "_ASR.raxml.ancestralStates", inplace=True): # Replace 'Node' with '>Node'
        print (line.replace("Node", ">Node")),
    os.system("seqkit fx2tab " + file_prefix + gene + "_ASR.raxml.ancestralStates > " + file_prefix + gene + "_ASR_ancestralStates.tabular")
    os.system("seqkit fx2tab " + file_prefix + gene + '_MAFFT_' + reference_for_ungap + '_ungap.fasta > ' + file_prefix + gene + '_MAFFT_' + reference_for_ungap + '_ungap.tabular')
    os.system("cat " + file_prefix + gene + "_ASR_ancestralStates.tabular " + file_prefix + gene + '_MAFFT_' + reference_for_ungap + '_ungap.tabular > ' + file_prefix + gene + '_full.tabular')
    os.system("seqkit tab2fx " + file_prefix + gene + '_full.tabular > ' + file_prefix + gene + '_full.fasta')

# Begin building fluctuation tables

# To load the FASTA alignment file

def read_fasta(alignment):  # To read the FASTA alignment file
    aa_dict = {}
    with open(alignment, mode='r') as handle:
        # Using Biopython's parse function to reduce memory footprint
        for record in SeqIO.parse(handle, 'fasta'):
            # Extract individual parts of the FASTA record
            identifier = record.id
            sequence = record.seq
            aa_dict[identifier] = sequence
    return aa_dict

# To retrieve clades by name

def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names

# To retrieve all edges in the tree

def all_edges(tree):

    alledges = []
    for parent in tree.find_clades(terminal=False, order='level'):
        for child in parent.clades:
            alledges.append(parent.name + '*' + child.name)

    return alledges


analyzed_four_fold_positions = []

report_everything_about_selected_positions_results = pd.DataFrame(columns = ['Protein','Degenerate_index','Edge_name','Edge_remark','Ancestral_node','Descendant_node_or_species','Ancestral_character','Descendant_character','Character_remark','Branch_length','L-strand_ancestral','L-strand_descendent'])

for gene in chosen_genes:
    
    # Initialize lists for FASTA information

    record_x_toward_seq_dataframe = []
    sequence_records = []
    alignment_record_name_list = []
    
    # Load input FASTA into dataframe
    
    alignfile = file_prefix + gene + '_full.fasta'
    for record in SeqIO.parse(alignfile,"fasta"):
        alignment_record_name_list.append(record.name)
        record_x_toward_seq_dataframe = list(record.seq)
        record_x_toward_seq_dataframe_UPPER = [x.upper() for x in record_x_toward_seq_dataframe] 
        sequence_records.append(record_x_toward_seq_dataframe_UPPER)

    sequence_dataframe = pd.DataFrame(sequence_records,index=alignment_record_name_list)
    sequence_dataframe_joined = sequence_dataframe.apply(''.join, axis=1)
    reference_sequence_concatenate = sequence_dataframe_joined.at[selected_accession]
    reference_sequence_concatenate_Seq = Seq(reference_sequence_concatenate)
    translated_reference_sequence_concatenate_Seq = reference_sequence_concatenate_Seq.translate(table=codon_table)
    translated_reference_sequence_concatenate_STR = str(translated_reference_sequence_concatenate_Seq)

    # Determine degenerate sites

    indices_by_amino_acid = []
    indices_by_nucleotide = []
    degenerate_count = 0
    for i in range(len(translated_reference_sequence_concatenate_STR)):
        amino_acid = translated_reference_sequence_concatenate_STR[i]
        identity_test_codon_1 = sequence_dataframe[i*3].value_counts(normalize = True)
        identity_test_codon_2 = sequence_dataframe[i*3+1].value_counts(normalize = True)
        if ((identity_test_codon_1 == 1).any()) & ((identity_test_codon_2 == 1).any()):
            if (identity_test_codon_1[identity_test_codon_1 == 1].index[0] == 'A') & (identity_test_codon_2[identity_test_codon_2 == 1].index[0] == 'G') & (amino_acid in 'S'):
                j = i*3+2
                nucleotide_possibilities_at_degenerate_list = sequence_dataframe[j].value_counts(normalize=True).index.tolist()
                nucleotide_possibilities_at_degenerate_string = ''.join(nucleotide_possibilities_at_degenerate_list)
                nucleotide_possibilities_at_degenerate_string = sorted(re.sub('[^ACGT]', '', nucleotide_possibilities_at_degenerate_string))
                if ('A' not in nucleotide_possibilities_at_degenerate_string) & ('G' not in nucleotide_possibilities_at_degenerate_string):
                    degenerate_count += 1
                    indices_by_amino_acid.append(i)
                    indices_by_nucleotide.append(j)
                    degenerate_index = gene + '_' + str(degenerate_count)
                    label_of_site = gene, amino_acid, i+1, (i+1)*3,degenerate_index,nucleotide_possibilities_at_degenerate_string
                    analyzed_four_fold_positions.append(label_of_site)
                    print (label_of_site)

    four_fold_dataframe = sequence_dataframe[indices_by_nucleotide]
    four_fold_dataframe_joined = four_fold_dataframe.apply(''.join, axis=1)

    # Write degenerate sites to FASTA

    ofile = open(file_prefix + 'degenerate_nucleotide_sites_' + gene + '.fasta', "w")
    for seqi in range(len(four_fold_dataframe_joined)):
        ofile.write(">" + alignment_record_name_list[seqi] + "\n" + four_fold_dataframe_joined[seqi] + "\n")
    ofile.close()
    
    alignment_file = file_prefix + 'degenerate_nucleotide_sites_' + gene + '.fasta'
    
    # Read the ancestral tree and record all edges

    from Bio import Phylo
    my_tree = Phylo.read(ancestral_tree, 'newick')
    edges = all_edges(my_tree)
    clades_dict = lookup_by_names(my_tree)

    # Read the four-fold degenerate FASTA file and record character values for each position

    all_values = read_fasta(alignment_file)

    report_everything_about_selected_positions_results_GENE = pd.DataFrame(columns = ['Protein','Alignment_position','Edge_name','Edge_remark','Ancestral_node','Descendant_node_or_species', 'Ancestral_character','Descendant_character','Character_remark','Branch_length','L-strand_ancestral','L-strand_descendent','Degenerate_index'])
    length_of_report_array_count = 0

    for i in range(len(indices_by_nucleotide)):
        query_label = str(i+1)
        query = i
        for edge in edges:
            parent_edge = edge.split('*')[0]
            child_edge = edge.split('*')[1]
            edge_len = str(clades_dict[child_edge].branch_length)
            parent_value = all_values[parent_edge][query]
            child_value = all_values[child_edge][query]
            if edge.count('_') == 0:
                edge_remark = 'Internal'
            elif edge.count('_') > 0:
                edge_remark = 'All'
            if parent_value == child_value:
                value_remark = 'Conserved'
            else:
                value_remark = 'Fluctuating'
            if gene != 'ND6':
                Lanc = parent_value
                Ldes = child_value
            elif gene == 'ND6':
                Lanc = str(Seq(parent_value).reverse_complement())
                Ldes = str(Seq(child_value).reverse_complement())
            other_degenerate_index = gene + '_' + query_label
            to_append_report = ([gene,query_label, edge, edge_remark, parent_edge, child_edge, parent_value, child_value,
                             value_remark, edge_len, Lanc, Ldes,other_degenerate_index])
            report_everything_about_selected_positions_results_GENE.at[length_of_report_array_count]=to_append_report
            length_of_report_array_count += 1
    
    report_everything_about_selected_positions_results_GENE.to_csv(file_prefix + 'fluctuation_analysis_degenerate_nucleotide_sites_' + gene + '.csv')
    report_everything_about_selected_positions_results
    report_everything_about_selected_positions_results = pd.concat([report_everything_about_selected_positions_results,report_everything_about_selected_positions_results_GENE], sort=False)

analyzed_four_fold_positions_DF = pd.DataFrame(analyzed_four_fold_positions, columns = ['Protein','Amino_acid','Amino_acid_position','CDS_nucleotide_position','Degenerate_index','Nucleotide_possibilities_at_degenerate (ACGT)'])

report_everything_after_merge = pd.merge(report_everything_about_selected_positions_results, analyzed_four_fold_positions_DF, on = ['Protein','Degenerate_index'], how = "left")
report_everything_after_merge.to_csv(file_prefix + 'fluctuation_analysis_after_merge.csv')

# Integrate TSS values

substitutions = report_everything_after_merge[report_everything_after_merge['L-strand_ancestral'] != report_everything_after_merge['L-strand_descendent']]
substitutions = substitutions[(substitutions['L-strand_ancestral'].isin(['A','C','G','T'])) & (substitutions['L-strand_descendent'].isin(['A','C','G','T']))] 
substitutions['Alignment_position'] = substitutions['Alignment_position'].astype(str)
substitutions['Degenerate_index'] = substitutions['Protein'] + '_' + substitutions['Alignment_position']

TSS_DF = substitutions['Degenerate_index'].value_counts()
TSS_DF = TSS_DF.reset_index()
TSS_DF.rename(columns = {'Degenerate_index':'TSS'}, inplace = True)
TSS_DF.rename(columns = {'index':'Degenerate_index'}, inplace = True)
TSS_DF.to_csv(file_prefix + 'TSS.csv')

FINAL_A = substitutions.set_index('Degenerate_index').join(TSS_DF.set_index('Degenerate_index'))
FINAL_A = FINAL_A.sort_values(by='TSS', ascending=True)

FINAL_A.to_csv(file_prefix + 'fluctuation_analysis_with_TSS_no_odd_characters_or_gaps.csv')

FINAL_B = analyzed_four_fold_positions_DF.set_index('Degenerate_index').join(TSS_DF.set_index('Degenerate_index'))
FINAL_B.to_csv(file_prefix + 'analyzed_degenerate_positions_all_proteins_w_TSS.csv')
