#!/usr/bin/env python
# coding: utf-8 

# Funding for this project was received from the European Research Council and the Sigrid Jusélius Foundation.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.dunn@helsinki.fi
# License: GPLv3
# Version: A
# Please cite DOI: 10.1101/2021.05.12.443844

# Load dependencies

import pandas as pd
import numpy as np
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo
from Bio import Entrez

genbank_file_to_use = 'NC_012920.gbk'
prefix = 'NC_012920_'
start_position_in_genome = 1
possible_nucleotides = 'ACGT'

codon_table = 2 # 2 is vertebrate mitochondrial

# Determine protein content and make lists

CDS_list = []
accession_error_list = []
accession_name_list = []
for seq_record in SeqIO.parse(genbank_file_to_use, "genbank"):
    accession = str(seq_record.id)
    accession_name_list.append(accession)
    for feat in range(len(seq_record.features)):
        feature_type = seq_record.features[feat].type
        if feature_type == 'CDS':
            try:
                prot_id = str(seq_record.features[feat].qualifiers['gene'][0])
                CDS_list.append(prot_id)
            except KeyError:
                accession_error_list.append(accession)

print('Protein set detected: ',CDS_list)

#Set up gene dataframe

gene_array = pd.DataFrame(index=CDS_list,columns=['Coding_sequence','Start_nucleotide','End_nucleotide','Strand'])
gene_array.index.name = 'Protein_coding_sequence'

# Move gene sequences to gene dataframe

accession_name_list = []
for seq_record in SeqIO.parse(genbank_file_to_use, "genbank"):
    accession = str(seq_record.id)
    gb_entire_sequence = list(seq_record.seq)
    gb_entire_sequence_UPPER = [x.upper() for x in gb_entire_sequence]
    gb_entire_sequence_joined = ''.join(gb_entire_sequence_UPPER)
    for feat in range(len(seq_record.features)):
        feature_type = seq_record.features[feat].type
        if feature_type == 'CDS':
            try:
                feature_start_zero_based_numbering = seq_record.features[feat].location.nofuzzy_start
                feature_end_zero_based_numbering = seq_record.features[feat].location.nofuzzy_end
                feature_strand = seq_record.features[feat].strand
                gene_id = str(seq_record.features[feat].qualifiers['gene'][0])
                sequence_slice = gb_entire_sequence_joined[feature_start_zero_based_numbering:feature_end_zero_based_numbering]
                if feature_strand == 1:
                    gene_array.at[gene_id,'Coding_sequence'] = sequence_slice
                if feature_strand == -1:
                    sequence_slice_BP = Seq(sequence_slice)
                    sequence_slice_rvscomp = str(sequence_slice_BP.reverse_complement())
                    gene_array.at[gene_id,'Coding_sequence'] = sequence_slice_rvscomp
                gene_array.at[gene_id,'Start_nucleotide'] = feature_start_zero_based_numbering + 1
                gene_array.at[gene_id,'End_nucleotide'] = feature_end_zero_based_numbering
                gene_array.at[gene_id,'Strand'] = feature_strand
                
            except KeyError:
                pass

def codon_change(input_sequence,test_protein,start_position_in_genome,end_position_in_genome,coding_strand):
    
    possible_codon_mutant_list = []
    number_of_codons = len(input_sequence) // 3
    for i in range(number_of_codons):
        
        start_position = ((i+1)*3 - 3)
        stop_position = ((i+1)*3)
        selected_codon = input_sequence[start_position:stop_position]
        
        # Codon position 1
  
        for nuc in possible_nucleotides:
            WTposition1 = selected_codon[0]
            positions2and3 = selected_codon[1:3]
            mutant_nuc = nuc+positions2and3
            WT_nuc_seq = Seq(selected_codon)
            mutant_nuc_seq = Seq(mutant_nuc)
            WT_aa_seq = WT_nuc_seq.translate(table=codon_table)
            mutant_aa_seq = mutant_nuc_seq.translate(table=codon_table)
            if WT_aa_seq == '*':
                WT_aa_seq = 'X'
            if mutant_aa_seq == '*':
                mutant_aa_seq = 'X'
            if coding_strand == 1:
                nucleotide_position1 = str(start_position_in_genome - 1 + i*3+1)
                nucleotide_change1 = WTposition1 + nucleotide_position1 + nuc
                if WTposition1 != nuc:
                    info_to_append = nucleotide_change1,test_protein+'_'+str(WT_aa_seq)+str(i+1)+str(mutant_aa_seq)
                    possible_codon_mutant_list.append(info_to_append)
            else:
                nucleotide_position1 = str(end_position_in_genome - i*3)
                nucleotide_change1 = str(Seq(WTposition1).complement()) + nucleotide_position1 + str(Seq(nuc).complement()) ###
                if str(Seq(WTposition1).complement()) != str(Seq(nuc).complement()):
                    info_to_append = nucleotide_change1,test_protein+'_'+str(WT_aa_seq)+str(i+1)+str(mutant_aa_seq)
                    possible_codon_mutant_list.append(info_to_append)
            
        # Codon position 2
        
        for nuc in possible_nucleotides:
            position1 = selected_codon[0]
            WTposition2 = selected_codon[1]
            position3 = selected_codon[2]
            mutant_nuc = position1+nuc+position3
            WT_nuc_seq = Seq(selected_codon)
            mutant_nuc_seq = Seq(mutant_nuc)
            WT_aa_seq = WT_nuc_seq.translate(table=codon_table)
            mutant_aa_seq = mutant_nuc_seq.translate(table=codon_table)
            if WT_aa_seq == '*':
                WT_aa_seq = 'X'
            if mutant_aa_seq == '*':
                mutant_aa_seq = 'X'
            if coding_strand == 1:
                nucleotide_position2 = str(start_position_in_genome - 1 + i*3+2)
                nucleotide_change2 = WTposition2 + nucleotide_position2 + nuc
                if WTposition2 != nuc:
                    info_to_append = nucleotide_change2,test_protein+'_'+str(WT_aa_seq)+str(i+1)+str(mutant_aa_seq)
                    possible_codon_mutant_list.append(info_to_append)
            else:
                nucleotide_position2 = str(end_position_in_genome - i*3 - 1)
                nucleotide_change2 = str(Seq(WTposition2).complement()) + nucleotide_position2 + str(Seq(nuc).complement()) ###
                if str(Seq(WTposition2).complement()) != str(Seq(nuc).complement()):
                    info_to_append = nucleotide_change2,test_protein+'_'+str(WT_aa_seq)+str(i+1)+str(mutant_aa_seq)
                    possible_codon_mutant_list.append(info_to_append)
            
        # Codon position 3
        
        for nuc in possible_nucleotides:
            positions1and2 = selected_codon[0:2]
            WTposition3 = selected_codon[2]
            mutant_nuc = positions1and2+nuc
            WT_nuc_seq = Seq(selected_codon)
            mutant_nuc_seq = Seq(mutant_nuc)
            WT_aa_seq = WT_nuc_seq.translate(table=codon_table)
            mutant_aa_seq = mutant_nuc_seq.translate(table=codon_table)
            if WT_aa_seq == '*':
                WT_aa_seq = 'X'
            if mutant_aa_seq == '*':
                mutant_aa_seq = 'X'
            if coding_strand == 1:
                nucleotide_position3 = str(start_position_in_genome - 1 + i*3+3)
                nucleotide_change3 = WTposition3 + nucleotide_position3 + nuc
                if WTposition3 != nuc:
                    info_to_append = nucleotide_change3,test_protein+'_'+str(WT_aa_seq)+str(i+1)+str(mutant_aa_seq)
                    possible_codon_mutant_list.append(info_to_append)
            else:
                nucleotide_position3 = str(end_position_in_genome - i*3 - 2)
                nucleotide_change3 = str(Seq(WTposition3).complement()) + nucleotide_position3 + str(Seq(nuc).complement()) ###
                if str(Seq(WTposition3).complement()) != str(Seq(nuc).complement()):
                    info_to_append = nucleotide_change3,test_protein+'_'+str(WT_aa_seq)+str(i+1)+str(mutant_aa_seq)
                    possible_codon_mutant_list.append(info_to_append)
        
    return(possible_codon_mutant_list)

# Define final dataframe

final_dataframe = pd.DataFrame(columns = ['Nucleotide_change','Amino_acid_change','Protein'])

# Run through sequences and produce results

for i, protein in gene_array.iterrows():
    coding_sequence_to_send = gene_array.loc[i,'Coding_sequence']
    start_position_of_coding = gene_array.loc[i,'Start_nucleotide']
    end_position_of_coding = gene_array.loc[i,'End_nucleotide']
    strnd = gene_array.loc[i,'Strand']
    list_of_codon_change_output = codon_change(coding_sequence_to_send,i,start_position_of_coding,end_position_of_coding,strnd)
    globals()[i] = pd.DataFrame(list_of_codon_change_output,columns=['Nucleotide_change','Amino_acid_change'])
    globals()[i]['Protein'] = i
    globals()[i].to_csv(prefix + i + '.csv',index=False)
    final_dataframe = pd.concat([final_dataframe,globals()[i]], axis=0)

# Final dataframe

final_dataframe.to_csv(prefix + 'final_dataframe_all_single_nucleotide_substitutions_as_protein_changes.csv',index=False)
print(final_dataframe)
