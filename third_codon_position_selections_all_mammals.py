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
from Bio import SeqIO
from Bio.Seq import Seq

# Features set by user

genbank_file_to_use = 'mito_synonymous_mammals_WO_A_punctatus.gb'
file_prefix = 'codon_third_positions_APR_27_2021_'
chosen_genes_set = {'ND1','ND2','COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CYTB'} # vertebrates
chosen_genes = list(chosen_genes_set)
codon_table = 2 ## vertebrate mitochondrial table

# Determine protein content and make lists

chosen_genes = list(chosen_genes_set)

CDS_list = []
gene_list = []
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
            
        if feature_type == 'gene':
            try:
                gene_id = str(seq_record.features[feat].qualifiers['gene'][0])
                gene_list.append(gene_id)
            except KeyError:
                accession_error_list.append(accession)
CDS_list_rem_dups = sorted(list(set(CDS_list)))
gene_list_rem_dups = sorted(list(set(gene_list)))
accession_error_list_rem_dups = sorted(list(set(accession_error_list)))
print('Protein set detected: ',CDS_list_rem_dups)

# Set up protein dataframe

protein_array = pd.DataFrame(index=accession_name_list,columns=CDS_list_rem_dups)
protein_array.index.name = 'Accession'

# Move translation products to protein dataframe

accession_name_list = []
for seq_record in SeqIO.parse(genbank_file_to_use, "genbank"):
    accession = str(seq_record.id)
    accession_name_list.append(accession)
    for feat in range(len(seq_record.features)):
        feature_type = seq_record.features[feat].type
        if feature_type == 'CDS':
            try:
                feature_translation = str(seq_record.features[feat].qualifiers['translation'][0])
                gene_id = str(seq_record.features[feat].qualifiers['gene'][0])
                protein_array.at[accession,gene_id] = feature_translation
            except KeyError:
                pass
            

# Count and record instances of each protein in the dataframe

protein_instances_list = []
for prot in CDS_list_rem_dups:
    protein_to_append = prot, len(protein_array)-protein_array[prot].isnull().sum()
    protein_instances_list.append(protein_to_append) 

protein_instances_list_DF = pd.DataFrame(protein_instances_list,columns = ['Protein','Instances'])

#Set up gene dataframe

gene_array = pd.DataFrame(index=accession_name_list,columns=CDS_list_rem_dups)
gene_array.index.name = 'Accession'

# Move gene sequences to gene dataframe

accession_name_list = []
for seq_record in SeqIO.parse(genbank_file_to_use, "genbank"):
    accession = str(seq_record.id)
    accession_name_list.append(accession)
    gb_entire_sequence = list(seq_record.seq)
    gb_entire_sequence_lower = [x.lower() for x in gb_entire_sequence]
    gb_entire_sequence_joined = ''.join(gb_entire_sequence_lower)
    for feat in range(len(seq_record.features)):
        feature_type = seq_record.features[feat].type
        if feature_type == 'gene':
            try:
                feature_start_zero_based_numbering = seq_record.features[feat].location.nofuzzy_start
                feature_end_zero_based_numbering = seq_record.features[feat].location.nofuzzy_end
                feature_strand = seq_record.features[feat].strand
                gene_id = str(seq_record.features[feat].qualifiers['gene'][0])
                sequence_slice = gb_entire_sequence_joined[feature_start_zero_based_numbering:feature_end_zero_based_numbering]
                if feature_strand == 1:
                    gene_array.at[accession,gene_id] = sequence_slice
                if feature_strand == -1:
                    sequence_slice_BP = Seq(sequence_slice)
                    sequence_slice_rvscomp = str(sequence_slice_BP.reverse_complement())
                    gene_array.at[accession,gene_id] = sequence_slice_rvscomp
            except KeyError:
                pass
            
# Count and record instances of each gene in the dataframe

gene_instances_list = []
for gene in gene_list_rem_dups:
    gene_to_append = gene, len(gene_array)-gene_array[gene].isnull().sum()
    gene_instances_list.append(gene_to_append) 
gene_instances_list_DF = pd.DataFrame(gene_instances_list,columns = ['Gene','Instances'])

# Generate and save a description dataframe

description_array = pd.DataFrame(index=accession_name_list,columns=['Organism'])
description_array.index.name = 'Accession'

accession_name_list = []
for seq_record in SeqIO.parse(genbank_file_to_use, "genbank"):
    accession = str(seq_record.id)
    for feat in range(len(seq_record.features)):
        feature_type = seq_record.features[feat].type
        if feature_type == 'source':
            try:
                organism = str(seq_record.features[feat].qualifiers['organism'][0])
                description_array.at[accession,'Organism'] = organism
            except KeyError:
                pass

description_array['Sequence_name'] = description_array.index + '_' + description_array['Organism']
description_array['Sequence_name'] = description_array['Sequence_name'].replace('[ .]', '_', regex=True)
description_array['Sequence_name'] = description_array['Sequence_name'].replace('[^a-zA-Z0-9_]', '', regex=True)

# Delete undesired columns based upon chosen genes

protein_array = protein_array[chosen_genes]
gene_array = gene_array[chosen_genes]

# Protein dataframe columns rename
    
protein_array_column_names_list = list(protein_array.columns.values.tolist())
protein_array_column_names_list_rename = [name + '_prot' for name in protein_array_column_names_list]
for i in range(len(protein_array_column_names_list)):
    protein_array.rename(columns = {protein_array_column_names_list[i]:protein_array_column_names_list_rename[i]}, inplace = True) 
    
# Dataframes concatenate sequences and form combined array

gene_array['Concatenate_Coding_DNA'] = ''

for gene in chosen_genes:
    
    gene_array['Concatenate_Coding_DNA'] = gene_array['Concatenate_Coding_DNA'] + gene_array[gene]

proteins_to_concatenate = []
for gene in chosen_genes:
    proteins_to_concatenate.append(gene + '_prot')

protein_array['Concatenate_Protein'] = ''
for protein in proteins_to_concatenate:
    protein_array['Concatenate_Protein'] = protein_array['Concatenate_Protein'] + protein_array[protein]
    
combined_array = description_array.join(protein_array).join(gene_array)

# Remove any entries from gene dataframe with empty cell

combined_array.replace('', np.nan, inplace=True)
combined_array.dropna(inplace=True)

# Remove duplicates from concatenates columns

combined_array.drop_duplicates(subset=['Concatenate_Coding_DNA'], keep='first',inplace=True)
combined_array.drop_duplicates(subset=['Concatenate_Protein'], keep='first',inplace=True)

# Choose sequence_name as index

combined_array.reset_index(drop=True, inplace=True)

# Save combined array

combined_array.to_csv(file_prefix + 'combined_input_array.csv', index=True)

print(str(len(combined_array)) + ' samples to be analyzed after pre-processing.')

for gene in chosen_genes:
    ofile = open(file_prefix + gene + '.fasta', "w")
    for seqi in range(len(combined_array)):
        ofile.write(">" + combined_array.at[seqi,'Sequence_name'] + "\n" + combined_array.at[seqi,gene] + "\n")
    ofile.close()
    
# Free memory

del combined_array
del gene_array
del protein_array

output_list = []

for gene in chosen_genes:
    
    # Initialize lists for FASTA information

    record_x_toward_seq_dataframe = []
    sequence_records = []
    alignment_record_name_list = []
    
    # Load input FASTA into dataframe
    
    alignfile = file_prefix + gene + '.fasta'
    
    for record in SeqIO.parse(alignfile,"fasta"):
        
        if record.name != 'NC_044125_1_Anolis_punctatus':

            record_string = str(record.seq)
            translated_reference_sequence_concatenate_Seq = record.seq.translate(table=codon_table)
            translated_reference_sequence_concatenate_STR = str(translated_reference_sequence_concatenate_Seq)

            for i in range(len(translated_reference_sequence_concatenate_STR)):
                first = i*3
                second = i*3+1
                third = i*3+2

                amino_acid = translated_reference_sequence_concatenate_STR[i]


                if codon_table == 2 and record_string[first].upper() == 'T' and record_string[second].upper() == 'C':
                    amino_acid_to_write = 'S1'
                
                elif codon_table == 2 and record_string[first].upper() == 'A' and record_string[second].upper() == 'G' and record_string[third].upper() == 'T':
                    amino_acid_to_write = 'S2'
                
                elif codon_table == 2 and record_string[first].upper() == 'A' and record_string[second].upper() == 'G' and record_string[third].upper() == 'C':
                    amino_acid_to_write = 'S2'
                
                elif codon_table == 2 and record_string[first].upper() == 'C' and record_string[second].upper() == 'T':
                    amino_acid_to_write = 'L2'
                
                elif codon_table == 2 and record_string[first].upper() == 'T' and record_string[second].upper() == 'T' and record_string[third].upper() == 'A':
                    amino_acid_to_write = 'L1'
                
                elif codon_table == 2 and record_string[first].upper() == 'T' and record_string[second].upper() == 'T' and record_string[third].upper() == 'G':
                    amino_acid_to_write = 'L1'
                
                else:
                    amino_acid_to_write = amino_acid

                if amino_acid_to_write in ['A','C','D','E','F','G','H','I','K','L1','L2','M','N','P','Q','R','S1','S2','T','V','W','Y']:
                    nucleotide_at_third_position = record_string[third]
                    output_data = record.name,gene,amino_acid_to_write,i+1,(i+1)*3,nucleotide_at_third_position.upper()
                    output_list.append(output_data)
                    
third_codon_instances_by_AA_and_acc = pd.DataFrame(output_list,columns = ['Sequence_record','Gene','Amino_acid','Amino_acid_position','Nucleotide_position_third_codon_position','Nucleotide'])

third_codon_instances_by_AA_and_acc

counts_of_amino_acids = third_codon_instances_by_AA_and_acc['Amino_acid'].value_counts()
counts_of_amino_acids.to_csv(file_prefix + 'count_of_amino_acids_tested_across_input_file.csv')

third_position_by_AA_DF_norm = pd.DataFrame(index = ['A','C','G','T'], columns = ['A','C','D','E','F','G','H','I','K','L1','L2','M','N','P','Q','R','S1','S2','T','V','W','Y'])
third_position_by_AA_DF_notnorm = pd.DataFrame(index = ['A','C','G','T'], columns = ['A','C','D','E','F','G','H','I','K','L1','L2','M','N','P','Q','R','S1','S2','T','V','W','Y'])

for amino_acid in ['A','C','D','E','F','G','H','I','K','L1','L2','M','N','P','Q','R','S1','S2','T','V','W','Y']:
    AA_dataframe = third_codon_instances_by_AA_and_acc[third_codon_instances_by_AA_and_acc['Amino_acid'] == amino_acid]
    nucleotide_fraction = AA_dataframe['Nucleotide']
    nucleotide_fractionACGT = nucleotide_fraction[nucleotide_fraction.isin(['A','C','G','T'])]
    third_position_by_AA_DF_norm[amino_acid] = nucleotide_fractionACGT.value_counts(normalize=True)
    third_position_by_AA_DF_notnorm[amino_acid] = nucleotide_fractionACGT.value_counts(normalize=False)

third_position_by_AA_DF_norm.to_csv(file_prefix + 'nucleotide_frequencies_at_third_codon_position_by_AA.csv')
third_position_by_AA_DF_notnorm.to_csv(file_prefix + 'nucleotide_hits_at_third_codon_position_by_AA.csv')



