#!/usr/bin/env python
# coding: utf-8 

# Funding for this project was received from the European Research Council and the Sigrid Jusélius Foundation.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.dunn@helsinki.fi
# License: GPLv3
# Version: A

# Import dependencies

import pandas as pd
import numpy as np

# Prepare dataframe

HelixMTdb_raw = pd.read_csv('HelixMTdb_20200327.tsv', sep='\t')

HelixMTdb_raw['alleles'].replace('\[','',regex=True,inplace=True)
HelixMTdb_raw['alleles'].replace(']','',regex=True,inplace=True)
HelixMTdb_raw['alleles'].replace('\"','',regex=True,inplace=True)
HelixMTdb_raw['locus'].replace('chrM:','',regex=True,inplace=True)

HelixMTdb_raw[['reference_nucleotide_TEMP','mutant_nucleotide_TEMP','X']] = HelixMTdb_raw.alleles.str.split(',',expand=True)

HelixMTdb_raw_remove_indel = HelixMTdb_raw[HelixMTdb_raw['mutant_nucleotide_TEMP'].str.len() == 1]

HelixMTdb_raw_remove_indel['Nucleotide_change'] = HelixMTdb_raw_remove_indel['reference_nucleotide_TEMP'] + HelixMTdb_raw_remove_indel['locus'].map(str) + + HelixMTdb_raw_remove_indel['mutant_nucleotide_TEMP']

HelixMTdb_output_DF = HelixMTdb_raw_remove_indel[['Nucleotide_change', 'gene', 'reference_nucleotide_TEMP', 'mutant_nucleotide_TEMP','counts_hom','counts_het']]

# Joining substitution set with amino acid changes

possible_amino_acid_substitutions = pd.read_csv('NC_012920_final_dataframe_all_single_nucleotide_substitutions_as_protein_changes.csv')
total_array = possible_amino_acid_substitutions.merge(HelixMTdb_output_DF,on=['Nucleotide_change'],how='left')

total_array['counts_hom'] = total_array['counts_hom'].fillna(0)
total_array['counts_het'] = total_array['counts_het'].fillna(0)

total_array['counts_total'] = total_array['counts_hom'] + total_array['counts_het']
total_array['het_percent'] = total_array['counts_het'] / total_array['counts_total'] * 100
total_array[['prot_to_del','change']] = total_array['Amino_acid_change'].str.split('_',expand=True)
total_array['reference_amino_acid'] = total_array['change'].str[:1]
total_array['mutant_amino_acid'] = total_array['change'].str[-1]
total_array['change_to_amino_acid'] = np.where(total_array['reference_amino_acid'] == total_array['mutant_amino_acid'],'No','Yes')


del total_array['prot_to_del']
del total_array['reference_nucleotide_TEMP']
del total_array['mutant_nucleotide_TEMP']
del total_array['gene']

total_array['reference_nucleotide'] = total_array['Nucleotide_change'].str[:1]
total_array['mutant_nucleotide'] = total_array['Nucleotide_change'].str[-1]

total_array_group_by_AAchange = total_array.groupby(by='Amino_acid_change').sum()
del total_array_group_by_AAchange['het_percent']
del total_array_group_by_AAchange['counts_total']
total_array_group_by_AAchange['counts_total'] = total_array_group_by_AAchange['counts_hom'] + total_array_group_by_AAchange['counts_het']
total_array_group_by_AAchange['merged_het_percent'] = total_array_group_by_AAchange['counts_het'] / total_array_group_by_AAchange['counts_total'] * 100

total_array_group_by_AAchange = total_array_group_by_AAchange.reset_index()
total_array_group_by_AAchange[['protein','change']] = total_array_group_by_AAchange['Amino_acid_change'].str.split('_',expand=True)

total_array_group_by_AAchange['reference_amino_acid'] = total_array_group_by_AAchange['change'].str[:1]
total_array_group_by_AAchange['mutant_amino_acid'] = total_array_group_by_AAchange['change'].str[-1]
total_array_group_by_AAchange['change_to_amino_acid'] = np.where(total_array_group_by_AAchange['reference_amino_acid'] == total_array_group_by_AAchange['mutant_amino_acid'],'No','Yes')

total_array.to_csv('HelixMTdb_results_NOT_grouped_by_AA_change_single_nuc_sub_from_human_ref.csv')
total_array_group_by_AAchange.to_csv('HelixMTdb_results_grouped_by_AA_change_single_nuc_sub_from_human_ref.csv')
