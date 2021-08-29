#!/usr/bin/env python
# coding: utf-8

# Funding for this project was received from the European Research Council and the Sigrid Jusélius Foundation.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.dunn@helsinki.fi
# License: GPLv3
# Version: B
# Please cite DOI: 10.1101/2021.05.12.443844

import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
import math

# Read TSS data by P3 type

two_fold_AG = pd.read_csv('two_fold_AG_P3_ALL.csv')
two_fold_AG = two_fold_AG.rename(columns = {'Unnamed: 7': 'Amino_merge'}, inplace = False)

two_fold_CT = pd.read_csv('two_fold_CT_P3_ALL.csv')
two_fold_CT = two_fold_CT.rename(columns = {'Unnamed: 8': 'Amino_merge'}, inplace = False)

four_fold = pd.read_csv('four_fold_P3_ALL.csv')
four_fold = four_fold.rename(columns = {'Unnamed: 8': 'Amino_merge'}, inplace = False)

# Save TSS data by amino acid

four_fold_amino_acids = ['G','P','V','S1','A','T','R']

for amino_acid_four_fold in four_fold_amino_acids:
    subset_four_fold = four_fold[four_fold['Amino_acid'] == amino_acid_four_fold]
    subset_four_fold.to_csv(amino_acid_four_fold + '_degenerate.csv')

two_fold_AG_amino_acids = ['E','M','Q','W','K']

for amino_acid_two_fold_AG in two_fold_AG_amino_acids:
    subset_two_fold_AG = two_fold_AG[two_fold_AG['Amino_acid'] == amino_acid_two_fold_AG]
    subset_two_fold_AG.to_csv(amino_acid_two_fold_AG + '_degenerate.csv')

two_fold_CT_amino_acids = ['Y','N','D','H','I','F','C','S2']

for amino_acid_two_fold_CT in two_fold_CT_amino_acids:
    subset_two_fold_CT = two_fold_CT[two_fold_CT['Amino_acid'] == amino_acid_two_fold_CT]
    subset_two_fold_CT.to_csv(amino_acid_two_fold_CT + '_degenerate.csv')

# P-value and D-stat for AA versus TSS (four-fold degenerate)

four_fold_AA_vs_TSS_KS_Dstat_DF = pd.DataFrame(index = four_fold_amino_acids,columns = four_fold_amino_acids)
four_fold_AA_vs_TSS_KS_Pvalue_DF = pd.DataFrame(index = four_fold_amino_acids,columns = four_fold_amino_acids)

for i in four_fold_amino_acids:
    for j in four_fold_amino_acids:
        test_column_1 = four_fold[four_fold['Amino_acid'] == i]['TSS']
        test_column_2 = four_fold[four_fold['Amino_acid'] == j]['TSS']
        if len(test_column_1) > 0 and len(test_column_2) > 0:
            teststatistic,pvalue = ks_2samp(test_column_1, test_column_2)
            
            four_fold_AA_vs_TSS_KS_Dstat_DF.at[i,j] = teststatistic
            four_fold_AA_vs_TSS_KS_Pvalue_DF.at[i,j] = pvalue

four_fold_AA_vs_TSS_KS_Dstat_DF.to_csv('four_fold_compare_AA_to_TSS_by_KS_test_Dstat.csv')
four_fold_AA_vs_TSS_KS_Pvalue_DF.to_csv('four_fold_compare_AA_to_TSS_by_KS_test_uncorrected_Pvalue.csv')
four_fold_AA_vs_TSS_KS_Pvalue_corrected_DF = four_fold_AA_vs_TSS_KS_Pvalue_DF * 21 # Twenty-one TSS comparisons performed across seven amino acid choices
four_fold_AA_vs_TSS_KS_Pvalue_corrected_DF.to_csv('four_fold_compare_AA_to_TSS_by_KS_test_corrected_Pvalue.csv')

# P-value and D-stat for AA versus TSS (two-fold AG)

two_fold_AG_AA_vs_TSS_KS_Dstat_DF = pd.DataFrame(index = two_fold_AG_amino_acids,columns = two_fold_AG_amino_acids)
two_fold_AG_AA_vs_TSS_KS_Pvalue_DF = pd.DataFrame(index = two_fold_AG_amino_acids,columns = two_fold_AG_amino_acids)

for i in two_fold_AG_amino_acids:
    for j in two_fold_AG_amino_acids:
        test_column_1 = two_fold_AG[two_fold_AG['Amino_acid'] == i]['TSS']
        test_column_2 = two_fold_AG[two_fold_AG['Amino_acid'] == j]['TSS']
        if len(test_column_1) > 0 and len(test_column_2) > 0:
            teststatistic,pvalue = ks_2samp(test_column_1, test_column_2)
            
            two_fold_AG_AA_vs_TSS_KS_Dstat_DF.at[i,j] = teststatistic
            two_fold_AG_AA_vs_TSS_KS_Pvalue_DF.at[i,j] = pvalue

two_fold_AG_AA_vs_TSS_KS_Dstat_DF.to_csv('two_fold_AG_compare_AA_to_TSS_by_KS_test_Dstat.csv')
two_fold_AG_AA_vs_TSS_KS_Pvalue_DF.to_csv('two_fold_AG_compare_AA_to_TSS_by_KS_test_uncorrected_Pvalue.csv')
two_fold_AG_AA_vs_TSS_KS_Pvalue_corrected_DF = two_fold_AG_AA_vs_TSS_KS_Pvalue_DF * 10 # Ten TSS comparisons performed across seven amino acid choices
two_fold_AG_AA_vs_TSS_KS_Pvalue_corrected_DF.to_csv('two_fold_AG_compare_AA_to_TSS_by_KS_test_corrected_Pvalue.csv')

# P-value and D-stat for AA versus TSS (two-fold CT)

two_fold_CT_AA_vs_TSS_KS_Dstat_DF = pd.DataFrame(index = two_fold_CT_amino_acids,columns = two_fold_CT_amino_acids)
two_fold_CT_AA_vs_TSS_KS_Pvalue_DF = pd.DataFrame(index = two_fold_CT_amino_acids,columns = two_fold_CT_amino_acids)

for i in two_fold_CT_amino_acids:
    for j in two_fold_CT_amino_acids:
        test_column_1 = two_fold_CT[two_fold_CT['Amino_acid'] == i]['TSS']
        test_column_2 = two_fold_CT[two_fold_CT['Amino_acid'] == j]['TSS']
        if len(test_column_1) > 0 and len(test_column_2) > 0:
            teststatistic,pvalue = ks_2samp(test_column_1, test_column_2)
            
            two_fold_CT_AA_vs_TSS_KS_Dstat_DF.at[i,j] = teststatistic
            two_fold_CT_AA_vs_TSS_KS_Pvalue_DF.at[i,j] = pvalue

two_fold_CT_AA_vs_TSS_KS_Dstat_DF.to_csv('two_fold_CT_compare_AA_to_TSS_by_KS_test_Dstat.csv')
two_fold_CT_AA_vs_TSS_KS_Pvalue_DF.to_csv('two_fold_CT_compare_AA_to_TSS_by_KS_test_uncorrected_Pvalue.csv')
two_fold_CT_AA_vs_TSS_KS_Pvalue_corrected_DF = two_fold_CT_AA_vs_TSS_KS_Pvalue_DF * 28 # Twenty-eight TSS comparisons performed across seven amino acid choices
two_fold_CT_AA_vs_TSS_KS_Pvalue_corrected_DF.to_csv('two_fold_CT_compare_AA_to_TSS_by_KS_test_corrected_Pvalue.csv')

# Load processed HelixMTdb results not collapsed by amino acid change

helixMTdb_NOT_grouped_by_AA = pd.read_csv('HelixMTdb_results_NOT_grouped_by_AA_change_single_nuc_sub_from_human_ref.csv')

# Separate HelixMTdb into synonymous substitutions and non-synonymous substitutions

helixMTdb_NOT_grouped_by_AA_syn = helixMTdb_NOT_grouped_by_AA[helixMTdb_NOT_grouped_by_AA['reference_amino_acid'] == helixMTdb_NOT_grouped_by_AA['mutant_amino_acid']]
helixMTdb_NOT_grouped_by_AA_NONsyn = helixMTdb_NOT_grouped_by_AA[helixMTdb_NOT_grouped_by_AA['reference_amino_acid'] != helixMTdb_NOT_grouped_by_AA['mutant_amino_acid']]
del helixMTdb_NOT_grouped_by_AA

# Place synonymous variants into classes by frequency


helixMTdb_NOT_grouped_by_AA_syn['Amino_merge'] = helixMTdb_NOT_grouped_by_AA_syn['Amino_acid_change'].str[:-1].astype(str)
helixMTdb_NOT_grouped_by_AA_syn['counts_total_logbase2'] = np.log2(helixMTdb_NOT_grouped_by_AA_syn['counts_total'])
helixMTdb_NOT_grouped_by_AA_syn['counts_total_logbase10'] = np.log10(helixMTdb_NOT_grouped_by_AA_syn['counts_total'])
helixMTdb_NOT_grouped_by_AA_syn['allele_frequency'] = helixMTdb_NOT_grouped_by_AA_syn['counts_total'] / 195983 # (total samples in HelixMTdb)

helixMTdb_NOT_grouped_by_AA_syn['allele_status'] = ''
helixMTdb_NOT_grouped_by_AA_syn['allele_status'] = np.where(helixMTdb_NOT_grouped_by_AA_syn['allele_frequency'] == 0, 'Absent', helixMTdb_NOT_grouped_by_AA_syn['allele_status'])
helixMTdb_NOT_grouped_by_AA_syn['allele_status'] = np.where(helixMTdb_NOT_grouped_by_AA_syn['allele_frequency'] >= .05, 'Common', helixMTdb_NOT_grouped_by_AA_syn['allele_status'])
helixMTdb_NOT_grouped_by_AA_syn['allele_status'] = np.where((helixMTdb_NOT_grouped_by_AA_syn['allele_frequency'] < .05) & (helixMTdb_NOT_grouped_by_AA_syn['allele_frequency'] >= .01), 'Low-frequency', helixMTdb_NOT_grouped_by_AA_syn['allele_status'])
helixMTdb_NOT_grouped_by_AA_syn['allele_status'] = np.where((helixMTdb_NOT_grouped_by_AA_syn['allele_frequency'] < .01) & (helixMTdb_NOT_grouped_by_AA_syn['allele_frequency'] >= .0001), 'Rare', helixMTdb_NOT_grouped_by_AA_syn['allele_status'])
helixMTdb_NOT_grouped_by_AA_syn['allele_status'] = np.where((helixMTdb_NOT_grouped_by_AA_syn['allele_frequency'] < .0001) & (helixMTdb_NOT_grouped_by_AA_syn['allele_frequency'] > 0), 'Ultra-rare',helixMTdb_NOT_grouped_by_AA_syn['allele_status'])

# Place non-synonymous variants into classes by frequency


helixMTdb_NOT_grouped_by_AA_NONsyn['Amino_merge'] = helixMTdb_NOT_grouped_by_AA_NONsyn['Amino_acid_change'].str[:-1].astype(str)
helixMTdb_NOT_grouped_by_AA_NONsyn['counts_total_logbase2'] = np.log2(helixMTdb_NOT_grouped_by_AA_NONsyn['counts_total'])
helixMTdb_NOT_grouped_by_AA_NONsyn['counts_total_logbase10'] = np.log10(helixMTdb_NOT_grouped_by_AA_NONsyn['counts_total'])
helixMTdb_NOT_grouped_by_AA_NONsyn['allele_frequency'] = helixMTdb_NOT_grouped_by_AA_NONsyn['counts_total'] / 195983 # (total samples in HelixMTdb)

helixMTdb_NOT_grouped_by_AA_NONsyn['allele_status'] = ''
helixMTdb_NOT_grouped_by_AA_NONsyn['allele_status'] = np.where(helixMTdb_NOT_grouped_by_AA_NONsyn['allele_frequency'] == 0, 'Absent', helixMTdb_NOT_grouped_by_AA_NONsyn['allele_status'])
helixMTdb_NOT_grouped_by_AA_NONsyn['allele_status'] = np.where(helixMTdb_NOT_grouped_by_AA_NONsyn['allele_frequency'] >= .05, 'Common', helixMTdb_NOT_grouped_by_AA_NONsyn['allele_status'])
helixMTdb_NOT_grouped_by_AA_NONsyn['allele_status'] = np.where((helixMTdb_NOT_grouped_by_AA_NONsyn['allele_frequency'] < .05) & (helixMTdb_NOT_grouped_by_AA_NONsyn['allele_frequency'] >= .01), 'Low-frequency', helixMTdb_NOT_grouped_by_AA_NONsyn['allele_status'])
helixMTdb_NOT_grouped_by_AA_NONsyn['allele_status'] = np.where((helixMTdb_NOT_grouped_by_AA_NONsyn['allele_frequency'] < .01) & (helixMTdb_NOT_grouped_by_AA_NONsyn['allele_frequency'] >= .0001), 'Rare', helixMTdb_NOT_grouped_by_AA_NONsyn['allele_status'])
helixMTdb_NOT_grouped_by_AA_NONsyn['allele_status'] = np.where((helixMTdb_NOT_grouped_by_AA_NONsyn['allele_frequency'] < .0001) & (helixMTdb_NOT_grouped_by_AA_NONsyn['allele_frequency'] > 0), 'Ultra-rare',helixMTdb_NOT_grouped_by_AA_NONsyn['allele_status'])

# Generate CSV of synonymous variants based upon counts total

try: 
    del helixMTdb_NOT_grouped_by_AA_syn ['Unnamed: 0']
except:
    pass
helixMTdb_NOT_grouped_by_AA_syn.reset_index(inplace=True,drop=True)
del helixMTdb_NOT_grouped_by_AA_syn['Protein']

helixMTdb_NOT_grouped_by_AA_syn.to_csv('HelixMTdb_synonymous_variant_set.csv')

helixMTdb_NOT_grouped_by_AA_syn_more_than_zero_count = helixMTdb_NOT_grouped_by_AA_syn[helixMTdb_NOT_grouped_by_AA_syn['counts_total'] > 0]

helixMTdb_NOT_grouped_by_AA_syn_more_than_zero_count.to_csv('HelixMTdb_synonymous_variant_set_more_than_zero.csv')

helixMTdb_NOT_grouped_by_AA_syn_more_than_ten_count = helixMTdb_NOT_grouped_by_AA_syn[helixMTdb_NOT_grouped_by_AA_syn['counts_total'] > 9]

helixMTdb_NOT_grouped_by_AA_syn_more_than_ten_count.to_csv('HelixMTdb_synonymous_variant_set_more_than_nine.csv')

# Generate CSV of non-synonymous variants based upon counts total

try: 
    del helixMTdb_NOT_grouped_by_AA_NONsyn ['Unnamed: 0']
except:
    pass
helixMTdb_NOT_grouped_by_AA_NONsyn.reset_index(inplace=True,drop=True)
del helixMTdb_NOT_grouped_by_AA_NONsyn['Protein']

helixMTdb_NOT_grouped_by_AA_NONsyn.to_csv('HelixMTdb_NONsynonymous_variant_set.csv')

helixMTdb_NOT_grouped_by_AA_NONsyn_more_than_zero_count = helixMTdb_NOT_grouped_by_AA_NONsyn[helixMTdb_NOT_grouped_by_AA_NONsyn['counts_total'] > 0]

helixMTdb_NOT_grouped_by_AA_NONsyn_more_than_zero_count.to_csv('HelixMTdb_NONsynonymous_variant_set_more_than_zero.csv')

helixMTdb_NOT_grouped_by_AA_NONsyn_more_than_ten_count = helixMTdb_NOT_grouped_by_AA_NONsyn[helixMTdb_NOT_grouped_by_AA_NONsyn['counts_total'] > 9]

helixMTdb_NOT_grouped_by_AA_NONsyn_more_than_ten_count.to_csv('HelixMTdb_NONsynonymous_variant_set_more_than_nine.csv')

# Generate CSV of I-P3s based upon allele frequency and four-fold degeneracy

four_fold_w_HelixMTdb_counts = four_fold.set_index('Amino_merge').join(helixMTdb_NOT_grouped_by_AA_syn.set_index('Amino_merge'))
four_fold_w_HelixMTdb_counts.sort_values(by=['allele_frequency','TSS'],ascending=True,inplace=True)
four_fold_w_HelixMTdb_counts.to_csv('four_fold_synonymous_w_HelixMTdb_counts.csv')

# Make dataframe of above with counts total > 0 (four-fold degeneracy)

four_fold_w_HelixMTdb_counts_more_than_zero = four_fold_w_HelixMTdb_counts[four_fold_w_HelixMTdb_counts['counts_total'] > 0]
four_fold_w_HelixMTdb_counts_more_than_zero['mutation_twochar'] = four_fold_w_HelixMTdb_counts_more_than_zero['reference_nucleotide'].map(str) + four_fold_w_HelixMTdb_counts_more_than_zero['mutant_nucleotide'].map(str)
four_fold_w_HelixMTdb_counts_more_than_zero['mutation_type'] = 'Transversion'
four_fold_w_HelixMTdb_counts_more_than_zero['mutation_type'] = np.where(four_fold_w_HelixMTdb_counts_more_than_zero['mutation_twochar'] == 'AG','Transition',four_fold_w_HelixMTdb_counts_more_than_zero['mutation_type'])
four_fold_w_HelixMTdb_counts_more_than_zero['mutation_type'] = np.where(four_fold_w_HelixMTdb_counts_more_than_zero['mutation_twochar'] == 'GA','Transition',four_fold_w_HelixMTdb_counts_more_than_zero['mutation_type'])
four_fold_w_HelixMTdb_counts_more_than_zero['mutation_type'] = np.where(four_fold_w_HelixMTdb_counts_more_than_zero['mutation_twochar'] == 'CT','Transition',four_fold_w_HelixMTdb_counts_more_than_zero['mutation_type'])
four_fold_w_HelixMTdb_counts_more_than_zero['mutation_type'] = np.where(four_fold_w_HelixMTdb_counts_more_than_zero['mutation_twochar'] == 'TC','Transition',four_fold_w_HelixMTdb_counts_more_than_zero['mutation_type'])
four_fold_w_HelixMTdb_counts_more_than_zero.sort_values(by='counts_total',ascending=True, inplace=True)
four_fold_w_HelixMTdb_counts_more_than_zero.to_csv('four_fold_synonymous_w_HelixMTdb_counts_more_than_zero_sorted_by_counts_total.csv')

# Generate CSV of I-P3s based upon allele frequency and two-fold AG degeneracy

two_fold_AG_w_HelixMTdb_counts = two_fold_AG.set_index('Amino_merge').join(helixMTdb_NOT_grouped_by_AA_syn.set_index('Amino_merge'))
two_fold_AG_w_HelixMTdb_counts.sort_values(by=['allele_frequency','TSS'],ascending=True,inplace=True)
two_fold_AG_w_HelixMTdb_counts.to_csv('two_fold_AG_synonymous_w_HelixMTdb_counts.csv')

# Make dataframe of above with counts total > 0 (two-fold AG degeneracy)

two_fold_AG_w_HelixMTdb_counts_more_than_zero = two_fold_AG_w_HelixMTdb_counts[two_fold_AG_w_HelixMTdb_counts['counts_total'] > 0]
two_fold_AG_w_HelixMTdb_counts_more_than_zero.sort_values(by='counts_total', ascending=True, inplace=True)
two_fold_AG_w_HelixMTdb_counts_more_than_zero.to_csv('two_fold_AG_synonymous_w_HelixMTdb_counts_more_than_zero_sorted_by_counts_total.csv')

# Generate CSV of I-P3s based upon allele frequency and two-fold CT degeneracy

two_fold_CT_w_HelixMTdb_counts = two_fold_CT.set_index('Amino_merge').join(helixMTdb_NOT_grouped_by_AA_syn.set_index('Amino_merge'))
two_fold_CT_w_HelixMTdb_counts.sort_values(by=['allele_frequency','TSS'],ascending=True,inplace=True)
two_fold_CT_w_HelixMTdb_counts.to_csv('two_fold_CT_synonymous_w_HelixMTdb_counts.csv')

# Make dataframe of above with counts total > 0 (two-fold CT degeneracy)

two_fold_CT_w_HelixMTdb_counts_more_than_zero = two_fold_CT_w_HelixMTdb_counts[two_fold_CT_w_HelixMTdb_counts['counts_total'] > 0]
two_fold_CT_w_HelixMTdb_counts_more_than_zero.sort_values(by='counts_total', ascending=True, inplace=True)
two_fold_CT_w_HelixMTdb_counts_more_than_zero.to_csv('two_fold_CT_synonymous_w_HelixMTdb_counts_more_than_zero_sorted_by_counts_total.csv')

# P-value and D-stat for four-fold class versus TSS [Kolmogorov-Smirnov analysis]

variant_classes = ['Absent','Ultra-rare','Rare','Low-frequency','Common']

four_fold_KS_Dstat_DF = pd.DataFrame(index = variant_classes,columns = variant_classes)
four_fold_KS_Pvalue_DF = pd.DataFrame(index = variant_classes,columns = variant_classes)

for i in variant_classes:
    for j in variant_classes:
        test_column_1 = four_fold_w_HelixMTdb_counts[four_fold_w_HelixMTdb_counts['allele_status'] == i]['TSS']
        test_column_2 = four_fold_w_HelixMTdb_counts[four_fold_w_HelixMTdb_counts['allele_status'] == j]['TSS']
        if len(test_column_1) > 0 and len(test_column_2) > 0:
            teststatistic,pvalue = ks_2samp(test_column_1, test_column_2)
            
            four_fold_KS_Dstat_DF.at[i,j] = teststatistic
            four_fold_KS_Pvalue_DF.at[i,j] = pvalue

four_fold_KS_Dstat_DF.to_csv('four_fold_compare_variant_class_to_TSS_by_KS_test_Dstat.csv')
four_fold_KS_Pvalue_DF.to_csv('four_fold_compare_variant_class_to_TSS_by_KS_test_uncorrected_Pvalue.csv')
four_fold_KS_Pvalue_corrected_DF = four_fold_KS_Pvalue_DF * 6 # six TSS comparisons performed across four classes COMMON NOT INCLUDED DUE TO LOW NUMBERS
four_fold_KS_Pvalue_corrected_DF.to_csv('four_fold_compare_variant_class_to_TSS_by_KS_test_corrected_Pvalue.csv')

# P-value and D-stat for two-fold AG class versus TSS [Kolmogorov-Smirnov analysis]

variant_classes = ['Absent','Ultra-rare','Rare','Low-frequency','Common']

two_fold_AG_KS_Dstat_DF = pd.DataFrame(index = variant_classes,columns = variant_classes)
two_fold_AG_KS_Pvalue_DF = pd.DataFrame(index = variant_classes,columns = variant_classes)

for i in variant_classes:
    for j in variant_classes:
        test_column_1 = two_fold_AG_w_HelixMTdb_counts[two_fold_AG_w_HelixMTdb_counts['allele_status'] == i]['TSS']
        test_column_2 = two_fold_AG_w_HelixMTdb_counts[two_fold_AG_w_HelixMTdb_counts['allele_status'] == j]['TSS']
        if len(test_column_1) > 0 and len(test_column_2) > 0:
            teststatistic,pvalue = ks_2samp(test_column_1, test_column_2)
            
            two_fold_AG_KS_Dstat_DF.at[i,j] = teststatistic
            two_fold_AG_KS_Pvalue_DF.at[i,j] = pvalue

two_fold_AG_KS_Dstat_DF.to_csv('two_fold_AG_compare_variant_class_to_TSS_by_KS_test_Dstat.csv')
two_fold_AG_KS_Pvalue_DF.to_csv('two_fold_AG_compare_variant_class_to_TSS_by_KS_test_uncorrected_Pvalue.csv')
two_fold_AG_KS_Pvalue_corrected_DF = two_fold_AG_KS_Pvalue_DF * 6 # six TSS comparisons performed across four classes COMMON NOT INCLUDED DUE TO LOW NUMBERS
two_fold_AG_KS_Pvalue_corrected_DF.to_csv('two_fold_AG_compare_variant_class_to_TSS_by_KS_test_corrected_Pvalue.csv')

# P-value and D-stat for two-fold CT class versus TSS [Kolmogorov-Smirnov analysis]

variant_classes = ['Absent','Ultra-rare','Rare','Low-frequency','Common']

two_fold_CT_KS_Dstat_DF = pd.DataFrame(index = variant_classes,columns = variant_classes)
two_fold_CT_KS_Pvalue_DF = pd.DataFrame(index = variant_classes,columns = variant_classes)

for i in variant_classes:
    for j in variant_classes:
        test_column_1 = two_fold_CT_w_HelixMTdb_counts[two_fold_CT_w_HelixMTdb_counts['allele_status'] == i]['TSS']
        test_column_2 = two_fold_CT_w_HelixMTdb_counts[two_fold_CT_w_HelixMTdb_counts['allele_status'] == j]['TSS']
        if len(test_column_1) > 0 and len(test_column_2) > 0:
            teststatistic,pvalue = ks_2samp(test_column_1, test_column_2)
            
            two_fold_CT_KS_Dstat_DF.at[i,j] = teststatistic
            two_fold_CT_KS_Pvalue_DF.at[i,j] = pvalue

two_fold_CT_KS_Dstat_DF.to_csv('two_fold_CT_compare_variant_class_to_TSS_by_KS_test_Dstat.csv')
two_fold_CT_KS_Pvalue_DF.to_csv('two_fold_CT_compare_variant_class_to_TSS_by_KS_test_uncorrected_Pvalue.csv')
two_fold_CT_KS_Pvalue_corrected_DF = two_fold_CT_KS_Pvalue_DF * 6 # six TSS comparisons performed across four classes COMMON NOT INCLUDED DUE TO LOW NUMBERS
two_fold_CT_KS_Pvalue_corrected_DF.to_csv('two_fold_CT_compare_variant_class_to_TSS_by_KS_test_corrected_Pvalue.csv')
