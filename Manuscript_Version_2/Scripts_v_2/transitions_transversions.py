#!/usr/bin/env python
# coding: utf-8

# Funding for this project was received from the European Research Council and the Sigrid Jusélius Foundation.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.dunn@helsinki.fi
# License: GPLv3
# Version: A
# Please cite DOI: 10.1101/2021.05.12.443844

# Import dependencies

import pandas as pd
import numpy as np

helix_single_sub = pd.read_csv('HelixMTdb_results_NOT_grouped_by_AA_change_single_nuc_sub_from_human_ref.csv')

# Set pyrimidines, purines, and detect transitions and transversions

helix_single_sub['ref_nuc_type'] = np.where((helix_single_sub['reference_nucleotide'] == 'A') | (helix_single_sub['reference_nucleotide'] == 'G'), 'Pur', 'Pyr')
helix_single_sub['mut_nuc_type'] = np.where((helix_single_sub['mutant_nucleotide'] == 'A') | (helix_single_sub['mutant_nucleotide'] == 'G'), 'Pur', 'Pyr')
helix_single_sub['transition_or_transversion'] = np.where(helix_single_sub['ref_nuc_type'] == helix_single_sub['mut_nuc_type'], 'Transition', 'Transversion')

for amino_acid in 'RTAVSGPL':

    helix_syn_transition_at_least_one = helix_single_sub[(helix_single_sub['transition_or_transversion'] == 'Transition') & (helix_single_sub['change_to_amino_acid'] == 'No') & (helix_single_sub['reference_amino_acid'] == amino_acid) & (helix_single_sub['counts_total'] > 0)].sort_values('counts_total',ascending=True)

    helix_syn_transversion_at_least_one = helix_single_sub[(helix_single_sub['transition_or_transversion'] == 'Transversion') & (helix_single_sub['change_to_amino_acid'] == 'No') & (helix_single_sub['reference_amino_acid'] == amino_acid) & (helix_single_sub['counts_total'] > 0)].sort_values('counts_total',ascending=True)
    
    helix_syn_transition_at_least_one.to_csv(amino_acid + '_full_DF_syn_changes_with_HelixMTdb_samples_by_AA_TRANSITION.csv')
    helix_syn_transversion_at_least_one.to_csv(amino_acid + '_full_DF_syn_changes_with_HelixMTdb_samples_by_AA_TRANSVERSION.csv')

# Syn changes with or without human samples in HelixMTdb by amino acid

syn_changes_with_without_human_samples_by_AA = []

for amino_acid in 'ACDEFGHIKLMNPQRSTVWY':
    total_syn_change = len(helix_single_sub[(helix_single_sub['reference_amino_acid'] == amino_acid) & (helix_single_sub['change_to_amino_acid'] == 'No')])
    no_human_syn_change = len(helix_single_sub[(helix_single_sub['reference_amino_acid'] == amino_acid) & (helix_single_sub['change_to_amino_acid'] == 'No') & (helix_single_sub['counts_total'] == 0)])
    yes_human_syn_change = len(helix_single_sub[(helix_single_sub['reference_amino_acid'] == amino_acid) & (helix_single_sub['change_to_amino_acid'] == 'No') & (helix_single_sub['counts_total'] > 0)])
    f_syn_change_without_human = no_human_syn_change/total_syn_change
    syn_changes_with_without_human_samples_by_AA.append([amino_acid, total_syn_change, no_human_syn_change, yes_human_syn_change, f_syn_change_without_human])

syn_changes_with_without_human_samples_by_AA_DF = pd.DataFrame(syn_changes_with_without_human_samples_by_AA,columns = ['Amino acid','Total synonymous changes','Syn changes no samples HelixMTdb','Syn changes with human samples HelixMTdb','Fraction nucleotide subs with no human samples'])
syn_changes_with_without_human_samples_by_AA_DF = syn_changes_with_without_human_samples_by_AA_DF.sort_values(by = 'Fraction nucleotide subs with no human samples',ascending=False).reset_index(drop=True)
syn_changes_with_without_human_samples_by_AA_DF.to_csv('Syn_changes_with_or_without_HelixMTdb_samples_by_AA.csv')

# How do transitions and transversions relate to counts in HelixMTdb

transitions_transversions = []

for amino_acid in 'ACDEFGHIKLMNPQRSTVWY':

    no_humans_transition_syn_change = len(helix_single_sub[(helix_single_sub['reference_amino_acid'] == amino_acid)  & (helix_single_sub['change_to_amino_acid'] == 'No') & (helix_single_sub['counts_total'] == 0) & (helix_single_sub['transition_or_transversion'] == 'Transition')])

    no_humans_transversion_syn_change = len(helix_single_sub[(helix_single_sub['reference_amino_acid'] == amino_acid) & (helix_single_sub['change_to_amino_acid'] == 'No') & (helix_single_sub['counts_total'] == 0) & (helix_single_sub['transition_or_transversion'] == 'Transversion')])

    yes_humans_transition_syn_change = len(helix_single_sub[(helix_single_sub['reference_amino_acid'] == amino_acid) & (helix_single_sub['change_to_amino_acid'] == 'No') & (helix_single_sub['counts_total'] > 0) & (helix_single_sub['transition_or_transversion'] == 'Transition')])

    yes_humans_transversion_syn_change = len(helix_single_sub[(helix_single_sub['reference_amino_acid'] == amino_acid) & (helix_single_sub['change_to_amino_acid'] == 'No') & (helix_single_sub['counts_total'] > 0) & (helix_single_sub['transition_or_transversion'] == 'Transversion')])
    
    total_syn_change_no_human = no_humans_transition_syn_change + no_humans_transversion_syn_change
    
    transitions_transversions.append([amino_acid, no_humans_transition_syn_change, no_humans_transversion_syn_change, yes_humans_transition_syn_change, yes_humans_transversion_syn_change,total_syn_change_no_human])

transitions_transversions_DF = pd.DataFrame(transitions_transversions,columns = ['Amino acid','Syn changes no human samples transition','Syn changes no human samples transversion','Syn changes yes human samples transition','Syn changes yes human samples transversion','Total synonymous changes no human samples'])

transitions_transversions_DF['f Syn changes no humans transition'] = transitions_transversions_DF['Syn changes no human samples transition']/transitions_transversions_DF['Total synonymous changes no human samples']
transitions_transversions_DF['f Syn changes no humans transversion'] = transitions_transversions_DF['Syn changes no human samples transversion']/transitions_transversions_DF['Total synonymous changes no human samples']

transitions_transversions_DF.to_csv('Syn_changes_with_or_without_HelixMTdb_samples_by_transition_transversion.csv')

