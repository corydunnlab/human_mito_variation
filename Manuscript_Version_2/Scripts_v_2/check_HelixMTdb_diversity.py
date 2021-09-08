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

# Load appropriate CSV

helixMTdb_single_subs = pd.read_csv('HelixMTdb_results_NOT_grouped_by_AA_change_single_nuc_sub_from_human_ref.csv',index_col=0)

# Obtain a list of synonymous substitutions at amino acids with four-fold degenerate third codon positions

four_fold_AA_mito = ['A','G','P','R','T','V']
helixMTdb_AGPRTV = helixMTdb_single_subs[helixMTdb_single_subs['reference_amino_acid'].str.contains('|'.join(four_fold_AA_mito))]
helixMTdb_four_fold = helixMTdb_AGPRTV[helixMTdb_AGPRTV['change_to_amino_acid'] == 'No']
list_of_four_fold_AA_mito = list(set(helixMTdb_four_fold['Amino_acid_change'].tolist()))

# How many synonymous substitutions of those possible (amino acids with four-fold degenerate third codon positions) were seen in the HelixMTdb? L and S avoided.

list_of_sites_and_destinations_from_reference_four_fold = []

for site in list_of_four_fold_AA_mito:
    count_of_destination_bases = len(helixMTdb_four_fold[(helixMTdb_four_fold.Amino_acid_change == site) & (helixMTdb_four_fold.counts_total > 0)])
    to_append = site,count_of_destination_bases
    list_of_sites_and_destinations_from_reference_four_fold.append(to_append)
    
sorted_list_of_sites_and_destinations_from_reference_four_fold = sorted(list_of_sites_and_destinations_from_reference_four_fold, key=lambda k: [k[1]])
destinations_for_four_fold_sites = pd.DataFrame(sorted_list_of_sites_and_destinations_from_reference_four_fold, \
                                                columns = ['Amino_acid_site', 'Destinations_from_reference'])
destinations_for_four_fold_sites.to_csv('destinations_for_four_fold_sites.csv')

# Obtain a list of synonymous substitutions at amino acids with two-fold degenerate third codon positions

two_fold_AA_mito = ['M','Q','E','F','N','K','D','W','I','Y','C','H']
helixMTdb_MQEFNKDWIYCH = helixMTdb_single_subs[helixMTdb_single_subs['reference_amino_acid'].str.contains('|'.join(two_fold_AA_mito))]
helixMTdb_two_fold = helixMTdb_MQEFNKDWIYCH[helixMTdb_MQEFNKDWIYCH['change_to_amino_acid'] == 'No']
list_of_two_fold_AA_mito = list(set(helixMTdb_two_fold['Amino_acid_change'].tolist()))

# How many synonymous substitutions of those possible (amino acids with two-fold degenerate third codon positions) were seen in the HelixMTdb?

list_of_sites_and_destinations_from_reference_two_fold = []

for site in list_of_two_fold_AA_mito:
    count_of_destination_bases = len(helixMTdb_two_fold[(helixMTdb_two_fold.Amino_acid_change == site) & (helixMTdb_two_fold.counts_total > 0)])
    to_append = site,count_of_destination_bases
    list_of_sites_and_destinations_from_reference_two_fold.append(to_append)
    
sorted_list_of_sites_and_destinations_from_reference_two_fold = sorted(list_of_sites_and_destinations_from_reference_two_fold, key=lambda k: [k[1]])
destinations_for_two_fold_sites = pd.DataFrame(sorted_list_of_sites_and_destinations_from_reference_two_fold, \
                                                columns = ['Amino_acid_site', 'Destinations_from_reference'])
destinations_for_two_fold_sites.to_csv('destinations_for_two_fold_sites.csv')

