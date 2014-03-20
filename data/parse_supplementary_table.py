#!/usr/bin/env python
# encoding: utf-8
"""
parse_supplementary_table.py

Created by Rod Docking on 2014-03-19.
"""

import csv

def tidy_field(field, revised, var_type):
    """
    Tidy field values in a standardized way
    """
    
    # For factorial variables, keep all non-whitespace entries
    if var_type == 'factor':
        field = field if field.strip() else 'NA'
    # For boolean types, set TRUE if text is entered, otherwise FALSE
    elif var_type == 'bool':
        field = 'TRUE' if field.strip() else 'FALSE'
    elif var_type == 'custom':
        if revised == 'RNAseq_available':
            field = 'TRUE' if field == 'Yes' else 'FALSE'
        
    return field

def standardize_nomenclature(row, column_names):
    """
    Standardize and simplify nomenclature
    """
    
    # Set up a new dict with simplified column names
    tidied_row = {}
    for original, revised, var_type in column_names:
        tidied_row[revised] = row.pop(original)
        tidied_row[revised] = tidy_field(tidied_row[revised], 
                                         revised, var_type)
    return tidied_row

def main():
    """
    Parse Supplementary Table 1 from the TCGA AML paper
    """
    
    # Simplify variable names
    # This is a list of tuples, with each member consisting of the 
    # original and revised variable names, and the 
    # type: either 'bool', 'factor', or 'custom' 
    # (for fields that don't fit the other types)
    # NOTE: 'factor' here isn't meant to imply that these variables
    #  will be treated as factorial variables in R
    column_names = [
        ('TCGA Patient ID', 'TCGA_patient_id', 'factor'),
        ('RNAseq data?', 'RNAseq_available', 'custom'),
        ('Expired?  4.30.13', 'Expired', 'bool'),
        ('Sex', 'Sex', 'factor'), 
        ('Race', 'Race', 'factor'),
        ('FAB', 'FAB_subtype', 'factor'),
        ('Age', 'Age', 'factor'),
        ('%BM Blast', 'BM_blast_pct', 'factor'),
        ('WBC', 'White_blood_cell_count', 'factor'),
        ('%PB Blast', 'PB_blast_pct', 'factor'),
        ('Subclones deteceted by WGS', 'WGS_subclones_detected', 'factor'),
        ('Cytogenetics', 'Cytogenetics', 'factor'),
        ('Gene Fusions by RNA-Seq', 'RNAseq_gene_fusions', 'factor'),
        ('Inferred genomic rearrangement (from RNA-Seq fusion)',
         'RNAseq_inferred_genomic_rearrangement', 'factor'),
        ('Cytogenetic Classification', 'Cytogenetic_classification', 'factor'),
        ('RISK (Cyto)', 'Cytogenetic_risk', 'factor'),
        ('Molecular Classification', 'Molecular_classification', 'factor'),
        ('RISK (Molecular)', 'Molecular_risk', 'factor'),
        ('SVs (from WGS)', 'SVs_from_WGS', 'factor'),
        ('EFS months    4.30.13', 'event_free_survival_months', 'factor'),
        ('OS months  4.30.13', 'overall_survival_months','factor'),
        ('PML-RARA', 'PML_RARA_present', 'bool'),
        ('MLL-partner', 'MLL_partner', 'factor'),
        ('MYH11-CBFB', 'MYH11_CBFB', 'bool'),
        ('RUNX1-RUNX1T1', 'RUNX1_RUNX1T1', 'bool'),
        ('MLLT10-partner', 'MLLT10_partner', 'factor'),
        ('BCR-ABL', 'BCR_ABL', 'bool'),
        ('GPR128-TFG', 'GPR128_TFG', 'bool'),
        ('NUP98-NSD1', 'NUP98_NSD1', 'bool'),
        ('MLL-PTD', 'MLL_PTD', 'factor'),
        ('Other in -frame fusions', 'other_in_frame_fusions', 'factor'),
        ('FLT3', 'FLT3', 'factor'),
        ('NPM1', 'NPM1', 'factor'),
        ('DNMT3A', 'DNMT3A', 'factor')
    ]
    modified_names = [t[1] for t in column_names]

    # Open the CSV file for reading
    supp_csv_file = "SuppTable01.update.2013.05.13.csv"
    csv_reader = csv.DictReader(open(supp_csv_file, 'r'), 
                            delimiter=',', quotechar='"')

    # Parse the CSV file row-by-row and write output
    with open('experimental_design.csv','wb') as out_handle:
        csv_writer = csv.DictWriter(out_handle, delimiter=',', 
                            fieldnames=modified_names)
        csv_writer.writeheader()
        for row in csv_reader:
            tidied_row = standardize_nomenclature(row, column_names)
            if not tidied_row['TCGA_patient_id'] == 'NA':
                csv_writer.writerow(tidied_row)


if __name__ == '__main__':
    main()

