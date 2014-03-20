#!/usr/bin/env python
# encoding: utf-8
"""
parse_supplementary_table.py

Created by Rod Docking on 2014-03-19.
"""

import csv

def standardize_nomenclature(row, original_names, modified_names):
    """
    Standardize and simplify nomenclature
    """
    
    # Set up a new dict with simplified column names
    tidied_row = {}
    for old, new in zip(original_names, modified_names):
        tidied_row[new] = row.pop(old)

    # Simplify row entries
    tidied_row['RNAseq_available'] = 'TRUE' if \
        tidied_row['RNAseq_available'] == 'Yes' else 'FALSE'
    tidied_row['Expired'] = 'TRUE' if tidied_row['Expired'] == '*' else 'FALSE'
    
    return tidied_row

def main():
    """
    Parse Supplementary Table 1 from the TCGA AML paper
    """
    
    # Simplify variable names
    original_names = ['TCGA Patient ID', 'RNAseq data?', 'Expired?  4.30.13',
                      'Sex', 'Race', 'FAB', 'Age', '%BM Blast', 'WBC',
                      '%PB Blast', 'Subclones deteceted by WGS', 
                      'Cytogenetics', 'Gene Fusions by RNA-Seq',
                      'Inferred genomic rearrangement (from RNA-Seq fusion)',
                      'Cytogenetic Classification', 'RISK (Cyto)',
                      'Molecular Classification', 'RISK (Molecular)',
                      'SVs (from WGS)']
    modified_names = ['TCGA_patient_id', 'RNAseq_available', 'Expired',
                      'Sex', 'Race', 'FAB', 'Age', 'BM_blast_pct',
                      'White_blood_cell_count', 'PB_blast_pct',
                      'WGS_subclones_detected', 'Cytogenetics',
                      'RNAseq_gene_fusions', 
                      'RNAseq_inferred_genomic_rearrangement',
                      'Cytogenetic_classification', 'Cytogenetic_risk',
                      'Molecular_classification', 'Molecular_risk',
                      'SVs_from_WGS']

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
            tidied_row = standardize_nomenclature(row, original_names, 
                                                  modified_names)
            if tidied_row['TCGA_patient_id']:
                csv_writer.writerow(tidied_row)


if __name__ == '__main__':
    main()

