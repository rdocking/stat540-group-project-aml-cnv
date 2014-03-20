#!/usr/bin/env python
# encoding: utf-8
"""
parse_supplementary_table.py

Created by Rod Docking on 2014-03-19.
"""

import csv

def standardize_nomenclature(row, column_names):
    """
    Standardize and simplify nomenclature
    """
    
    # Set up a new dict with simplified column names
    tidied_row = {}
    for original, modified in column_names:
        tidied_row[modified] = row.pop(original)

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
    # This is a list of tuples, with each member consisting of the 
    # original and revised variable names
    column_names = [
        ('TCGA Patient ID', 'TCGA_patient_id'),
        ('RNAseq data?', 'RNAseq_available'),
        ('Expired?  4.30.13', 'Expired'),
        ('Sex', 'Sex'), 
        ('Race', 'Race'),
        ('FAB', 'FAB_subtype'),
        ('Age', 'Age'),
        ('%BM Blast', 'BM_blast_pct'),
        ('WBC', 'White_blood_cell_count'),
        ('%PB Blast', 'PB_blast_pct'),
        ('Subclones deteceted by WGS', 'WGS_subclones_detected'),
        ('Cytogenetics', 'Cytogenetics'),
        ('Gene Fusions by RNA-Seq', 'RNAseq_gene_fusions'),
        ('Inferred genomic rearrangement (from RNA-Seq fusion)',
         'RNAseq_inferred_genomic_rearrangement'),
        ('Cytogenetic Classification', 'Cytogenetic_classification'),
        ('RISK (Cyto)', 'Cytogenetic_risk'),
        ('Molecular Classification', 'Molecular_classification'),
        ('RISK (Molecular)', 'Molecular_risk'),
        ('SVs (from WGS)', 'SVs_from_WGS')
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
            if tidied_row['TCGA_patient_id']:
                csv_writer.writerow(tidied_row)


if __name__ == '__main__':
    main()

