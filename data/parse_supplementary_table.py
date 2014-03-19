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
    
    for old, new in zip(original_names, modified_names):
        row[new] = row.pop(old)

    # Simplify row entries
    row['RNAseq_available'] = 'TRUE' if \
        row['RNAseq_available'] == 'Yes' else 'FALSE'
    row['Expired'] = 'TRUE' if row['Expired'] == '*' else 'FALSE'
    
    return row

def main():
    """
    Parse Supplementary Table 1 from the TCGA AML paper
    """
    
    # Open the CSV file for reading
    supp_csv_file = "SuppTable01.update.2013.05.13.csv"
    csv_reader = csv.DictReader(open(supp_csv_file, 'r'), 
                              delimiter=',', quotechar='"')

    # Simplify variable names
    original_names = ['TCGA Patient ID', 'RNAseq data?', 'Expired?  4.30.13',
                      'Sex', ]
    modified_names = ['TCGA_patient_id', 'RNAseq_available', 'Expired',
                      'Sex']

    # Print CSV header
    print ','.join(modified_names)
    # Parse the CSV file row-by-row
    for row in csv_reader:
        # Clean up nomenclature
        row = standardize_nomenclature(row, original_names, modified_names)
        # Skip if no TCGA ID
        if not row['TCGA_patient_id']:
            continue
        print ",".join(row[f] for f in modified_names)

if __name__ == '__main__':
    main()

