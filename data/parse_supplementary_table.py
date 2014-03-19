#!/usr/bin/env python
# encoding: utf-8
"""
parse_supplementary_table.py

Created by Rod Docking on 2014-03-19.
"""

import csv

def main():
    """
    Parse Supplementary Table 1 from the TCGA AML paper
    """
    
    supp_csv_file = "SuppTable01.update.2013.05.13.csv"
    csv_reader = csv.DictReader(open(supp_csv_file, 'r'), 
                              delimiter=',', quotechar='"')
    for row in csv_reader:
        print row['TCGA Patient ID']

if __name__ == '__main__':
    main()

