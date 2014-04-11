Clinical Data Import
====================

The scripts and R Markdown documents in this directory dealt with importing, inspecting, and cleaning the clinical data from the paper. These documents and scripts were largely written by R. Docking, with a substantial contribution from R. Johnston to the first document.

Import and Cleaning of Clinical Data
------------------------------------

The first R Markdown document, [import_and_clean_clinical_data.md](import_and_clean_clinical_data.md), describes the initial download, inspection, and interpretation of the clinical data from the TCGA paper.

This document was written over the course of a few days, and was updated following discussion with professors, TAs, and group members. The script also relies on a Python script, [parse_supplementary_table.py](parse_supplementary_table.py), which was used to parse the supplemental data CSV file from the TCGA paper into something more R-friendly. 

Import and Describe Count Data
------------------------------

The second R Markdown document in this sub-directory, [import_and_describe_count_data.md](import_and_describe_count_data.md), was written in response to [issue #7](https://github.com/rdocking/stat540-group-project-aml-cnv/issues/7), where it was suggested that the group use raw count data, rather than RPKM values.

We ultimately decided to stick with the RPKM data rather than the count data, but it was a useful exercise to determine just how the count values were derived.
