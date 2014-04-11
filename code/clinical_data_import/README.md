Clinical Data Import
====================

The scripts and R Markdown documents in this directory dealt with importing, inspecting, and cleaning the clinical data from the paper. These documents and scripts were largely written by R. Docking and R. Johnston.

Import and Cleaning of Clinical Data
------------------------------------

The first R Markdown document, [import_and_clean_clinical_data.md](import_and_clean_clinical_data.md), describes the initial download, inspection, and interpretation of the clinical data from the TCGA paper.

This document was written over the course of a few days, and was updated following discussion with professors, TAs, and group members. The script also relies on a Python script, [parse_supplementary_table.py](parse_supplementary_table.py), which was used to parse the supplemental data CSV file from the TCGA paper into something more R-friendly. 

Import and Describe Count Data
------------------------------

The second R Markdown document in this sub-directory, [import_and_describe_count_data.md](import_and_describe_count_data.md), was written in response to [issue #7](https://github.com/rdocking/stat540-group-project-aml-cnv/issues/7), where it was suggested that the group use raw count data, rather than RPKM values.

We ultimately decided to stick with the RPKM data rather than the count data, but it was a useful exercise to determine just how the count values were derived.

Clean RNA-seq Data
------------------

The following code filters the RNA-seq RPKM and read count data and saves the output to use for all downstream analyses. 

> The filters include: 
  - Remove genes without HUGO gene IDs ("?" as gene name)
  - Change naming scheme of sample IDs to match metadata
  - Remove rows with RPKM/Read count values = 0 across all samples
  - Add 1 to all RPKM/Read counts (to avoid negative values post-log2 transformation)

- Clean RNA-seq RPKM data: [MD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/clean_rna_seq_rpkm_data.md) [RMD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/clean_rna_seq_rpkm_data.Rmd) [HTML](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/clean_rna_seq_rpkm_data.html)

- Clean RNA-seq read count data: [MD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/clean_rna_seq_read_count_data.md) [RMD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/clean_rna_seq_read_count_data.Rmd) [HTML](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/clean_rna_seq_read_count_data.html)
