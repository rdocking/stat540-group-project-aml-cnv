Data files
==========

This directory contains: 
- [RNA-seq data](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/data/laml.rnaseq.179_v1.0_gaf2.0_rpkm_matrix.txt.tcgaID.txt) GAF 2.0 normalized RPKM values
- [Experimental design](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/data/experimental_design.csv) The  covariates associated with each of the 200 samples, the output of our [script](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/parse_supplementary_table.py)
- [Supplementary Table 1](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/data/SuppTable01.update.2013.05.13.csv) This is the raw clinical data, the input for our [script](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/parse_supplementary_table.py)
- [Clinical data](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/data/clinical_patient_laml.tsv) A back up file, since our code uses `wget` to access the file from the TCGA website. We were originally going to use this data set for our experimental design, but found supplementary table 1 from the paper had more annotations.

The generation and cleaning of the experimental design file is described in [import_and_clean_clinical_data.md](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/import_and_clean_clinical_data.Rmd)
