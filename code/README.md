Code files
==========

This directory contains:

- Import and clean RNA-seq data and metadata:
  - Import and clean clinical data: [MD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code//clinical_data_import/import_and_clean_clinical_data.md) [RMD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/import_and_clean_clinical_data.Rmd) [HTML](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code//clinical_data_import/import_and_clean_clinical_data.html)
  - Parse clinical data: [PY](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/parse_supplementary_table.py)
  - Clean RNA-seq RPKM data: [MD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/clean_rna_seq_rpkm_data.md) [RMD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/clean_rna_seq_rpkm_data.Rmd) [HTML](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/clean_rna_seq_rpkm_data.html)
  - Description of RNA-seq read count data: [MD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/import_and_describe_count_data.md) [RMD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/import_and_describe_count_data.Rmd) [HTML](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/import_and_describe_count_data.html)
  - Clean RNA-seq read count data: [MD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/clean_rna_seq_read_count_data.md) [RMD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/clean_rna_seq_read_count_data.Rmd) [HTML](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/clinical_data_import/clean_rna_seq_read_count_data.html)
  
- Differential expression analysis:  
  - Differential Expression I - Cytogenetic Risk:
      - Using RPKM data (chosen method): [MD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/diff_expr_analysis/diff_expr_rna_seq_rpkm.md) [RMD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/diff_expr_analysis/diff_expr_rna_seq_rpkm.Rmd) [HTML](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/diff_expr_analysis/diff_expr_rna_seq_rpkm.html)
      - Using read count data: [MD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/diff_expr_analysis/diff_expr_rna_seq_read_count.md) [RMD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/diff_expr_analysis/diff_expr_rna_seq_read_count.Rmd) [HTML](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/diff_expr_analysis/diff_expr_rna_seq_read_count.html)

  - Differential Expression II - Karyotypic Events: 
      - Using RPKM data (chosen method): [MD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/diff_expr_analysis/Bayly_rna_seq_diff_exp_analysis.md) [RMD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/diff_expr_analysis/Bayly_rna_seq_diff_exp_analysis.Rmd) [HTML](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/diff_expr_analysis/Bayly_rna_seq_diff_exp_analysis.html)
      - Using read count data: [MD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/Bayly_rna_seq_diff_exp_analysis.md) [RMD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/Bayly_rna_seq_diff_exp_analysis.Rmd) [HTML](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/Bayly_rna_seq_diff_exp_analysis.html) [R](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/Bayly_rna_seq_diff_exp_analysis.R)

  - Exploratory analysis using sample correlation heatmaps:
      - Using RPKM RNA-seq data: [MD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/diff_expr_analysis/Bayly_correlation_heatmaps.md) [RMD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/diff_expr_analysis/Bayly_correlation_heatmaps.Rmd) [HTML](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/diff_expr_analysis/Bayly_correlation_heatmaps.html)
      - Using read count RNA-seq data: [MD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/rna_seq_count_data_corr_heatmap.md) [RMD](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/rna_seq_count_data_corr_heatmap.Rmd) [HTML](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/rna_seq_count_data_corr_heatmap.html)


- SVM predictions: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/svm_exploratory)

- PCA feature selection: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/pca_exploratory)

- Random forest predictions: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/rf_exploratory)
