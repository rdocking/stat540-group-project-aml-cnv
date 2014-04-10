Using RNA-Seq Data to Predict Large-scale Copy-number Alterations in AML
========================================================================

Group project for STAT540, Winter 2014.

Overview
--------
Here is a general diagram of our workflow:  

<img src="method-workflow.png" height=400>

### Project Proposal
* Proposal: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/Proposal.md)

### Inputs
* RPKM data: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/data/aml.rnaseq.gaf2.0_rpkm_cleaned.txt)
* Patient metadata: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/data/experimental_design_cleaned.txt)

### Outputs
* Poster text: **link to be added**
* Poster figures: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/tree/master/poster)

### Analysis scripts
1. Linear regression - cytogenetic risk: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/diff_expr_rna_seq_rpkm.md)
2. Linear regression - CNA events: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/Bayly_rna_seq_diff_exp_analysis.md)
3. PCA & SVM analysis:
 * Cytogenetic risk: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/pca_exploratory/pca_SVM_analysis_Cytogenic_risk.md)
 * del5: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/pca_exploratory/pca_SVM_analysis_del_5_final.Rmd)
 * del7: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/pca_exploratory/pca_SVM_analysis_del_7_final.Rmd)
 * trisomy8: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/pca_exploratory/pca_SVM_analysis_trsomy_8_final.Rmd)
4. Random forest analysis: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/rf_exploratory/rf_exploratory.md)
5. SVM analysis: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/svm_exploratory/svm_predictions_new.md)

### Bibliography
1. Genomic and Epigenomic Landscapes of Adult De Novo Acute Myeloid Leukemia. New England Journal of Medicine 368, 2059–2074 (2013). [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/23634996)
2. Grimwade, D. et al. Refinement of cytogenetic classification in acute myeloid leukemia: determination of prognostic significance of rare recurring chromosomal abnormalities among 5876 younger adult patients treated in the United Kingdom Medical Research Council trials. Blood 116, 354–365 (2010). [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/20385793)
3. Prebet, T. et al. Secondary Philadelphia chromosome after non-myeloablative peripheral blood stem cell transplantation for a myelodysplastic syndrome in transformation. Bone Marrow Transplant 33, 247–249 (2004). [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/14716291)
4. Law, C. W., Chen, Y., Shi, W. & Smyth, G. K. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology 15, R29 (2014). [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24485249)

* Missing one reference from poster
* Just the ones in the poster, or additional?