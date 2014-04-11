Using RNA-Seq Data to Predict Large-scale Copy-number Alterations in AML
========================================================================

Group project for STAT540, Winter 2014.

Overview
--------
Here is a general diagram of our workflow:  

<img src="proposal-method-workflow.png" height=400>

### Project Proposal
* Proposal: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/Proposal.md)

### Inputs
* RPKM data: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/data/aml.rnaseq.gaf2.0_rpkm_cleaned.txt)
* Patient metadata: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/data/experimental_design_cleaned.txt)

### Poster
* Poster pdf: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/poster/STAT540_AML_poster_FINAL.pdf)
* Poster text: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/poster/PosterText.md)
* Poster figures: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/tree/master/poster/images)
* Poster tables: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/tree/master/poster/tables)

### Analysis scripts
0. Import and clean data: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/tree/master/code/clinical_data_import)
1. Linear regression - Cytogenetic risk: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/diff_expr_rna_seq_rpkm.md)
2. Linear regression - CNA events: [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/Bayly_rna_seq_diff_exp_analysis.md)
3. SVM (PCA-based feature selection): [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/pca_exploratory)
4. SVM (lm and correlation-based feature selection): [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/svm_exploratory)
5. Random forest (lm and correlation-based feature selection): [here](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/code/rf_exploratory)


### Bibliography
1. National Comprehensive Cancer Network (2014). *NCCN Clinical Practice Guidelines in Oncology - Acute Myeloid Leukemia.* Version 2. [NCCN](http://www.nccn.org/default.aspx)

2. Grimwade, D et al. (2010). Renement of cytogenetic classication in acute myeloid leukemia: determination of prognostic signicance of rare recurring chromosomal abnormalities among 5876 younger adult patients treated in the United Kingdom Medical Research Council trials. *Blood*. Jul 22;116(3):354-65. [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/20385793)

3. Prebet, T et al. (2004). Secondary Philadelphia chromosome after non-myeloablative peripheral blood stem cell transplantation for a myelodysplastic syndrome in transformation. *Bone Marrow Transplantation*. Jan;33(2):247-9. [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/14716291)

4. The Cancer Genome Atlas Research Network (2013). Genomic and Epigenomic Landscapes of Adult De Novo Acute Myeloid Leukemia. *New England Journal of Medicine*. May 30;368(22):2059-74.[PubMed](http://www.ncbi.nlm.nih.gov/pubmed/23634996)

5. Law, CW et al. (2014). Voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biology*. Feb 3;15(2):R29. [PubMed](http://www.ncbi.nlm.nih.gov/pubmed/24485249)