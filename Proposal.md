Using RNA-Seq Data to Predict Large-scale Copy-number Alterations in AML
========================================================================
  
  Group project for STAT540, Winter 2014.

Project Description
-------------------
  
Acute Myeloid Leukemia (AML) is characterized by recurrent, large-scale chromosomal abnormalities. Detection of these abnormalities is typically performed through karyotying or dedicated approaches such as array-CGH. However, these methods are technically limiting: karyotyping has low resolution, and array-based methods are limited to what is probed on the array, and thus, structural variants such as insertions or gene fusions cannot be captured. 

The aim of this project is to infer large-scale chromosomal abnormalities in AML using RNA-seq data. RNA-seq data provides an unbiased view of gene expression, and allows the detection of expressed structural variants that may be functionally relevant to the tumour phenotype. Since changes in gene expression linked to these events are not expected to be simple, we will employ a range of machine learning approaches to predict these events including principal component analysis and unsupervised clustering.

The main rationale for using RNA-Seq, rather than the more-straightforward array-based and karyotyping approaches, is to increase the clinical utility of RNA-Seq based diagnostic testing in AML. For a clinical test to be useful in guiding AML treatment, the results need to be rapidly available. Furthermore, additional tests add to the time and expense of diagnosis. The larger goal of the [AML Personalized Medicine Project](http://bccancerfoundation.com/blog/november-21-2011/personalized-medicine-project-pmp) is to provide better diagnostic tests for AML using next-generation sequencing technologies. Being able to provide accurate inference of large-scale chromosomal abnormalities from RNA-Seq data alone would be a significant advance in the field, and ultimately lead to better patient care.

Data Sets
---------
  
  We will analyse a publicly available AML dataset from The Cancer Genome Atlas (TCGA), as published in the [New England Journal of Medicine in 2013](http://www.ncbi.nlm.nih.gov/pubmed/23634996). This data set includes clinically annotated samples from a total of 200 AML patients, representing all the well-described morphologic and cytogenic subtypes of the disease. We will focus on the following available data types:
  
- RNA-seq data from 172 patients with AML, using Illumina HiSeq2000 paired-end 75 bp sequencing
- SNP-array data from the same patients, for both tumour and skin samples, using Affymetrix SNP array 6.0
- Clinical information describing karyotype results for all patients

The specific data sets to be used are available through the [TCGA Data Portal site](https://tcga-data.nci.nih.gov/docs/publications/laml_2012/) for the AML marker paper. The main data files to be used are:
  
- [RNAseq GAF 2.0 normalized RPKM](https://tcga-data.nci.nih.gov/docs/publications/laml_2012/laml.rnaseq.179_v1.0_gaf2.0_rpkm_matrix.txt.tcgaID.txt.gz) (RNA-Seq RPKM data)
- [Polymorphisms identified using the Affymetrix SNP 6 platform](https://tcga-data.nci.nih.gov/docs/publications/laml_2012/LAML.Genome_Wide_SNP_6.Level_3.tgz) (SNP-array data)
- [Patient Clinical Data](https://tcga-data.nci.nih.gov/docs/publications/laml_2012/clinical_patient_laml.tsv)

Some example `R` code for importing and inspecting the main data files can be found at [data_import.R](https://github.com/rdocking/stat540-group-project-aml-cnv/blob/master/data_import.R). Further cleaning and interpretation will be required, but the data straightforwardly loads into R, and the main expression data has already been normalized into an analysis-ready format.

Initial Analysis Plan
---------------------
  
### Data Summary and Tidying
  
- Summarize the available data files
- Construct a simplified categorical variable grouping the patients into categories based on large-scale chromosomal abnormalities

### Data exploration and differential expression analysis

- Perform sample correlation (`heatmap()` function from `lattice`)
- Principal Component Analysis (PCA) to infer sub-groups present within RNA-seq data
- Find most variable components in the data using colour for different factors: Cytogenetic abnormality? Race? Sex? Age? etc? 	
- Differential expression analysis of RNA-seq data (`voom` from `R` package `limma`)
- Use results from PCA to inform model
- Hierarchical clustering of samples given differentially expressed genes (`heatmap()` function from `lattice`)

### Correlation between CNAs and RNA expression

- Do the predicted CNAs ([supplementary table S5](https://tcga-data.nci.nih.gov/docs/publications/laml_2012/SupplementalTable05.tsv) from the TCGA paper) lead to concurrent changes in gene expression?
- Individual sample scale: Circos plots? Manual checks of key chromosomal losses and gains?
- Across all samples: How can we analyse the correlation in an automated fashion? Map CNAs to genes, then call up and down-regulation and gain and loss for each gene predicted to have CNA?

### Machine learning

The idea for applying machine learning methods can be divided in two steps:
  
#### 1- Pre-processing
  
In this step, we aim to prepare the data for the main analysis by first doing sanity checks and the applying appropriate normalization for removing the systematic variations.

#### 2- Data Analysis

- Clustering (Unsupervised): We can use Independent Component Analysis (ICA) for extracting the biological significant dimensions from RNA-seq data. ICA assumes non-Gaussian expression variation and models the Micorarray observations as linear combination of its component. The components are chosen to be as independent as possible. 
- Recreate Figure 4A: "Unsupervised RNA expression patterns" with RNA abundance heatmaps and sample annotations, including AML FAB subtype

- Classification (supervised): We can apply Linear discriminant analysis or SVM  for classification of the data.
We can also do the main analysis on data with smaller number of features, that is basically reducing the dimensionality of the data in order to identify the salient features.  The dimensionality reduction step can be accomplished by using PCA. Afterwards, we will be able to compare the result of the "original data" and "data with smaller number of features" and conclude how it is necessary in RNA-seq.

Group Composition
-----------------

| Name  | Program | Supervisor  | Expertise |
| ------------- | ------------- | ----- | ------- |
| Fatemeh Dorri | Computer Science (PhD) | Dr. Sohrab Shah and Dr. Anne Condon | Machine learning, Data mining, Bioinformatics | 
| Rod Docking | Experimental Medicine (PhD) | Dr. Aly Karsan | Bioinformatics, Genomics |
| Lauren Chong | Bioinformatics (MSc) | Rotations (current: Dr. Ryan Morin) | Bioinformatics, R |
| Rebecca Johnston | Bioinformatics (MSc) | Rotations (current: Dr. Christian Steidl) | Molecular biology, cancer biology |
| Emily Hindalong | Bioinformatics (MSc) | Current Rotation: Dr. Paul Pavlidis | Computer Science, Software, Cognitive Systems |
| Carmen Bayly | GSAT | Current rotation: Dr. Corey Nislow | Biochemistry |

Deliverables and Responsibilities
---------------------------------

The main deliverables for the project will be a poster, as well as this git repository. Within this git repository, the work will be broken down into a series of [R Markdown](http://rmarkdown.rstudio.com) documents.  

Per the group meeting on 2014-04-02, here is the breakdown of the relevant R Markdown documents, their interconnection, and the group member in charge of each document:

|  Document Number   |   Subject   |   Owner   | Input | Output |
| ------------------ | ----------- | --------- | ----- | ------ | 
| **1**    | Experimental Design Sheet Generation     | Rod  | Raw Clinical Data | Cleaned Experimental design CSV |
| **2**    | Expression Matrix Generation | Rebecca | Raw count and RPKM Data | Cleaned count data | 
| **3**    | Differential Expression I - Cytogenetic Risk | Rebecca | Output of **1**, **2** | Table of differentially expressed genes, with associated q-values, logFC |
| **4**    | Differential Expression II - Karyotypic Events | Carmen | Output of **1**, **2** | Table of differentially expressed genes, with associated q-values, logFC |
| **5** | Differential Expression Summary | CB, RD, RJ | Output of **3**, **4** | Summary of the differential expression analysis |
| **6** | Principal Components Analysis | Fatemeh | Output of **1**, **2** | Principal Components |
| **7** | Support Vector Machine  | Emily | Output of **1, 2, 6** | SVM |
| **8** | Random Forest  | Lauren | Output of **1, 2, 6** | RF |
| **9** | ML Summary | FD, EH, LC | Output of **6, 7, 8** | Performance summary (sensitivity/sepcifity) for classifiers |

Note that the 'owner' in each case should not be the sole person working on each document, but will be the main group member responsible for its content. All steps should be modular enough so that, for example, a change in the expression matrix data (2) can be easily incorporated in subsequent steps.

### 1. Experimental Design Sheet Generation

The task here is to turn the raw clinical data (from the Supplemental material made available in the publication), into an analysis-ready CSV file describing the experimental design of the study.

The output CSV file will be used in all subsequent analysis steps.

This task can be viewed at [import_and_clean_clinical_data.md](code/import_and_clean_clinical_data.md)

### 2. Expression Matrix Generation

From the raw RPKM and count data available, summarize and clean the data into an analysis-ready format.

The output matrix data will be used in all subsequent analysis steps.

### 3. Differential Expression I - Cytogenetic Risk

Using `voom`, predict differentially expressed genes between the 'Good', 'Poor', and 'Intermediate' cytogenetic risk statuses.

### 4. Differential Expression II - Karyotpyic Events

Using `voom`, predict differentially expressed genes between samples with or without the karyotypic events of interest.

### 5. Differential Expression Summary

Compare and summarize the approaches described in **3, 4**.

### 6. Principal Components Analysis

Using the matrix of expression data, predict principal components.

### 7. SVM Analysis

Using the matrix of expression data (with perhaps a reduced data set produced in **6**), construct and test a SVM.

### 8. Random Forest Analysis

Using the matrix of expression data (with perhaps a reduced data set produced in **6**), construct and test a Random Forest model.

### 8. Machine Learning Summary

Summarize the work from documents **6, 7, 8**.
