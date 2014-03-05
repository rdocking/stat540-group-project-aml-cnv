Using RNA-Seq Data to Predict Large-scale Copy-number Alterations in AML
========================================================================

Group project for STAT540, Winter 2014.

Project Description
-------------------

Acute Myeloid Leukemia (AML) is characterized by recurrent, large-scale chromosomal abnormalities. Detection of these abnormalities is typically performed through karyotying or dedicated approaches such as array-CGH. However, these methods are technically limiting: karyotyping has low resolution, and array-based methods are limited to what is probed on the array, and thus, structural variants such as insertions or gene fusions cannot be captured.  

The aim of this project is to infer large-scale chromosomal abnormalities in AML using RNA-seq data. RNA-seq data provides an unbiased view of gene expression, and allows the detection of expressed structural variants that may be functionally relevant to the tumour phenotype. Since changes in gene expression linked to these events are not expected to be simple, we will employ a range of approaches to predict these events, including:

- unsupervised and supervised clustering of RNA-seq data
- machine learning

Data Sets
---------

We will analyse a publicly available AML dataset from The Cancer Genome Atlas (TCGA), as published in the [New England Journal of Medicine in 2013](http://www.ncbi.nlm.nih.gov/pubmed/23634996). This data set includes clinically annotated samples from a total of 200 AML patients, representing all the well-described morphologic and cytogenic subtypes of the disease. We will focus on the following available data types:

- RNA-seq data from 172 patients with AML, using Illumina HiSeq2000 paired-end 75 bp sequencing
- SNP-array data from the same patients, for both tumour and skin samples, using Affymetrix SNP array 6.0

The goal is to use the second data set as 'ground truth' by which to validate methods developed for the first data set.

Group Composition
-----------------

| Name  | Program | Supervisor  | Expertise |
| ------------- | ------------- | ----- | ------- |
| Fatemeh Dorri | Computer Science (PhD) | Dr. Sohrab Shah and Dr. Anne Condon | Machine learning, Data mining, Bioinformatics | 
| Rod Docking | Experimental Medicine (PhD) | Dr. Aly Karsan | Bioinformatics, Genomics |
| Lauren Chong | Bioinformatics (MSc) | Rotations (current: Dr. Ryan Morin) | Bioinformatics, R |
| Rebecca Johnston | Bioinformatics (MSc) | Rotations (current: Dr. Christian Steidl) | Molecular biology, cancer biology |
| Emily Hindalong | Bioinformatics (MSc) | Current Rotation: Dr. Paul Pavlidis | Computer Science, Software, Cognitive Systems |
| Carmen Bayly | GSAT | | Biochemistry |
