stat540-group-project-aml-cnv
=============================

Group project for STAT540, Winter 2014. Detection of CNV events in Acute Myeloid Leukemia using RNA-seq data

Project Description
-------------------

Acute Myeloid Leukemia (AML) is characterized by recurrent, large-scale chromosomal abnormalities.
Detection of these abnormalities is typically performed through karyoptying or dedicated approaches such as array-CGH. However these methods are technically limiting, as structural variants cannot be captured by array-based methods (is that correct? May also want to include a statement about why WGS/WES is not sufficient here?).  

This project will (take a novel approach to?) infer large-scale chromosomal abnormalities in AML using RNA-seq data. RNA-seq data provides an unbiased view of gene expression, and allows the detection of structural variants that may be functionally relevant to the tumour phenotype. Since changes in gene expression linked to these events are not expected to be simple, we will employ a range of approaches (TBD, but most likely including some sort of machine learning or clustering approaches) to predict these events. 

Data Sets
---------

We will analyse a publicly available AML dataset from The Cancer Genome Atlas (TCGA), as published in the [New England Journal of Medicine in 2013](http://www.ncbi.nlm.nih.gov/pubmed/23634996). This data set includes clinically annotated samples from a total of 200 AML patients. We will focus on the following available data types:

- RNA-seq data from 172 patients with AML
- SNP-array data from the same patients
- Karyotype information detailing large-scale chromosomal abnormalities (?)

The goal is to use the second two data sets as 'ground truth' by which to validate methods developed for the first data set.

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
