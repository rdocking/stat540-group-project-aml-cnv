STAT540 poster text
===================

## Title:
Using RNA-Seq data to predict large-scale copy number alterations in Acute Myeloid Leukemia

## Authors:
Carmen Bayly, Lauren Chong, Rod Docking, Fatemeh Dorri, Emily Hindalong, Rebecca Johnston.

## Research question:
- Hypothesis:
- Aim:

## Introduction:

Acute Myeloid Leukemia (AML) is a cancer of the blood cells characterized by recurrent, large-scale chromosomal abnormalities. Detection of these abnormalities is critical for proper treatment stratification(1, 2). The current standard of care involves testing of both cytogenetic and molecular markers, assessed through karyotyping and PCR-based tests, respectively. However, as new research identifies more clinically relevant markers, these serial tests are difficult to scale. Further, next-generation sequencing makes genome-wide testing possible, allowing for more complete testing and stratification.

The aim of this project is to infer large-scale chromosomal abnormalities in AML using RNA-seq data. RNA-seq is known to be a powerful technique for detection of SNVs, as well as AML-relevant fusion and internal tandem duplication events. 

However, detection of copy-number variants (CNVs) using RNA-seq data alone is quite challenging - copy-number gain or loss events are not reflected straightforwardly in observed RNA-seq read counts. Being able to provide accurate inference of large-scale chromosomal abnormalities from RNA-Seq data alone would be a significant advance in the field, and ultimately lead to better patient care.

## Methods:

### Data Sets

We have re-analyzed data made available through a recent publication by The Cancer Genome Atlas (3). The data-set consists of RPKM measurements made from RNA-seq libraries from 179 patients with AML, with matched clinical data noting CNV events observed through standard cytogenetics.

The specific aim of the project is to classify samples according to either their *cytogenetic risk status*, or to predict the *presence or absence of three recurrent CNV events*.

There are three levels of cytogenetic risk status according to current treatment guidelines: good, intermediate, and poor risk. These different classes correspond to the presence or absence of specific CNV events. In addition, three of the most commonly recurrent CNV events are trisomy 8 (+8), complete or partial deletion of chromosome 5 (-5/del(5q)), and complete or partial deletion of chromosome 7 (-7/del(7q)). We used the cytogenetic risk status made available by the authors, and scored each sample for presence of the three specific CNV events:

(Insert table showing breakdown of sample classes).

The data set was filtered to remove transcripts where the measured RPKM was 0 across all samples, and log2-transformed.

### Differential Expression

We used `voom` (4) to test for differential expression between:

1. Genes differentially expressed between Good, Intermediate, and Poor cytogenetic risk.
2. Genes differentially expressed in samples containing or lacking the three specific CNV events.

In both cases, we set a FDR of 1e-05, and examined the observed hits for overlap between categories, and potential biological relevance.

### Machine Learning

We applied several different machine-learning approaches to classify the samples. First, we used Principal Component Analysis (PCA) to ...

Next, we built a Support Vector Machine (SVM) model to...

Finally, we used Random Forest (RF) to...
 
## Results:


## Conclusions:


## Bibliography:

1. National Comprehensive Cancer Network. NCCN Clinical Practice Guidelines in Oncology - Acute Myeloid Leukemia. Version 2.2014. [NCCN.org](http://nccn.org) 

2. Grimwade, D. et al. Refinement of cytogenetic classification in acute myeloid leukemia: determination of prognostic significance of rare recurring chromosomal abnormalities among 5876 younger adult patients treated in the United Kingdom Medical Research Council trials. Blood 116, 354–365 (2010). [PMID: 20385793](http://www.ncbi.nlm.nih.gov/pubmed/?term=20385793)

3. The Cancer Genome Atlas Research Network. Genomic and Epigenomic Landscapes of Adult De Novo Acute Myeloid Leukemia. N. Engl. J. Med. 368, 2059–2074 (2013). [PMID: 23634996](http://www.ncbi.nlm.nih.gov/pubmed/?term=23634996)

4.	Law CW, Chen Y, Shi W, Smyth GK. Voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol. 2014 Feb 3;15(2):R29. 