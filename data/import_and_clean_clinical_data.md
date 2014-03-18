Clinical Data Import and Cleaning
=================================

Read in the clinical data-sheet from the [TCGA Publication Website](https://tcga-data.nci.nih.gov/docs/publications/laml_2012/):


```r
library(RCurl)
```

```
## Loading required package: bitops
```

```r
library(knitr)
d <- getURL("https://tcga-data.nci.nih.gov/docs/publications/laml_2012/clinical_patient_laml.tsv")
raw_clinical_data <- read.table(text = d, header = TRUE, sep = "\t")
```


This data-sheet has a lot more data than we'll actually need for this project:


```r
dim(raw_clinical_data)
```

```
## [1] 200  78
```


Also available, as an Excel sheet, is a Supplementary table listing much of the same information, but with more annotation:

*Note that we're using the updated version as of 2013-05-13. Since the data is only available as an Excel sheet, we've done some manipulation outside of R:*

- Download the Supplementary file from [the paper website](https://tcga-data.nci.nih.gov/docs/publications/laml_2012/SuppTable01.update.2013.05.13.xlsx)
- Open in Excel, Save as CSV


```r
supp_d <- read.csv("SuppTable01.update.2013.05.13.csv")
dim(supp_d)
```

```
## [1]  201 2006
```


As a first pass, we'll try to clean the `cytogenetic_abnormality` column from the raw clinical data into a simpler categorical variable:


```r
summary(raw_clinical_data$cytogenetic_abnormality)
```

```
##                                                                                     [Not Available] 
##                                                                                                  19 
##                                         Complex (greater than or equal to 3 distinct abnormalities) 
##                                                                                                   4 
##                          Complex (greater than or equal to 3 distinct abnormalities)|del (5q) / 5q- 
##                                                                                                   4 
##           Complex (greater than or equal to 3 distinct abnormalities)|del (5q) / 5q-|del (7q) / 7q- 
##                                                                                                   4 
## Complex (greater than or equal to 3 distinct abnormalities)|del (5q) / 5q-|del (7q) / 7q-|Trisomy 8 
##                                                                                                   2 
##                Complex (greater than or equal to 3 distinct abnormalities)|del (5q) / 5q-|Trisomy 8 
##                                                                                                   2 
##                Complex (greater than or equal to 3 distinct abnormalities)|del (7q) / 7q-|Trisomy 8 
##                                                                                                   1 
##                               Complex (greater than or equal to 3 distinct abnormalities)|t (15;17) 
##                                                                                                   1 
##                               Complex (greater than or equal to 3 distinct abnormalities)|Trisomy 8 
##                                                                                                   2 
##                                                                                      del (5q) / 5q- 
##                                                                                                   1 
##                                                                       del (5q) / 5q-|del (7q) / 7q- 
##                                                                                                   1 
##                                                                                      del (7q) / 7q- 
##                                                                                                   3 
##                                                                                            inv (16) 
##                                                                                                   4 
##                                                                                              Normal 
##                                                                                                 102 
##                                  Normal|Complex (greater than or equal to 3 distinct abnormalities) 
##                                                                                                   1 
##                   Normal|Complex (greater than or equal to 3 distinct abnormalities)|del (7q) / 7q- 
##                                                                                                   2 
##         Normal|Complex (greater than or equal to 3 distinct abnormalities)|del (7q) / 7q-|Trisomy 8 
##                                                                                                   1 
##                        Normal|Complex (greater than or equal to 3 distinct abnormalities)|Trisomy 8 
##                                                                                                   1 
##                                                                Normal|del (5q) / 5q-|del (7q) / 7q- 
##                                                                                                   1 
##                                                                               Normal|del (7q) / 7q- 
##                                                                                                   4 
##                                                                      Normal|del (7q) / 7q-|t (9;11) 
##                                                                                                   1 
##                                                           Normal|del (7q) / 7q-|Trisomy 8|t (15;17) 
##                                                                                                   1 
##                                                                                     Normal|inv (16) 
##                                                                                                   5 
##                                                                                    Normal|t (15;17) 
##                                                                                                   5 
##                                                                                     Normal|t (8;21) 
##                                                                                                   1 
##                                                                                     Normal|t (9;11) 
##                                                                                                   1 
##                                                                                    Normal|Trisomy 8 
##                                                                                                   4 
##                                                                                           t (15;17) 
##                                                                                                   9 
##                                                                                            t (8;21) 
##                                                                                                   6 
##                                                                                           Trisomy 8 
##                                                                                                   6 
##                                                                                 Trisomy 8|t (15;17) 
##                                                                                                   1
```


There are some immediate problems with this approach:

- 19 samples are missing data
- Several samples have many of the CNAs of interest
- The notation seems inconsistent between samples

*Note: Supplemenary Table 5 from the paper website might be an alternate approach here - this lists the inferred CNAs directly*

Alternate approaches possible:

1. Limit the analysis to things with a single, distinct karyotypic event (will lose lots of samples)
2. Make a 'best guess' at the most relevant alteration for each sample (messy)
3. Make several different binary variables for each of the main CNAs, e.g.:

| Patient | del(5q) | del(7q) | Trisomy 8 | t(15:17) | inv(16) | 
| ------- | ------- | ------- | --------- | -------- | ------- |
| 1 | T | F | F | F |T |
| 2 | F | F | F |T | F |

Thoughts? This would make the annotation simpler, but might complicate the analysis a bit.

### Followup to that question: 

#### Summarize Columns in Raw Clinical Data and Supplementary Table 1

To address the comments at [issue #3](https://github.com/rdocking/stat540-group-project-aml-cnv/issues/3), I'm going to try a few approaches to stratifying the data-set.

First, I'll work on trying to get a workable stratification out of the `raw_clinical_data` imported above. The main issues I've been having so far:

- In the `cytogenetic_abnormality` column, the entries are sometimes split by pipes, e.g. 'Normal|del (7q) / 7q-' vs. 'del (7q) / 7q-' - it's not clear that these should be treated equivalently.
- There is a second column, `cytogenetic_abnormality_other` that is even more inconsistently filled in (i.e. no data for 172 samples, uses 'no', 'No', and 'NO' for others)

This is looking a bit like a rabbit-hole. I'm going to try another tack and see if I can regenerate the sample counts from Table 1 in the main text of the paper:

![Table 1 Snippet](table_1_snippet.png)

From Supplementary Table 1, (loaded above as `supp_d`), the cytogenetic data is available in the following columns:

- `Cytogenetics` - karyotype using [cytogenetic nomenclature](http://www.radford.edu/~rsheehy/cytogenetics/Cytogenetic_Nomeclature.html). These seem consistently filled in (yay!), but are tricky to parse (boo!)
- `Gene Fusions by RNA-Seq` - predicted gene fusions from the RNA-Seq data. This has entries like 'PML(+)RARA(+) (In frame)', which corresponds to 't(15;17)' from the `Cytogenetics` column. This _should_ be pretty accurate for fusions, but not for CNAs
- `Inferred genomic rearrangement (from RNA-Seq fusion)` - translation of `Gene Fusions by RNA-Seq` into `Cytogenetics` nomenclature
- `Cytogenetic Classification` - A categorical variable summarizing the previous columns:


```r
summary(supp_d$Cytogenetic.Classification)
```

```
##                                           
##                                         1 
##                                  BCR-ABL1 
##                                         3 
##                                CBFB-MYH11 
##                                        12 
##                      Complex Cytogenetics 
##                                        24 
## Intermediate Risk Cytogenetic Abnormality 
##                                        22 
##              MLL translocation, poor risk 
##                                         5 
##                MLL translocation, t(9;11) 
##                                         2 
##                                      N.D. 
##                                         5 
##                          Normal Karyotype 
##                                        92 
##                                  PML-RARA 
##                                        18 
##         Poor Risk Cytogenetic Abnormality 
##                                        10 
##                             RUNX1-RUNX1T1 
##                                         7
```


  this is _not quite_ what we're after - we'd like to call individual karyotype-level events.
- `RISK (Cyto)` - another categorical classification of the samples by risk type:


```r
summary(supp_d$RISK..Cyto.)
```

```
##                      Good Intermediate         N.D.         Poor 
##            1           37          115            5           43
```


- Then two more classifications, `Molecular Classification` and `RISK (Molecular)`
  
#### Annotate Samples with Karyotypic events of interest  
  
