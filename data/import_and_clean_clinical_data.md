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

Thoughts?
