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

Per conversation with group members and Shaun, I'll clean the supplemental data in `SuppTable01.update.2013.05.13.csv`. Note that I'm doing this in Python (with ` parse_supplementary_table.py`) for expediency.

Manipulations are as follows:

- Select a relevant subset of variables
- For variables that can be converted to boolean types, do so
- Convert column names to be more R-friendly (e.g., `'Expired?  4.30.13'` became `Expired` and `'%BM Blast'` became `BM_blast_pct`)

The script should be run like so:

```
python parse_supplementary_table.py > experimental_design.csv
```


Now I can read in the modified CSV file to get a cleaner data frame:


```r
cleaned_data <- read.csv("experimental_design.csv")
```


This data frame is quite a bit smaller and simpler than the original:


```r
str(cleaned_data)
```

```
## 'data.frame':	200 obs. of  19 variables:
##  $ TCGA_patient_id                      : int  2803 2806 2870 2815 2872 2998 2914 2819 2875 2823 ...
##  $ RNAseq_available                     : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
##  $ Expired                              : logi  TRUE TRUE TRUE TRUE FALSE TRUE ...
##  $ Sex                                  : Factor w/ 2 levels "F","M": 1 2 2 2 2 1 1 1 2 1 ...
##  $ Race                                 : Factor w/ 13 levels "A","B","H","NH/A",..: 12 12 12 12 2 12 12 12 7 12 ...
##  $ FAB                                  : Factor w/ 9 levels "M0","M1","M2",..: 4 2 2 5 4 4 3 3 3 4 ...
##  $ Age                                  : int  61 46 76 49 42 68 22 52 43 61 ...
##  $ BM_blast_pct                         : int  44 90 73 81 88 85 55 67 40 73 ...
##  $ White_blood_cell_count               : num  1 29.4 34 57.1 2.1 29 51.8 4.1 4.3 86.4 ...
##  $ PB_blast_pct                         : int  NA 81 55 48 2 32 70 18 39 68 ...
##  $ WGS_subclones_detected               : int  NA NA NA NA NA 1 NA NA NA NA ...
##  $ Cytogenetics                         : Factor w/ 125 levels "37~49,XY,+Y,der(1)add(1)(p13)del(1)(q21q25),-5,der(7)inv(7)(p15q11.2)?inv(7)(q22q32),+17,add(17)(p13),+21,+mar[cp20]",..: 52 13 74 72 69 51 46 53 15 100 ...
##  $ RNAseq_gene_fusions                  : Factor w/ 59 levels " ","AF086125(+)USP22(+) (Out of Frame),CA7(-)USP22(+) (Out of Frame)",..: 49 56 35 39 44 50 37 54 54 46 ...
##  $ RNAseq_inferred_genomic_rearrangement: Factor w/ 52 levels " ","del17q11.2",..: 35 49 40 39 35 35 42 48 48 35 ...
##  $ Cytogenetic_classification           : Factor w/ 11 levels "BCR-ABL1","CBFB-MYH11",..: 9 11 2 2 9 9 2 11 11 9 ...
##  $ Cytogenetic_risk                     : Factor w/ 4 levels "Good","Intermediate",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ Molecular_classification             : Factor w/ 13 levels "BCR-ABL1","CBFB-MYH11",..: 11 13 2 2 11 11 2 13 13 11 ...
##  $ Molecular_risk                       : Factor w/ 4 levels "Good","Intermediate",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ SVs_from_WGS                         : Factor w/ 19 levels " ","1:166453112-3:53150925(CTX),12:128688979-12:128689118(DEL),14:84366862-14:84371916(DEL),15:32718038-15:32718374(DEL),15:6683662"| __truncated__,..: 1 1 1 1 1 7 1 1 1 1 ...
```



  
