Differential expression analysis for trisomy 8, deletion 5 and deletion 7
========================================================

> To knit .rmd file, read data files in using "../data"  
> To run chunks in Rstudio, read data files in using "./data"

This code performs differential expression analysis on cleaned RNA-seq read count data. In particular, it tests whether there is differential expression between different *deletions* (trisomy 8, del 5 and del 7) using `voom`.


## Load data and required libraries
Load RNA-seq data and the experimental design files:

```r
rDes <- read.table("../data/experimental_design_cleaned.txt", sep = "\t", header = TRUE, 
    row.names = 1)  #might need to fix pathname
rDat <- read.table("../data/aml.rnaseq.gaf2.0_read_count_cleaned.txt", sep = "\t", 
    header = TRUE, row.names = 1)  #might need to fix pathname
```


Load required libraries:

```r
library(lattice)  # if you don't already have this loaded ...
library(ggplot2)  # we'll make figures with both
library(reshape2)  # for the function melt
library(limma)
library(edgeR)
library(car)
library(RColorBrewer)
```


## Data inspection

```r
str(rDat, max.level = 0)
```

```
## 'data.frame':	20001 obs. of  179 variables:
##   [list output truncated]
```

```r
rDat[1:4, 1:4]
```

```
##                            X2803 X2805  X2806  X2807
## A1BG-AS|503538_calculated  794.1 431.6  893.2 1097.4
## A1BG|1_calculated         1141.2 405.4 1006.7 1123.7
## A1CF|29974_calculated        2.0   2.0    2.0    3.0
## A2LD1|87769_calculated     196.5 229.1  181.8  113.1
```

```r
head(names(rDat))
```

```
## [1] "X2803" "X2805" "X2806" "X2807" "X2808" "X2810"
```

```r
head(rownames(rDat), n = 10)
```

```
##  [1] "A1BG-AS|503538_calculated" "A1BG|1_calculated"        
##  [3] "A1CF|29974_calculated"     "A2LD1|87769_calculated"   
##  [5] "A2ML1|144568_calculated"   "A2M|2_calculated"         
##  [7] "A4GALT|53947_calculated"   "A4GNT|51146_calculated"   
##  [9] "AAA1|404744_calculated"    "AAAS|8086_calculated"
```

```r
tail(rownames(rDat), n = 10)
```

```
##  [1] "ZWINT|11130_calculated"      "ZXDA|7789_calculated"       
##  [3] "ZXDB|158586_calculated"      "ZXDC|79364_calculated"      
##  [5] "ZYG11B|79699_calculated"     "ZYX|7791_calculated"        
##  [7] "ZZEF1|23140_calculated"      "ZZZ3|26009_calculated"      
##  [9] "psiTPTE22|387590_calculated" "tAKR|389932_calculated"
```

```r
str(rDes, max.level = 0)
```

```
## 'data.frame':	179 obs. of  9 variables:
```

```r
head(rDes)
```

```
##      Sex Race FAB_subtype Age trisomy_8 del_5 del_7 Cytogenetic_risk
## 2803   F    W          M3  61     FALSE FALSE FALSE             Good
## 2805   M    W          M0  77     FALSE FALSE FALSE     Intermediate
## 2806   M    W          M1  46     FALSE FALSE FALSE             Good
## 2807   F    W          M1  68     FALSE FALSE FALSE     Intermediate
## 2808   M    W          M2  23     FALSE FALSE FALSE     Intermediate
## 2810   F    B          M2  76     FALSE FALSE FALSE             N.D.
##      Molecular_risk
## 2803           Good
## 2805   Intermediate
## 2806           Good
## 2807   Intermediate
## 2808   Intermediate
## 2810           N.D.
```


Both the design matrix and the data matrix have already been cleaned and filtered. 

RNA-seq data: there are 20001 transcripts (rows) for 179 patients (columns). Experimental design: there are 179 rows, representing information for each of the patients with RNA-seq data in the AML TCGA data set, and 179 variables.

### Differential expression analysis

I will use `voom` to perform differential expression analysis.

**Deletion**
Which genes are differentially expressed between trisomy 8, deletion 5, deletion 7?

```r
trisomy8 <- rDes$trisomy_8
table(trisomy8)
```

```
## trisomy8
## FALSE  TRUE 
##   160    19
```

```r

del5 <- rDes$del_5
table(del5)
```

```
## del5
## FALSE  TRUE 
##   163    16
```

```r

del7 <- rDes$del_7
table(del7)
```

```
## del7
## FALSE  TRUE 
##   158    21
```


Apply scale normalization:

```r
normFactor <- calcNormFactors(rDat)
```


Preparing model matrix:

```r
# Preparing Model matrices for Voom normalization-------------
modelMat_noInt <- model.matrix(~trisomy_8 * del_5 * del_7, rDes)
# The intercept represents a group with all samples without deletion of
# interest. this group changes depending on the deletion, since some samples
# share multiple deletions.
head(modelMat_noInt)
```

```
##      (Intercept) trisomy_8TRUE del_5TRUE del_7TRUE trisomy_8TRUE:del_5TRUE
## 2803           1             0         0         0                       0
## 2805           1             0         0         0                       0
## 2806           1             0         0         0                       0
## 2807           1             0         0         0                       0
## 2808           1             0         0         0                       0
## 2810           1             0         0         0                       0
##      trisomy_8TRUE:del_7TRUE del_5TRUE:del_7TRUE
## 2803                       0                   0
## 2805                       0                   0
## 2806                       0                   0
## 2807                       0                   0
## 2808                       0                   0
## 2810                       0                   0
##      trisomy_8TRUE:del_5TRUE:del_7TRUE
## 2803                                 0
## 2805                                 0
## 2806                                 0
## 2807                                 0
## 2808                                 0
## 2810                                 0
```


Now find genes differentially expressed between deletion type:

```r
dat.voomed <- voom(rDat, modelMat_noInt, plot = TRUE, lib.size = colSums(rDat))
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 

```r
fit_lm <- lmFit(dat.voomed, modelMat_noInt)
fit <- eBayes(fit_lm)
```


Investigating hits (with no interaction)

```r
ttfit_tris8 <- topTable(fit, number = Inf, coef = "trisomy_8TRUE", p.value = 1e-05)
ttfit_del5 <- topTable(fit, number = Inf, coef = "del_5TRUE", p.value = 1e-05)
ttfit_del7 <- topTable(fit, number = Inf, coef = "del_7TRUE", p.value = 1e-05)
```

there are:
6 hits for trisomy 8, 24 hits for del 5, and 141 hits for del 7.

Is there overlap between the hits?

```r
a <- length(intersect(rownames(ttfit_tris8), rownames(ttfit_del5)))
b <- length(intersect(rownames(ttfit_tris8), rownames(ttfit_del7)))
c <- length(intersect(rownames(ttfit_del5), rownames(ttfit_del7)))
sum(a, b, c)
```

```
## [1] 1
```

there is: 1 overlapping gene, at `intersect(rownames(ttfit_del5), rownames(ttfit_del7)))`

Investigating hits with interaction at 1e-5

```r
a <- nrow(ttfit_t8d5 <- topTable(fit, number = Inf, coef = "trisomy_8TRUE:del_5TRUE", 
    p.value = 1e-05))
b <- nrow(ttfit_t8d7 <- topTable(fit, number = Inf, coef = "trisomy_8TRUE:del_7TRUE", 
    p.value = 1e-05))
c <- nrow(ttfit_d5d7 <- topTable(fit, number = Inf, coef = "del_5TRUE:del_7TRUE", 
    p.value = 1e-05))
d <- nrow(ttfit_t8d5d7 <- topTable(fit, number = Inf, coef = "trisomy_8TRUE:del_5TRUE:del_7TRUE", 
    p.value = 1e-05))
sum(a, b, c, d)
```

```
## [1] 0
```

there are 0 genes with differential expression influenced by interaction effects (with cutoff of 1e-5).

Investigating hits with interaction at 1e-4

```r
a <- nrow(ttfit_t8d5 <- topTable(fit, number = Inf, coef = "trisomy_8TRUE:del_5TRUE", 
    p.value = 1e-04))
b <- nrow(ttfit_t8d7 <- topTable(fit, number = Inf, coef = "trisomy_8TRUE:del_7TRUE", 
    p.value = 1e-04))
c <- nrow(ttfit_d5d7 <- topTable(fit, number = Inf, coef = "del_5TRUE:del_7TRUE", 
    p.value = 1e-04))
d <- nrow(ttfit_t8d5d7 <- topTable(fit, number = Inf, coef = "trisomy_8TRUE:del_5TRUE:del_7TRUE", 
    p.value = 1e-04))
sum(a, b, c, d)
```

```
## [1] 0
```

there is 0 gene with differential expression influenced by interaction effects.

*Exploring hit genes from trisomy_8.*
topTable for Trisomy 8


```r
(ttfit_tris8)
```

```
##                           logFC AveExpr     t   P.Value adj.P.Val     B
## PPP2R2A|5520_calculated  0.6999   6.278 8.344 2.128e-14 4.255e-10 22.03
## NEIL2|252969_calculated  1.0329   3.084 7.861 3.809e-13 3.809e-09 18.63
## ZNF7|7553_calculated     0.6466   5.240 6.666 3.341e-10 1.758e-06 12.84
## COPS5|10987_calculated   0.5830   5.827 6.656 3.516e-10 1.758e-06 12.82
## INTS10|55174_calculated  0.6870   6.474 6.515 7.537e-10 3.015e-06 12.11
## WHSC1L1|54904_calculated 0.7291   7.255 6.328 2.034e-09 6.781e-06 11.17
```


Plotsmear of trisomy_8 hits


```r
# Create a DGEList object
voomTris8genes <- rownames(ttfit_tris8)
trisomy_8 <- rDes$trisomy_8
dgeGlmT8 <- DGEList(counts = rDat, group = as.numeric(trisomy_8))
plotSmear(dgeGlmT8, de.tags = voomTris8genes, ylab = "logFC", xlab = "AverageCounts", main = "Counts of genes differentially expressed in trisomy_8 samples")
abline(h = c(-1, 1), col = "blue")
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.png) 


Creating a boxplot with the 6 genes of interest (FDR 1e-5) for trisomy 8


```r
#subsetting and reforming transcripts of interest from main RPKM matrix
tris8Dat <- rDat[voomTris8genes,]
tris8Dat$Transcript <- rownames(tris8Dat) #creating a transcript group
tris8Dat <- melt(tris8Dat, id.vars = "Transcript", 
                   variable.name = "TCGA_patient_id",
                   value.name = "Counts")

#cleaup molten RPKM data
tris8Dat$TCGA_patient_id <- gsub("X", "", tris8Dat$TCGA_patient_id)
tris8Dat$Transcript <- gsub("[|].*$", "", tris8Dat$Transcript)

#preparing a small design matrix and merging
miniDes <- rDes[,c("trisomy_8", "del_5" , "del_7")]
miniDes$"TCGA_patient_id" <- rownames(rDes)
tris8Dat <- merge(tris8Dat, miniDes, by = "TCGA_patient_id") #merging

#plotting
ggplot(tris8Dat, aes(Transcript, Counts, colour = trisomy_8)) +
  geom_boxplot() +
  facet_wrap(~ Transcript, scales = "free") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14.png) 


*Exploring hit genes from del_5.*
TopTable for del_5


```r
(ttfit_del5)
```

```
##                                 logFC  AveExpr      t   P.Value adj.P.Val
## KIAA0087|9808_calculated       5.5886 -1.13041  8.362 1.912e-14 3.825e-10
## CCDC48|79825_calculated        4.2779  1.01960  8.172 5.974e-14 5.975e-10
## COL5A1|1289_calculated         4.2630  2.45496  7.587 1.889e-12 7.566e-09
## KIAA1257|57501_calculated      1.8314  2.56677  7.607 1.681e-12 7.566e-09
## SIGLECP3|284367_calculated     3.0916  2.88471  7.501 3.096e-12 1.032e-08
## CREB3L3|84699_calculated       4.5214 -1.48049  7.587 1.891e-12 7.566e-09
## EGF|1950_calculated            4.1183  0.21886  7.373 6.466e-12 1.847e-08
## MME|4311_calculated            3.9615  0.37304  7.089 3.237e-11 8.093e-08
## PFDN1|5201_calculated         -0.9278  5.26149 -6.971 6.259e-11 1.391e-07
## SALL4|57167_calculated         4.2145 -0.19362  6.851 1.214e-10 2.429e-07
## WDR55|54853_calculated        -1.0674  5.08447 -6.500 8.156e-10 1.483e-06
## ATP9A|10079_calculated         3.4796  2.38543  6.457 1.028e-09 1.629e-06
## CEACAM1|634_calculated         2.4744  2.85692  6.314 2.182e-09 2.728e-06
## HARS|3035_calculated          -0.8780  5.76861 -6.281 2.604e-09 3.064e-06
## PRSS2|5645_calculated          4.9984 -1.63022  6.346 1.843e-09 2.528e-06
## B3GALT5|10317_calculated       2.9013 -3.80784  6.451 1.059e-09 1.629e-06
## SYN2|6854_calculated           3.9954 -1.36087  6.341 1.896e-09 2.528e-06
## KIAA0141|9812|1of2_calculated -0.8955  6.61507 -6.198 4.014e-09 4.363e-06
## RBM22|55696_calculated        -0.6885  6.75538 -6.141 5.387e-09 5.319e-06
## ZNF793|390927_calculated       2.8784  3.10667  6.134 5.584e-09 5.319e-06
## XKR3|150165_calculated         3.7365 -1.01183  6.192 4.145e-09 4.363e-06
## STARD10|10809_calculated       1.6156  2.25264  6.089 7.061e-09 6.420e-06
## PPP2CA|5515_calculated        -0.8213  7.07684 -6.064 8.028e-09 6.981e-06
## LHX6|26468_calculated          4.1038 -0.09352  6.051 8.570e-09 7.142e-06
##                                    B
## KIAA0087|9808_calculated      21.943
## CCDC48|79825_calculated       20.959
## COL5A1|1289_calculated        17.823
## KIAA1257|57501_calculated     17.746
## SIGLECP3|284367_calculated    17.353
## CREB3L3|84699_calculated      17.085
## EGF|1950_calculated           16.384
## MME|4311_calculated           14.933
## PFDN1|5201_calculated         14.417
## SALL4|57167_calculated        13.600
## WDR55|54853_calculated        11.970
## ATP9A|10079_calculated        11.823
## CEACAM1|634_calculated        11.098
## HARS|3035_calculated          10.939
## PRSS2|5645_calculated         10.918
## B3GALT5|10317_calculated      10.777
## SYN2|6854_calculated          10.708
## KIAA0141|9812|1of2_calculated 10.536
## RBM22|55696_calculated        10.255
## ZNF793|390927_calculated      10.221
## XKR3|150165_calculated        10.042
## STARD10|10809_calculated       9.879
## PPP2CA|5515_calculated         9.874
## LHX6|26468_calculated          9.675
```


Plotsmear of del_5 hits

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.png) 


Creating a boxplot with the 24 genes of interest (FDR 1e-5) for del 5

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17.png) 


Creating a boxplot with the genes of interest (FDR 1e-6) for del 5
![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18.png) 


testing expression level cutoff: Removing genes with max counts of less than 10000
![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19.png) 


Examining impact of removing 3 samples containing apparent outliers, and then filtereing for max counts>10000. This is because going by adj.P.Val alone does seem to give false positives...

```
## X2908 
##    94
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20.png) 


PPP2CA seems to differentiate well. however, it is also the second highest row index on the topTable , so this might not be the best results filtering method.

*Exploring hit genes from del_7.*
TopTable for genes from del_7.
There are 141, so I will explore the head an tail.


```r
head(ttfit_del7)
```

```
##                            logFC AveExpr       t   P.Value adj.P.Val     B
## PDAP1|11333_calculated   -0.9866  5.5998 -11.241 2.148e-22 4.297e-18 39.91
## MKRN1|23608_calculated   -1.2127  6.9044 -10.516 2.436e-20 2.436e-16 35.50
## FAM169A|26049_calculated  2.6138  2.4384   9.917 1.160e-18 7.731e-15 31.56
## PAWR|5074_calculated      4.6390  0.2477   9.308 5.531e-17 2.766e-13 27.69
## YKT6|10652_calculated    -0.9353  5.6862  -8.991 4.028e-16 1.611e-12 26.02
## SUMF2|25870_calculated   -1.3228  6.8411  -8.804 1.281e-15 4.271e-12 24.95
```

```r
nrow(ttfit_del7)
```

```
## [1] 141
```

```r
tail(ttfit_del7)
```

```
##                               logFC AveExpr      t   P.Value adj.P.Val
## ZNF154|7710_calculated       1.8708  4.3456  5.656 6.227e-08 9.226e-06
## MYO1D|4642_calculated        2.2978  2.6546  5.643 6.640e-08 9.624e-06
## C7orf59|389541_calculated   -1.3052  5.0764 -5.639 6.787e-08 9.766e-06
## GALNT11|63917_calculated    -0.9325  5.4251 -5.634 6.960e-08 9.872e-06
## KIF7|374654|2of2_calculated  2.7038 -0.1821  5.652 6.352e-08 9.274e-06
## CDK5|1020_calculated        -1.2392  3.5607 -5.636 6.880e-08 9.829e-06
##                                 B
## ZNF154|7710_calculated      7.870
## MYO1D|4642_calculated       7.870
## C7orf59|389541_calculated   7.859
## GALNT11|63917_calculated    7.812
## KIF7|374654|2of2_calculated 7.798
## CDK5|1020_calculated        7.786
```


Plotsmear of del_7 hits
![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22.png) 


Creating a boxplot with genes of interest (FDR 1e-5) for del 7. 
Again, this is exploring the head and tail of the genelist.
![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23.png) 


![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24.png) 


Remaking topTable with a cutoff of 1e-8


```r
nrow(ttfit_del72 <- topTable(fit, number = Inf, coef = "del_7TRUE", p.value = 1e-08))
```

```
## [1] 29
```


This is a much more manageable size.

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-26.png) 


Another way to cut would be based on counts - here, excluding genes with a max count of >10000.

![plot of chunk unnamed-chunk-27](figure/unnamed-chunk-27.png) 


*Investigating hit shared by del_5 and del_7*
adj. P. Values hit shared by del_5 and del_7


```r
SharedGene <- intersect(rownames(ttfit_del5), rownames(ttfit_del7))
(ttfit_del5[SharedGene, ]$adj.P.Val)
```

```
## [1] 7.566e-09
```

```r
(ttfit_del7[SharedGene, ]$adj.P.Val)
```

```
## [1] 1.861e-08
```


Plotsmear of hit shared by `del_5` and `del_7`

![plot of chunk unnamed-chunk-29](figure/unnamed-chunk-291.png) ![plot of chunk unnamed-chunk-29](figure/unnamed-chunk-292.png) 


Boxplot of hit shared by `del_5` and `del_7`

![plot of chunk unnamed-chunk-30](figure/unnamed-chunk-30.png) 


What is this gene?
