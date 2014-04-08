Differential expression analysis for trisomy 8, deletion 5 and deletion 7
========================================================

> This code is for plotting hits of interest from the DEA of the 3key deletion genotypes (trisomy 8, del 5 and del 7).
> It performs differential expression analysis on cleaned RNA-seq read count data, and when given abrreviated hit names, 
> it finds the full names, prepares the appropriate data matrix for plotting and then boxplots them.

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
library(ggplot2)  # we'll make figures with both
library(reshape2)  # for the function melt
library(limma)
library(edgeR)
# library(RColorBrewer)
```


Both the design matrix and the data matrix have already been cleaned, filtered and inspected.

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


Finding genes differentially expressed between deletion type:

```r
dat.voomed <- voom(rDat, modelMat_noInt, plot = TRUE, lib.size = colSums(rDat) * 
    normFactor)
```

![plot of chunk Vooming and Limming and Baysing - oh my!](figure/Vooming_and_Limming_and_Baysing_-_oh_my_.png) 

```r
fit_lm <- lmFit(dat.voomed, modelMat_noInt)
fit <- eBayes(fit_lm)
```


Making topTables for the hits (with no interaction)

```r
ttfit_tris8 <- topTable(fit, number = Inf, coef = "trisomy_8TRUE", p.value = 1e-05)
ttfit_del5 <- topTable(fit, number = Inf, coef = "del_5TRUE", p.value = 1e-05)
ttfit_del7 <- topTable(fit, number = Inf, coef = "del_7TRUE", p.value = 1e-05)
```

there are:
4 hits for trisomy 8, 29 hits for del 5, and 190 hits for del 7.

*Retrieving boxplots for the 4 hits of trisomy 8*

```r
(ttfit_tris8)
```

```
##                           logFC AveExpr     t   P.Value adj.P.Val     B
## NEIL2|252969_calculated  0.9572   3.084 7.761 6.879e-13 1.376e-08 18.05
## PPP2R2A|5520_calculated  0.6299   6.278 7.323 8.630e-12 8.630e-08 16.31
## ZNF7|7553_calculated     0.5785   5.240 6.698 2.813e-10 1.875e-06 12.99
## WHSC1L1|54904_calculated 0.6583   7.255 6.578 5.388e-10 2.694e-06 12.42
```

```r

voomTris8genes <- rownames(ttfit_tris8)
tris8Dat <- rDat[voomTris8genes, ]
tris8Dat$Transcript <- rownames(tris8Dat)  #creating a transcript group
tris8Dat <- melt(tris8Dat, id.vars = "Transcript", variable.name = "TCGA_patient_id", 
    value.name = "Counts")

# cleaup molten counts data
tris8Dat$TCGA_patient_id <- gsub("X", "", tris8Dat$TCGA_patient_id)
tris8Dat$Transcript <- gsub("[|].*$", "", tris8Dat$Transcript)

# preparing a small design matrix and merging
miniDes <- rDes[, c("trisomy_8", "del_5", "del_7")]
miniDes$TCGA_patient_id <- rownames(rDes)
tris8Dat <- merge(tris8Dat, miniDes, by = "TCGA_patient_id")  #merging
```


Boxplot for trisomy 8 genes

```r
# plotting
ggplot(tris8Dat, aes(Transcript, log2(Counts), colour = trisomy_8)) + geom_boxplot() + 
    facet_wrap(~Transcript, scales = "free") + theme(axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), aspect.ratio = 2) + ggtitle("Hits for Trisomy 8")
```

![plot of chunk Trisomy 8: Boxplots](figure/Trisomy_8:_Boxplots.png) 


*Retrieving hit genes from del_5.*

Getting full names of genes for plotting (note they all use caps).

```r
(ttfit_del5)
```

```
##                                 logFC  AveExpr      t   P.Value adj.P.Val
## KIAA0087|9808_calculated       5.6456 -1.13041  8.577 5.195e-15 1.039e-10
## CCDC48|79825_calculated        4.1990  1.01960  8.035 1.367e-13 1.367e-09
## KIAA1257|57501_calculated      1.7930  2.56677  7.710 9.259e-13 4.630e-09
## CREB3L3|84699_calculated       4.5521 -1.48049  7.724 8.510e-13 4.630e-09
## COL5A1|1289_calculated         4.1928  2.45496  7.526 2.697e-12 1.079e-08
## EGF|1950_calculated            4.0692  0.21886  7.421 4.929e-12 1.643e-08
## SIGLECP3|284367_calculated     3.0554  2.88471  7.278 1.119e-11 2.796e-08
## MME|4311_calculated            4.0055  0.37304  7.286 1.065e-11 2.796e-08
## SALL4|57167_calculated         4.2178 -0.19362  6.968 6.380e-11 1.418e-07
## KDM3B|51780_calculated        -0.8282  7.80337 -6.908 8.884e-11 1.777e-07
## HARS|3035_calculated          -0.9414  5.76861 -6.803 1.584e-10 2.881e-07
## ATP9A|10079_calculated         3.4181  2.38543  6.453 1.050e-09 1.750e-06
## CEACAM1|634_calculated         2.3850  2.85692  6.369 1.643e-09 2.191e-06
## PRSS2|5645_calculated          4.9834 -1.63022  6.389 1.472e-09 2.102e-06
## KIAA0141|9812|1of2_calculated -0.9614  6.61507 -6.265 2.833e-09 3.148e-06
## PFDN1|5201_calculated         -0.9939  5.26149 -6.269 2.775e-09 3.148e-06
## RBM22|55696_calculated        -0.7448  6.75538 -6.235 3.318e-09 3.493e-06
## B3GALT5|10317_calculated       2.8209 -3.80784  6.406 1.351e-09 2.079e-06
## XKR3|150165_calculated         3.7480 -1.01183  6.271 2.741e-09 3.148e-06
## TCOF1|6949_calculated         -0.9526  6.68845 -6.145 5.295e-09 5.295e-06
## ZNF793|390927_calculated       2.7844  3.10667  6.106 6.461e-09 5.619e-06
## LHX6|26468_calculated          4.0936 -0.09352  6.112 6.267e-09 5.619e-06
## WDR55|54853_calculated        -1.1310  5.08447 -6.049 8.682e-09 6.946e-06
## SYN2|6854_calculated           3.9667 -1.36087  6.107 6.425e-09 5.619e-06
## STARD10|10809_calculated       1.5538  2.25264  6.004 1.094e-08 8.416e-06
## ELN|2006_calculated            3.8871  0.19153  5.983 1.215e-08 9.004e-06
## NFASC|23114_calculated         3.2288 -1.08132  6.062 8.135e-09 6.779e-06
## VPREB1|7441_calculated         4.3300 -0.87257  5.975 1.266e-08 9.043e-06
## CLNK|116449_calculated         3.6920 -0.31977  5.968 1.315e-08 9.069e-06
##                                    B
## KIAA0087|9808_calculated      23.210
## CCDC48|79825_calculated       20.209
## KIAA1257|57501_calculated     18.337
## CREB3L3|84699_calculated      17.879
## COL5A1|1289_calculated        17.491
## EGF|1950_calculated           16.673
## SIGLECP3|284367_calculated    16.136
## MME|4311_calculated           16.000
## SALL4|57167_calculated        14.229
## KDM3B|51780_calculated        14.176
## HARS|3035_calculated          13.598
## ATP9A|10079_calculated        11.807
## CEACAM1|634_calculated        11.372
## PRSS2|5645_calculated         11.164
## KIAA0141|9812|1of2_calculated 10.868
## PFDN1|5201_calculated         10.862
## RBM22|55696_calculated        10.716
## B3GALT5|10317_calculated      10.665
## XKR3|150165_calculated        10.460
## TCOF1|6949_calculated         10.272
## ZNF793|390927_calculated      10.084
## LHX6|26468_calculated          9.986
## WDR55|54853_calculated         9.772
## SYN2|6854_calculated           9.676
## STARD10|10809_calculated       9.495
## ELN|2006_calculated            9.364
## NFASC|23114_calculated         9.349
## VPREB1|7441_calculated         9.283
## CLNK|116449_calculated         9.209
```

```r
# higher
a <- grep("COL5A1", rownames(ttfit_del5))
b <- grep("CEACAM1", rownames(ttfit_del5))
c <- grep("CREB3L3", rownames(ttfit_del5))

# lower
d <- grep("RBM22", rownames(ttfit_del5))
e <- grep("KDM3B", rownames(ttfit_del5))
f <- grep("HARS", rownames(ttfit_del5))

# alltogether
del5_index <- c(a, b, c, d, e, f)

# list of rownames for plugging back into the system
voomDel5genes <- rownames(ttfit_del5[del5_index, ])
```


This is just the plumbing for graph-making.

```r
del5Dat <- rDat[voomDel5genes, ]  #subsetting transcripts of interest from rDat
del5Dat$Transcript <- rownames(del5Dat)  #creating a transcript group
del5Dat <- melt(del5Dat, id.vars = "Transcript", variable.name = "TCGA_patient_id", 
    value.name = "Counts")

# cleaup molten counts data
del5Dat$TCGA_patient_id <- gsub("X", "", del5Dat$TCGA_patient_id)  #removing the 'X' preceding patient ID names
del5Dat$Transcript <- gsub("[|].*$", "", del5Dat$Transcript)  #shortening transcript names

# miniDes already made while exploring trisomy 8.
del5Dat <- merge(del5Dat, miniDes, by = "TCGA_patient_id")  #merging
```


This is for the actual plot.

```r
ggplot(del5Dat, aes(Transcript, log2(Counts), colour = del_5)) + geom_boxplot() + 
    facet_wrap(~Transcript, scales = "free") + theme(axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), aspect.ratio = 2) + ggtitle("Hits for Del 5")
```

![plot of chunk Del 5: Boxplots](figure/Del_5:_Boxplots.png) 


*Exploring hit genes from del_7.*
TopTable for genes from del_7.
There are 190, so I will explore the head an tail.


```r
head(ttfit_del7)
```

```
##                            logFC AveExpr       t   P.Value adj.P.Val     B
## PDAP1|11333_calculated   -1.0883  5.5998 -11.688 1.149e-23 2.299e-19 42.75
## LUC7L2|51631_calculated  -1.0413  7.3370 -10.015 6.244e-19 6.245e-15 32.37
## MKRN1|23608_calculated   -1.3148  6.9044  -9.899 1.311e-18 8.739e-15 31.64
## FAM169A|26049_calculated  2.4923  2.4384   9.728 3.903e-18 1.952e-14 30.41
## C7orf42|55069_calculated -0.8228  7.0560  -9.535 1.330e-17 5.321e-14 29.40
## PAWR|5074_calculated      4.5181  0.2477   9.182 1.228e-16 4.095e-13 26.96
```

```r
nrow(ttfit_del7)
```

```
## [1] 190
```

```r

# higher
a <- grep("CDCP1", rownames(ttfit_del7))
b <- grep("FAM169A", rownames(ttfit_del7))
c <- grep("BEND4", rownames(ttfit_del7))

# lower
d <- grep("LUC7L2", rownames(ttfit_del7))
e <- grep("MKRN1", rownames(ttfit_del7))
f <- grep("SLC25A13", rownames(ttfit_del7))

# alltogether
(del7_index <- c(a, b, c, d, e, f))
```

```
## [1] 19  4 29  2  3  8
```

```r

# list of rownames for plugging back into the system
voomDel7genes <- rownames(ttfit_del7[del7_index, ])
```


This is just the plumbing for graph-making.

```r
del7Dat <- rDat[voomDel7genes, ]  #subsetting transcripts of interest from rDat
del7Dat$Transcript <- rownames(del7Dat)  #creating a transcript group
del7Dat <- melt(del7Dat, id.vars = "Transcript", variable.name = "TCGA_patient_id", 
    value.name = "Counts")

# cleaup molten counts data
del7Dat$TCGA_patient_id <- gsub("X", "", del7Dat$TCGA_patient_id)  #removing the 'X' preceding patient ID names
del7Dat$Transcript <- gsub("[|].*$", "", del7Dat$Transcript)  #shortening transcript names

# miniDes already made while exploring trisomy 8.
del7Dat <- merge(del7Dat, miniDes, by = "TCGA_patient_id")  #merging
```



```r
ggplot(del7Dat, aes(Transcript, log2(Counts), colour = del_7)) + geom_boxplot() + 
    facet_wrap(~Transcript, scales = "free") + theme(axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), aspect.ratio = 2) + ggtitle("Hits for Del 7")
```

![plot of chunk Del 7: Boxplots](figure/Del_7:_Boxplots.png) 


*Investigating hit shared by del_5 and del_7*
adj. P. Values hit shared by del_5 and del_7


```r
# Retrieving full hit name
(SharedGene <- intersect(rownames(ttfit_del5), rownames(ttfit_del7)))
```

```
## [1] "COL5A1|1289_calculated"
```

```r

# preparing hit for boxplot
sgDat <- rDat[SharedGene, ]  #subsetting transcripts of interest from rDat
sgDat$Transcript <- rownames(sgDat)  #creating a transcript group
sgDat <- melt(sgDat, id.vars = "Transcript", variable.name = "TCGA_patient_id", 
    value.name = "Counts")

# cleaup molten counts data
sgDat$TCGA_patient_id <- gsub("X", "", sgDat$TCGA_patient_id)  #removing the 'X' preceding patient ID names
sgDat$Transcript <- gsub("[|].*$", "", sgDat$Transcript)  #shortening transcript names

# miniDes already made while exploring trisomy 8.
sgDat <- merge(sgDat, miniDes, by = "TCGA_patient_id")  #merging
```



```r
ggplot(sgDat, aes(Transcript, log2(Counts), colour = del_5)) + geom_boxplot() + 
    facet_wrap(~Transcript, scales = "free") + theme(axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), aspect.ratio = 2) + ggtitle("Collagen, type V, alpha 1: del 5 grouping")
```

![plot of chunk Del 5 & Del 7: Boxplots using both groupings](figure/Del_5___Del_7:_Boxplots_using_both_groupings1.png) 

```r

ggplot(sgDat, aes(Transcript, log2(Counts), colour = del_7)) + geom_boxplot() + 
    facet_wrap(~Transcript, scales = "free") + theme(axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), aspect.ratio = 2) + ggtitle("Collagen, type V, alpha 1: del 7 grouping")
```

![plot of chunk Del 5 & Del 7: Boxplots using both groupings](figure/Del_5___Del_7:_Boxplots_using_both_groupings2.png) 


