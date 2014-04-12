Boxplotter for Differential expression analysis for trisomy 8, deletion 5 and deletion 7
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
# normFactor <- calcNormFactors(rDat)
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

```
## Error: object 'normFactor' not found
```

```r
fit_lm <- lmFit(dat.voomed, modelMat_noInt)
```

```
## Error: object 'dat.voomed' not found
```

```r
fit <- eBayes(fit_lm)
```

```
## Error: object 'fit_lm' not found
```


Making topTables for the hits (with no interaction)

```r
ttfit_tris8 <- topTable(fit, number = Inf, coef = "trisomy_8TRUE", p.value = 1e-05)
```

```
## Error: object 'fit' not found
```

```r
ttfit_del5 <- topTable(fit, number = Inf, coef = "del_5TRUE", p.value = 1e-05)
```

```
## Error: object 'fit' not found
```

```r
ttfit_del7 <- topTable(fit, number = Inf, coef = "del_7TRUE", p.value = 1e-05)
```

```
## Error: object 'fit' not found
```

there are:


```

Error in nrow(ttfit_tris8) : object 'ttfit_tris8' not found

```

 hits for trisomy 8, 

```

Error in nrow(ttfit_del5) : object 'ttfit_del5' not found

```

 hits for del 5, and 

```

Error in nrow(ttfit_del7) : object 'ttfit_del7' not found

```

 hits for del 7.

*Retrieving boxplots for the 4 hits of trisomy 8*

```r
(ttfit_tris8)
```

```
## Error: object 'ttfit_tris8' not found
```

```r

voomTris8genes <- rownames(ttfit_tris8)
```

```
## Error: object 'ttfit_tris8' not found
```

```r
tris8Dat <- rDat[voomTris8genes, ]
```

```
## Error: object 'voomTris8genes' not found
```

```r
tris8Dat$Transcript <- rownames(tris8Dat)  #creating a transcript group
```

```
## Error: object 'tris8Dat' not found
```

```r
tris8Dat <- melt(tris8Dat, id.vars = "Transcript", variable.name = "TCGA_patient_id", 
    value.name = "RPKM")
```

```
## Error: object 'tris8Dat' not found
```

```r

# cleaup molten RPKM data
tris8Dat$TCGA_patient_id <- gsub("X", "", tris8Dat$TCGA_patient_id)
```

```
## Error: object 'tris8Dat' not found
```

```r
tris8Dat$Transcript <- gsub("[|].*$", "", tris8Dat$Transcript)
```

```
## Error: object 'tris8Dat' not found
```

```r

# preparing a small design matrix and merging
miniDes <- rDes[, c("trisomy_8", "del_5", "del_7")]
miniDes$TCGA_patient_id <- rownames(rDes)
tris8Dat <- merge(tris8Dat, miniDes, by = "TCGA_patient_id")  #merging
```

```
## Error: object 'tris8Dat' not found
```


Boxplot for trisomy 8 genes

```r
# plotting
ggplot(tris8Dat, aes(Transcript, log2(RPKM), colour = trisomy_8)) + geom_boxplot() + 
    facet_wrap(~Transcript, scales = "free") + theme(axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), aspect.ratio = 2) + ggtitle("Hits for Trisomy 8")
```

```
## Error: object 'tris8Dat' not found
```


*Retrieving hit genes from del_5.*

Getting full names of genes for plotting (note they all use caps).

```r
(ttfit_del5)
```

```
## Error: object 'ttfit_del5' not found
```

```r
# higher
a <- grep("COL5A1", rownames(ttfit_del5))
```

```
## Error: object 'ttfit_del5' not found
```

```r
b <- grep("CEACAM1", rownames(ttfit_del5))
```

```
## Error: object 'ttfit_del5' not found
```

```r
c <- grep("CREB3L3", rownames(ttfit_del5))
```

```
## Error: object 'ttfit_del5' not found
```

```r

# lower
d <- grep("RBM22", rownames(ttfit_del5))
```

```
## Error: object 'ttfit_del5' not found
```

```r
e <- grep("KDM3B", rownames(ttfit_del5))
```

```
## Error: object 'ttfit_del5' not found
```

```r
f <- grep("HARS", rownames(ttfit_del5))
```

```
## Error: object 'ttfit_del5' not found
```

```r

# alltogether
(del5_index <- c(a, b, c, d, e, f))
```

```
## Error: object 'a' not found
```

```r

# list of rownames for plugging back into the system
voomDel5genes <- rownames(ttfit_del5[del5_index, ])
```

```
## Error: object 'ttfit_del5' not found
```


This is just the plumbing for graph-making.

```r
del5Dat <- rDat[voomDel5genes, ]  #subsetting transcripts of interest from rDat
```

```
## Error: object 'voomDel5genes' not found
```

```r
del5Dat$Transcript <- rownames(del5Dat)  #creating a transcript group
```

```
## Error: object 'del5Dat' not found
```

```r
del5Dat <- melt(del5Dat, id.vars = "Transcript", variable.name = "TCGA_patient_id", 
    value.name = "RPKM")
```

```
## Error: object 'del5Dat' not found
```

```r

# cleaup molten RPKM data
del5Dat$TCGA_patient_id <- gsub("X", "", del5Dat$TCGA_patient_id)  #removing the 'X' preceding patient ID names
```

```
## Error: object 'del5Dat' not found
```

```r
del5Dat$Transcript <- gsub("[|].*$", "", del5Dat$Transcript)  #shortening transcript names
```

```
## Error: object 'del5Dat' not found
```

```r

# miniDes already made while exploring trisomy 8.
del5Dat <- merge(del5Dat, miniDes, by = "TCGA_patient_id")  #merging
```

```
## Error: object 'del5Dat' not found
```


This is for the actual plot.

```r
ggplot(del5Dat, aes(Transcript, log2(RPKM), colour = del_5)) + geom_boxplot() + 
    facet_wrap(~Transcript, scales = "free") + theme(axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), aspect.ratio = 2) + ggtitle("Hits for Del 5")
```

```
## Error: object 'del5Dat' not found
```


*Exploring hit genes from del_7.*
TopTable for genes from del_7.
There are 

```

Error in nrow(ttfit_del7) : object 'ttfit_del7' not found

```

, so I will explore the head an tail.


```r
head(ttfit_del7)
```

```
## Error: object 'ttfit_del7' not found
```

```r
nrow(ttfit_del7)
```

```
## Error: object 'ttfit_del7' not found
```

```r

# higher
a <- grep("CDCP1", rownames(ttfit_del7))
```

```
## Error: object 'ttfit_del7' not found
```

```r
b <- grep("FAM169A", rownames(ttfit_del7))
```

```
## Error: object 'ttfit_del7' not found
```

```r
c <- grep("BEND4", rownames(ttfit_del7))
```

```
## Error: object 'ttfit_del7' not found
```

```r

# lower
d <- grep("LUC7L2", rownames(ttfit_del7))
```

```
## Error: object 'ttfit_del7' not found
```

```r
e <- grep("MKRN1", rownames(ttfit_del7))
```

```
## Error: object 'ttfit_del7' not found
```

```r
f <- grep("SLC25A13", rownames(ttfit_del7))
```

```
## Error: object 'ttfit_del7' not found
```

```r

# alltogether
(del7_index <- c(a, b, c, d, e, f))
```

```
## Error: object 'a' not found
```

```r

# list of rownames for plugging back into the system
voomDel7genes <- rownames(ttfit_del7[del7_index, ])
```

```
## Error: object 'ttfit_del7' not found
```


This is just the plumbing for graph-making.

```r
del7Dat <- rDat[voomDel7genes, ]  #subsetting transcripts of interest from rDat
```

```
## Error: object 'voomDel7genes' not found
```

```r
del7Dat$Transcript <- rownames(del7Dat)  #creating a transcript group
```

```
## Error: object 'del7Dat' not found
```

```r
del7Dat <- melt(del7Dat, id.vars = "Transcript", variable.name = "TCGA_patient_id", 
    value.name = "RPKM")
```

```
## Error: object 'del7Dat' not found
```

```r

# cleaup molten RPKM data
del7Dat$TCGA_patient_id <- gsub("X", "", del7Dat$TCGA_patient_id)  #removing the 'X' preceding patient ID names
```

```
## Error: object 'del7Dat' not found
```

```r
del7Dat$Transcript <- gsub("[|].*$", "", del7Dat$Transcript)  #shortening transcript names
```

```
## Error: object 'del7Dat' not found
```

```r

# miniDes already made while exploring trisomy 8.
del7Dat <- merge(del7Dat, miniDes, by = "TCGA_patient_id")  #merging
```

```
## Error: object 'del7Dat' not found
```



```r
ggplot(del7Dat, aes(Transcript, log2(RPKM), colour = del_7)) + geom_boxplot() + 
    facet_wrap(~Transcript, scales = "free") + theme(axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), aspect.ratio = 2) + ggtitle("Hits for Del 7")
```

```
## Error: object 'del7Dat' not found
```


*Investigating hit shared by del_5 and del_7*
adj. P. Values hit shared by del_5 and del_7


```r
# Retrieving full hit name
(SharedGene <- intersect(rownames(ttfit_del5), rownames(ttfit_del7)))
```

```
## Error: object 'ttfit_del7' not found
```

```r

# preparing hit for boxplot
sgDat <- rDat[SharedGene, ]  #subsetting transcripts of interest from rDat
```

```
## Error: object 'SharedGene' not found
```

```r
sgDat$Transcript <- rownames(sgDat)  #creating a transcript group
```

```
## Error: object 'sgDat' not found
```

```r
sgDat <- melt(sgDat, id.vars = "Transcript", variable.name = "TCGA_patient_id", 
    value.name = "RPKM")
```

```
## Error: object 'sgDat' not found
```

```r

# cleaup molten RPKM data
sgDat$TCGA_patient_id <- gsub("X", "", sgDat$TCGA_patient_id)  #removing the 'X' preceding patient ID names
```

```
## Error: object 'sgDat' not found
```

```r
sgDat$Transcript <- gsub("[|].*$", "", sgDat$Transcript)  #shortening transcript names
```

```
## Error: object 'sgDat' not found
```

```r

# miniDes already made while exploring trisomy 8.
sgDat <- merge(sgDat, miniDes, by = "TCGA_patient_id")  #merging
```

```
## Error: object 'sgDat' not found
```



```r
ggplot(sgDat, aes(Transcript, log2(RPKM), colour = del_5)) + geom_boxplot() + 
    facet_wrap(~Transcript, scales = "free") + theme(axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), aspect.ratio = 2) + ggtitle("Collagen, type V, alpha 1: del 5 grouping")
```

```
## Error: object 'sgDat' not found
```

```r

ggplot(sgDat, aes(Transcript, log2(RPKM), colour = del_7)) + geom_boxplot() + 
    facet_wrap(~Transcript, scales = "free") + theme(axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), aspect.ratio = 2) + ggtitle("Collagen, type V, alpha 1: del 7 grouping")
```

```
## Error: object 'sgDat' not found
```


