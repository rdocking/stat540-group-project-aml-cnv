RNA-seq differential expression analysis
========================================
> To knit .rmd file, read data files in using "../data"  
> To run chunks in Rstudio, read data files in using "./data"

Use `voom + limma` to perform differential expression analysis on the RNA-seq data.

Load required libraries:

```r
library(RColorBrewer)
library(reshape2)
library(plyr)
library(ggplot2)
library(limma)
library(edgeR)
```


Load RNA-seq data and the experimental design:

```r
rDat <- read.table("../data/laml.rnaseq.179_v1.0_gaf2.0_rpkm_matrix.txt.tcgaID.txt", 
    sep = "\t", header = TRUE, row.names = 1)

rDes <- read.csv("../data/experimental_design.csv")
```


Data inspection:

```r
str(rDat, max.level = 0)
```

```
## 'data.frame':	20442 obs. of  179 variables:
##   [list output truncated]
```

```r
rDat[1:4, 1:4]
```

```
##                        TCGA.AB.2803 TCGA.AB.2807 TCGA.AB.2963 TCGA.AB.2826
## ?|100132510_calculated       1.1250       0.3295       1.9477        2.735
## ?|100134860_calculated      13.0095      15.0336       5.7934        4.625
## ?|10357_calculated           0.2028       0.0831       0.2774        1.127
## ?|10431_calculated          36.1781      17.8495      34.9591       29.438
```

```r
head(names(rDat))
```

```
## [1] "TCGA.AB.2803" "TCGA.AB.2807" "TCGA.AB.2963" "TCGA.AB.2826"
## [5] "TCGA.AB.2867" "TCGA.AB.2818"
```

```r
head(rownames(rDat), n = 10)
```

```
##  [1] "?|100132510_calculated" "?|100134860_calculated"
##  [3] "?|10357_calculated"     "?|10431_calculated"    
##  [5] "?|114130_calculated"    "?|115669_calculated"   
##  [7] "?|120126_calculated"    "?|1231_calculated"     
##  [9] "?|127550_calculated"    "?|136157_calculated"
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
## 'data.frame':	200 obs. of  48 variables:
```




RNA-seq data: there are 20442 transcripts (rows) for 179 patients (columns). The row names are strange, I will attempt to change them.

Experimental design: there are 200 rows, representing information for each of the patients in the AML TCGA data set, and 179 variables, but I will only choose some variables for differential expression analysis. Note the sample ID naming scheme does not match across `rDat` and `rDes`, so I will need to fix this too.

### Clean rDat
Change sample names in `rDat` to match `rDes` by extracting substrings, i.e. extract the numbers in each sample name:

```r
names(rDat) <- regmatches(names(rDat), regexpr("(?<=AB.).*", names(rDat), perl = TRUE))
head(names(rDat))
```

```
## [1] "2803" "2807" "2963" "2826" "2867" "2818"
```


Now to fix the row names: some begin with "?" so I assume this means they cannot be associated with a gene. Also, each row name has the suffix "|", followed by an integer, then "_calculated". I will remove the rows with "?" in the row name first:

```r
# '?' only present at the start of the row name: rownames(rDat)[grep('[?]',
# rownames(rDat))]
length(grep("[?]", rownames(rDat)))
```

```
## [1] 123
```

```r
rDat <- rDat[grep("[?]", rownames(rDat), invert = TRUE), ]
dim(rDat)
```

```
## [1] 20319   179
```


In the row names, I attempted to remove substrings after and including the "|" symbol, but this was not allowed since I obtain non-unique row names. What does the integer in suffix of the row names mean?:

```r
# test <- tail(rownames(rDat)) gsub('[|].*$', '', test) rownames(rDat) <-
# gsub('[|].*$', '', rownames(rDat))
```



**Filtering:** 
Remove transcripts with RPKM = 0 across all samples:

```r
# Number of transcripts with RPKM = 0 for all samples
nrow(rDat[rowSums(rDat) == 0, ])
```

```
## [1] 318
```

```r
# Remove these transcripts
rDat <- rDat[rowSums(rDat) != 0, ]
dim(rDat)
```

```
## [1] 20001   179
```

This filter does not remove all RPKM values of 0. What do we do when RPKM = 0? Do we need to apply more filters for low RPKM values? Or do we add a small number to all RPKM values?

Potential additional filters:

```r
nrow(rDat)
```

```
## [1] 20001
```

```r
# 1. Remove rows where sum RPKM values across all samples < 5
nrow(rDat) - nrow(rDat[rowSums(rDat) < 5, ])
```

```
## [1] 17441
```

```r
# 2. Remove rows where at least one sample has RPKM value = 0; too stringent
nrow(rDat) - nrow(rDat[apply(rDat, 1, prod) != 0, ])
```

```
## [1] 7097
```

```r
# 3. Remove rows where more than 50 samples have RPKM values < 1; too
# stringent
nrow(rDat) - nrow(rDat[apply(rDat, 1, function(x) sum(abs(x) < 1) < 50), ])
```

```
## [1] 8690
```

I have actually run the differential expression analysis code right through using the 3rd filter (remove rows where more than 50 samples have RPKM values < 1), but I found almost no genes had logFC > 1, and only 8 genes were called as differentially expressed with FDR 1e-5, and they did not appear to be sex related. So I believe the 3rd filter is far too stringent.

Instead, I have chosen to apply the 1st filter, since it does not remove > 50% of the probes:

```r
# Remove transcripts where sum of RPKM values across all samples is < 5
rDat <- rDat[rowSums(rDat) > 5, ]
dim(rDat)
```

```
## [1] 17441   179
```

```r
head(rDat[1:5, 1:5])
```

```
##                             2803    2807   2963   2826   2867
## A1BG-AS|503538_calculated 7.3159 10.3701 2.9940 3.5681 5.3165
## A1BG|1_calculated         7.4715  7.5407 1.6134 1.9149 2.8357
## A2LD1|87769_calculated    1.7760  1.0395 1.5691 0.9025 1.4967
## A2ML1|144568_calculated   0.1024  0.0478 0.1885 0.0539 0.0829
## A2M|2_calculated          5.4296  5.0108 0.6181 0.6500 2.5318
```






### Clean rDes
I don't need all the variables stored in the experimental design file, so I will only keep the columns I want to test using a linear model for differential expression:

```r
rDes <- rDes[, c("TCGA_patient_id", "Sex", "Race", "FAB_subtype", "Age", "trisomy_8", 
    "del_5", "del_7", "Cytogenetic_risk", "Molecular_risk")]
str(rDes)
```

```
## 'data.frame':	200 obs. of  10 variables:
##  $ TCGA_patient_id : int  2803 2806 2870 2815 2872 2998 2914 2819 2875 2823 ...
##  $ Sex             : Factor w/ 2 levels "F","M": 1 2 2 2 2 1 1 1 2 1 ...
##  $ Race            : Factor w/ 13 levels "A","B","H","NH/A",..: 12 12 12 12 2 12 12 12 7 12 ...
##  $ FAB_subtype     : Factor w/ 9 levels "M0","M1","M2",..: 4 2 2 5 4 4 3 3 3 4 ...
##  $ Age             : int  61 46 76 49 42 68 22 52 43 61 ...
##  $ trisomy_8       : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
##  $ del_5           : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
##  $ del_7           : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
##  $ Cytogenetic_risk: Factor w/ 4 levels "Good","Intermediate",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ Molecular_risk  : Factor w/ 4 levels "Good","Intermediate",..: 1 1 1 1 1 1 1 1 1 1 ...
```

Much better, now I only have 10 variables (columns).

Now I need to subset the 200 patients in `rDes` by the 179 patients in `rDat`:

```r
rDes <- rDes[rDes$TCGA_patient_id %in% names(rDat), ]
dim(rDes)
```

```
## [1] 179  10
```


Next, I need to ensure `rDat` and `rDes` are in the same order:

```r
rDat <- rDat[, order(names(rDat))]
head(rDat[1:5, 1:5])
```

```
##                             2803   2805   2806    2807   2808
## A1BG-AS|503538_calculated 7.3159 3.1563 5.8948 10.3701 5.0019
## A1BG|1_calculated         7.4715 2.1048 4.7194  7.5407 3.3064
## A2LD1|87769_calculated    1.7760 1.6495 1.1761  1.0395 1.0664
## A2ML1|144568_calculated   0.1024 0.1126 0.1380  0.0478 0.0864
## A2M|2_calculated          5.4296 0.8499 0.5144  5.0108 3.9734
```

```r
rDes <- rDes[order(rDes$TCGA_patient_id), ]
rDes$TCGA_patient_id <- as.character(rDes$TCGA_patient_id)
head(rDes)
```

```
##     TCGA_patient_id Sex Race FAB_subtype Age trisomy_8 del_5 del_7
## 1              2803   F    W          M3  61     FALSE FALSE FALSE
## 40             2805   M    W          M0  77     FALSE FALSE FALSE
## 2              2806   M    W          M1  46     FALSE FALSE FALSE
## 41             2807   F    W          M1  68     FALSE FALSE FALSE
## 42             2808   M    W          M2  23     FALSE FALSE FALSE
## 153            2810   F    B          M2  76     FALSE FALSE FALSE
##     Cytogenetic_risk Molecular_risk
## 1               Good           Good
## 40      Intermediate   Intermediate
## 2               Good           Good
## 41      Intermediate   Intermediate
## 42      Intermediate   Intermediate
## 153             N.D.           N.D.
```

```r
identical(names(rDat), rDes$TCGA_patient_id)
```

```
## [1] TRUE
```



### Explore rDat
Check the RNA-seq data using a sample-to-sample correlation heatmap:

```r
heatmap(cor(rDat), Rowv = NA, symm = TRUE, col = brewer.pal(n = 9, name = "Blues"))
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 

Ok there are definitely groups of samples with correlated expression. I need to reorder the samples based on different variables to find if there are trends.

Check the box plots:

```r
rDatMelt <- melt(rDat, variable.name = "Sample", value.name = "RPKM")
```

```
## Using  as id variables
```

```r
head(rDatMelt)
```

```
##   Sample   RPKM
## 1   2803 7.3159
## 2   2803 7.4715
## 3   2803 1.7760
## 4   2803 0.1024
## 5   2803 5.4296
## 6   2803 0.2918
```

```r
ggplot(rDatMelt, aes(Sample, RPKM)) + geom_boxplot()
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.png) 

This is very messy, but it shows that log2 transformation is necessary:


```r
rDatMelt <- melt(log2(rDat), variable.name = "Sample", value.name = "RPKM")
```

```
## Using  as id variables
```

```r
head(rDatMelt)
```

```
##   Sample    RPKM
## 1   2803  2.8710
## 2   2803  2.9014
## 3   2803  0.8286
## 4   2803 -3.2877
## 5   2803  2.4408
## 6   2803 -1.7769
```

```r
ggplot(rDatMelt, aes(Sample, RPKM)) + geom_boxplot()
```

```
## Warning: Removed 162722 rows containing non-finite values (stat_boxplot).
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17.png) 

Ok now we can see the results! No samples appear to stick out in terms of expression. But post-log2 transformation many RPKM values become negative. I'm not sure if this is a problem?


### Differential expression analysis
I will use `voom` to perform differential expression analysis. From my experience, `voom` makes the most stringent calls for differential expression. 

**Sex**
First I will do a simple differential expression analysis: which genes are differentially expressed between males and females?

```r
sex <- rDes$Sex
table(sex)
```

```
## sex
##  F  M 
## 85 94
```


Apply scale normalization. I'm not sure if I need to do this given we have RPKM values?

```r
normFactor <- calcNormFactors(rDat)
```


Linear modelling: Use `voom` to convert RPKM to log2-RPKM ready for linear modelling:

```r
design <- model.matrix(~sex)
# The intercept represents females
head(design)
```

```
##   (Intercept) sexM
## 1           1    0
## 2           1    1
## 3           1    1
## 4           1    0
## 5           1    1
## 6           1    0
```

```r
rDatVoom <- voom(rDat, design, plot = TRUE)
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20.png) 

Looking good so far.

Now find genes differentially expressed between males and females:

```r
fit <- lmFit(rDatVoom, design)
fit <- eBayes(fit)
voomSex <- topTable(fit, coef = "sexM", p.value = 1e-05, n = Inf)
nrow(voomSex)
```

```
## [1] 30
```

```r
(voomSexgenes <- rownames(voomSex))
```

```
##  [1] "RPS4Y1|6192_calculated"       "XIST|7503_calculated"        
##  [3] "DDX3Y|8653_calculated"        "TSIX|9383_calculated"        
##  [5] "KDM5D|8284_calculated"        "EIF1AY|9086_calculated"      
##  [7] "PRKY|5616_calculated"         "UTY|7404_calculated"         
##  [9] "CYorf15A|246126_calculated"   "ZFY|7544_calculated"         
## [11] "CYorf15B|84663_calculated"    "USP9Y|8287_calculated"       
## [13] "TTTY15|64595_calculated"      "TTTY14|83869_calculated"     
## [15] "NCRNA00185|55410_calculated"  "TMSB4Y|9087_calculated"      
## [17] "TTTY10|246119_calculated"     "BCORP1|286554_calculated"    
## [19] "PRKX|5613_calculated"         "KDM5C|8242_calculated"       
## [21] "ZRSR2|8233_calculated"        "RPS4X|6191_calculated"       
## [23] "KDM6A|7403_calculated"        "ZFX|7543_calculated"         
## [25] "SRY|6736_calculated"          "NCRNA00183|554203_calculated"
## [27] "PNPLA4|8228_calculated"       "XGPY2|100132596_calculated"  
## [29] "RPS4Y2|140032_calculated"     "EIF1AX|1964_calculated"
```

Ok XIST and TSIX are popping up, this is a promising result! Plus it turns out the top hit, RPS4Y1|6192_calculated is "Ribosomal Protein S4, Y-Linked 1"! So I must be doing something right :)

Can I plot the differentially expressed genes using a DGEList object and plotSmear function from `edgeR`?

```r
# Create a DGEList object
dgeGlm <- DGEList(counts = rDat, group = sex)
plotSmear(dgeGlm, de.tags = voomSexgenes, ylab = "logFC", xlab = "AverageRPKM")
abline(h = c(-1, 1), col = "blue")
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22.png) 

I'm concerned that some of the differentially expressed transcripts being called have negative logCPM values. I'm also concerned that I'm not applying the correct tests given I'm working with RPKM values instead of CPM values :S

Let's plot the top hits to see if I'm on the right track:

```r
rDatvoomSex <- rDat[voomSexgenes[1:6], ]
dim(rDatvoomSex)
```

```
## [1]   6 179
```

```r
rDatvoomSex$Transcript <- rownames(rDatvoomSex)
rDatvoomSex <- melt(rDatvoomSex, id.vars = "Transcript", 
                   variable.name = "TCGA_patient_id",
                   value.name = "RPKM")
rDatvoomSex$Transcript <- gsub("[|].*$", "", rDatvoomSex$Transcript)
head(rDatvoomSex)
```

```
##   Transcript TCGA_patient_id     RPKM
## 1     RPS4Y1            2803   0.0000
## 2       XIST            2803 174.3137
## 3      DDX3Y            2803   0.0052
## 4       TSIX            2803  64.3488
## 5      KDM5D            2803   0.0000
## 6     EIF1AY            2803   0.0000
```

```r
rDatvoomSex <- merge(rDatvoomSex, rDes, by = "TCGA_patient_id")
head(rDatvoomSex)
```

```
##   TCGA_patient_id Transcript     RPKM Sex Race FAB_subtype Age trisomy_8
## 1            2803     RPS4Y1   0.0000   F    W          M3  61     FALSE
## 2            2803       XIST 174.3137   F    W          M3  61     FALSE
## 3            2803      DDX3Y   0.0052   F    W          M3  61     FALSE
## 4            2803       TSIX  64.3488   F    W          M3  61     FALSE
## 5            2803      KDM5D   0.0000   F    W          M3  61     FALSE
## 6            2803     EIF1AY   0.0000   F    W          M3  61     FALSE
##   del_5 del_7 Cytogenetic_risk Molecular_risk
## 1 FALSE FALSE             Good           Good
## 2 FALSE FALSE             Good           Good
## 3 FALSE FALSE             Good           Good
## 4 FALSE FALSE             Good           Good
## 5 FALSE FALSE             Good           Good
## 6 FALSE FALSE             Good           Good
```

```r
ggplot(rDatvoomSex, aes(Transcript, RPKM, colour = Sex)) +
  geom_boxplot()
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23.png) 

So the top genes differentially expressed between males and females have no expression in one of the sexes. Fair enough. But this is not a very exciting result.

**Cytogenetic risk**
Now to explore another variable, "Cytogenetic_risk". 

```r
levels(rDes$Cytogenetic_risk)
```

```
## [1] "Good"         "Intermediate" "N.D."         "Poor"
```

```r
table(rDes$Cytogenetic_risk)
```

```
## 
##         Good Intermediate         N.D.         Poor 
##           33          101            3           42
```


Which transcripts are differentially expressed between "Good", "Intermediate", and "Poor" cytogenetic risk?

First, remove samples where cytogenetic risk could not be determined "N.D":

```r
# CRGIP = Cytogenetic Response Good + Intermediate + Poor
rDesCRGIP <- droplevels(subset(rDes, Cytogenetic_risk != "N.D."))
str(rDesCRGIP)
```

```
## 'data.frame':	176 obs. of  10 variables:
##  $ TCGA_patient_id : chr  "2803" "2805" "2806" "2807" ...
##  $ Sex             : Factor w/ 2 levels "F","M": 1 2 2 1 2 2 1 2 1 2 ...
##  $ Race            : Factor w/ 13 levels "A","B","H","NH/A",..: 12 12 12 12 12 12 2 12 12 12 ...
##  $ FAB_subtype     : Factor w/ 9 levels "M0","M1","M2",..: 4 1 2 2 3 5 3 5 1 5 ...
##  $ Age             : int  61 77 46 68 23 81 25 78 39 49 ...
##  $ trisomy_8       : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
##  $ del_5           : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
##  $ del_7           : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
##  $ Cytogenetic_risk: Factor w/ 3 levels "Good","Intermediate",..: 1 2 1 2 2 2 2 3 3 1 ...
##  $ Molecular_risk  : Factor w/ 3 levels "Good","Intermediate",..: 1 2 1 2 2 2 2 3 3 1 ...
```

```r
dim(rDesCRGIP)
```

```
## [1] 176  10
```

```r
rDatCRGIP <- rDat[, rDesCRGIP$TCGA_patient_id]
dim(rDatCRGIP)
```

```
## [1] 17441   176
```

```r
identical(names(rDatCRGIP), rDesCRGIP$TCGA_patient_id)
```

```
## [1] TRUE
```

```r
cytoRisk <- rDesCRGIP$Cytogenetic_risk
```



```r
normFactor <- calcNormFactors(rDatCRGIP)
design <- model.matrix(~cytoRisk)
colnames(design)
```

```
## [1] "(Intercept)"          "cytoRiskIntermediate" "cytoRiskPoor"
```

```r
# The intercept represents 'Good' cytogenetic risk
head(design)
```

```
##   (Intercept) cytoRiskIntermediate cytoRiskPoor
## 1           1                    0            0
## 2           1                    1            0
## 3           1                    0            0
## 4           1                    1            0
## 5           1                    1            0
## 6           1                    1            0
```

```r
rDatCRGIPvoom <- voom(rDatCRGIP, design, plot = TRUE)
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-26.png) 

```r
fit <- lmFit(rDatCRGIPvoom, design)
fit <- eBayes(fit)
voomCR <- topTable(fit, coef = c("cytoRiskIntermediate", "cytoRiskPoor"), p.value = 1e-05, 
    n = Inf)
nrow(voomCR)
```

```
## [1] 774
```

```r
head(voomCR)
```

```
##                         cytoRiskIntermediate cytoRiskPoor AveExpr      F
## HOXA3|3200_calculated                  3.421        2.797  6.9596 111.89
## HOXA4|3201_calculated                  3.190        2.464  0.7177 107.69
## CPNE8|144402_calculated                3.618        3.816  2.8490 107.28
## HOXA7|3204_calculated                  4.034        3.614  2.0014 106.74
## HOXA6|3203_calculated                  4.103        3.515  3.3951 100.59
## LPO|4025_calculated                   -3.435       -3.453  5.6259  84.13
##                           P.Value adj.P.Val
## HOXA3|3200_calculated   3.938e-32 6.869e-28
## HOXA4|3201_calculated   2.576e-31 1.723e-27
## CPNE8|144402_calculated 3.091e-31 1.723e-27
## HOXA7|3204_calculated   3.952e-31 1.723e-27
## HOXA6|3203_calculated   6.689e-30 2.333e-26
## LPO|4025_calculated     2.116e-26 6.151e-23
```

```r
voomCRgenes <- rownames(voomCR)
```






