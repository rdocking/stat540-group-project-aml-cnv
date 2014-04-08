Random forest exploratory analysis
========================================================

Load the data and metadata. In this analysis, I will be using the cleaned RPKM data.

Load libraries:

```r
library(kernlab)
library(cvTools)
```

```
## Loading required package: lattice
## Loading required package: robustbase
```

```r
library(randomForest)
```

```
## randomForest 4.6-7
## Type rfNews() to see new features/changes/bug fixes.
```

```r
library(edgeR)
```

```
## Loading required package: limma
```

```r
library(caret)
```

```
## Warning: package 'caret' was built under R version 3.0.3
```

```
## Loading required package: ggplot2
```

```r
library(plyr)
library(limma)
library(VennDiagram)
```

```
## Loading required package: grid
```



```r
rDes <- read.delim("../../data/experimental_design_cleaned.txt")
rownames(rDes) <- rDes$TCGA_patient_id
rDat <- read.delim("../../data/aml.rnaseq.gaf2.0_rpkm_cleaned.txt", row.names = 1, 
    check.names = FALSE)

all.results <- list()
```



## 1) Functions for cross-validation

```r
# Function to select features using linear models input.dat: training data
# set input.labels: true outcomes for training data
fs.lm <- function(input.dat, input.labels) {
    norm.factor <- calcNormFactors(input.dat)
    design <- model.matrix(~input.labels)
    colnames(design) <- c("Intercept", "Label")
    dat.voomed <- voom(input.dat, design, lib.size = colSums(input.dat) * norm.factor)
    fit <- lmFit(dat.voomed, design)
    ebFit <- eBayes(fit)
    hits <- topTable(ebFit, n = Inf, coef = "Label")
    train.features <- rownames(hits)[1:25]
    return(train.features)
}
```




```r
# Function to choose features using correlations to outcomes input.dat:
# training data set input.labels: true outcomes for training data
fs.corr <- function(input.dat, input.levels) {
    gene.corrs <- apply(input.dat, 1, function(x) return(suppressWarnings(cor(x, 
        as.numeric(input.levels), method = "spearman"))))
    gene.corrs <- gene.corrs[order(abs(gene.corrs), na.last = NA, decreasing = TRUE)]
    return(names(gene.corrs[1:25]))
}
```



```r
# Function to run k-fold cross validation with a random forest all.dat: all
# data used in the analysis all.labels: true outcomes for the data K: number
# of folds to use in CV fs.method: the strategy to use for feature selection
rf.cv <- function(all.dat, all.labels, all.levels, K = 5, fs.method = "lm", 
    plot.varimp = TRUE) {
    set.seed(540)
    folds <- cvFolds(ncol(all.dat), K = K)
    
    conf_matrix <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("true0", 
        "true1"), c("pred0", "pred1")))
    feature.list <- list()
    
    for (f in 1:K) {
        train.samples <- folds$subsets[folds$which != f, ]
        test.samples <- folds$subsets[folds$which == f, ]
        
        if (fs.method == "lm") {
            train.features <- fs.lm(all.dat[, train.samples], all.labels[train.samples])
        } else if (fs.method == "corr") {
            train.features <- fs.corr(all.dat[, train.samples], all.levels[train.samples])
        }
        feature.list[[paste("Fold", f, sep = "")]] <- train.features
        
        train.dat <- data.frame(class = all.labels[train.samples], t(all.dat[train.features, 
            train.samples]))
        test.dat <- data.frame(t(all.dat[train.features, test.samples]))
        test.labels <- all.labels[test.samples]
        
        fit.rf <- randomForest(class ~ ., train.dat)
        if (plot.varimp) {
            varImpPlot(fit.rf)
        }
        
        pred.rf <- predict(fit.rf, newdata = test.dat, type = "response")
        
        results <- table(factor(test.labels, levels = c(0, 1)), factor(pred.rf, 
            levels = c(0, 1)), dnn = c("obs", "pred"))
        
        conf_matrix[1, 1] <- conf_matrix[1, 1] + results[1, 1]
        conf_matrix[1, 2] <- conf_matrix[1, 2] + results[1, 2]
        conf_matrix[2, 1] <- conf_matrix[2, 1] + results[2, 1]
        conf_matrix[2, 2] <- conf_matrix[2, 2] + results[2, 2]
    }
    
    rf.sens <- conf_matrix[2, 2]/sum(conf_matrix[2, ])
    rf.spec <- conf_matrix[1, 1]/sum(conf_matrix[1, ])
    rf.acc <- (conf_matrix[1, 1] + conf_matrix[2, 2])/sum(conf_matrix)
    
    return(list(acc = rf.acc, sens = rf.sens, spec = rf.spec, feature.list = feature.list))
}
```



## 2) Train RF to predict "Poor" cytogenetic risk, use linear models for feature selection

Set up the data. Remove samples where the cytogenetic risk category is not determined, and set categories to poor and not poor:

```r
rfDes <- rDes[rDes$Cytogenetic_risk != "N.D.", ]
rf.labels <- mapvalues(rfDes$Cytogenetic_risk, c("Good", "Intermediate", "Poor"), 
    c(0, 0, 1), warn_missing = TRUE)
rf.labels <- factor(rf.labels)
rf.levels <- mapvalues(rfDes$Cytogenetic_risk, c("Good", "Intermediate", "Poor"), 
    c(3, 2, 1), warn_missing = TRUE)
rf.levels <- factor(rf.levels)
rfDat <- rDat[, rownames(rfDes)]
```


Run a 5-fold cross-validation for the data:

```r
cv.risk.res <- rf.cv(rfDat, rf.labels, rf.levels, K = 5, fs.method = "lm")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-71.png) ![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-72.png) ![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-73.png) ![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-74.png) ![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-75.png) 

```r
cv.risk.res[1:3]
```

```
## $acc
## [1] 0.8523
## 
## $sens
## [1] 0.5714
## 
## $spec
## [1] 0.9403
```

```r
all.results[["lm.poor"]] <- cv.risk.res[1:3]
fts <- cv.risk.res[[4]]

plot.new()
venn.plot <- venn.diagram(fts, filename = NULL, fill = c("red", "blue", "green", 
    "yellow", "purple"), margin = 0.1)
grid.draw(venn.plot)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-76.png) 

```r

(common.fts.risk <- intersect(fts[[1]], intersect(fts[[2]], intersect(fts[[3]], 
    intersect(fts[[4]], fts[[5]])))))
```

```
## [1] "SCD|6319_calculated"     "STYXL1|51657_calculated"
## [3] "PDAP1|11333_calculated"  "LUC7L2|51631_calculated"
## [5] "GSTK1|373156_calculated"
```



## 3) Predict the different cytogenetic mutations, use linear models for feature selection

First look at trisomy 8:

```r
cv.trisomy8.res <- rf.cv(rDat, factor(as.numeric(rDes$trisomy_8)), rf.levels, 
    K = 5, fs.method = "lm")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-81.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-82.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-83.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-84.png) ![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-85.png) 

```r
cv.trisomy8.res[1:3]
```

```
## $acc
## [1] 0.9553
## 
## $sens
## [1] 0.6842
## 
## $spec
## [1] 0.9875
```

```r
all.results[["lm.trisomy8"]] <- cv.trisomy8.res[1:3]
fts <- cv.trisomy8.res[[4]]

plot.new()
venn.plot <- venn.diagram(fts, filename = NULL, fill = c("red", "blue", "green", 
    "yellow", "purple"), margin = 0.1)
grid.draw(venn.plot)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-86.png) 

```r

(common.fts.trisomy8 <- intersect(fts[[1]], intersect(fts[[2]], intersect(fts[[3]], 
    intersect(fts[[4]], fts[[5]])))))
```

```
## [1] "NEIL2|252969_calculated"   "PPP2R2A|5520_calculated"  
## [3] "ZNF7|7553_calculated"      "KIAA1967|57805_calculated"
## [5] "R3HCC1|203069_calculated"  "TSNARE1|203062_calculated"
## [7] "C8orf55|51337_calculated"  "ZFP41|286128_calculated"
```


Next look at deletions of chromosome 5:

```r
cv.del5.res <- rf.cv(rDat, factor(as.numeric(rDes$del_5)), rf.levels, K = 5, 
    fs.method = "lm")
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-91.png) ![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-92.png) ![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-93.png) ![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-94.png) ![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-95.png) 

```r
cv.del5.res[1:3]
```

```
## $acc
## [1] 0.9721
## 
## $sens
## [1] 0.8125
## 
## $spec
## [1] 0.9877
```

```r
all.results[["lm.del5"]] <- cv.del5.res[1:3]
fts <- cv.del5.res[[4]]

pdf("../../poster/rf_venn.pdf")
plot.new()
venn.plot <- venn.diagram(fts, filename = NULL, fill = c("red", "blue", "green", 
    "yellow", "purple"), margin = 0.1)
grid.draw(venn.plot)
dev.off()
```

```
## pdf 
##   2
```

```r

(common.fts.del5 <- intersect(fts[[1]], intersect(fts[[2]], intersect(fts[[3]], 
    intersect(fts[[4]], fts[[5]])))))
```

```
##  [1] "KDM3B|51780_calculated"        "EIF4EBP3|8637_calculated"     
##  [3] "PCBD2|84105_calculated"        "KIAA0141|9812|1of2_calculated"
##  [5] "PFDN1|5201_calculated"         "RBM22|55696_calculated"       
##  [7] "ZMAT2|153527_calculated"       "CSNK1A1|1452_calculated"      
##  [9] "PPP2CA|5515_calculated"        "WDR55|54853_calculated"       
## [11] "CATSPER3|347732_calculated"    "HARS|3035_calculated"
```


Next look at deletions of chromosome 7:

```r
cv.del7.res <- rf.cv(rDat, factor(as.numeric(rDes$del_7)), rf.levels, K = 5, 
    fs.method = "lm")
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-101.png) ![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-102.png) ![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-103.png) ![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-104.png) ![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-105.png) 

```r
cv.del7.res[1:3]
```

```
## $acc
## [1] 0.9497
## 
## $sens
## [1] 0.7143
## 
## $spec
## [1] 0.981
```

```r
all.results[["lm.del7"]] <- cv.del7.res[1:3]
fts <- cv.del7.res[[4]]

plot.new()
venn.plot <- venn.diagram(fts, filename = NULL, fill = c("red", "blue", "green", 
    "yellow", "purple"), margin = 0.1)
grid.draw(venn.plot)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-106.png) 

```r

(common.fts.del7 <- intersect(fts[[1]], intersect(fts[[2]], intersect(fts[[3]], 
    intersect(fts[[4]], fts[[5]])))))
```

```
## [1] "LUC7L2|51631_calculated"   "PDAP1|11333_calculated"   
## [3] "MKRN1|23608_calculated"    "SLC25A13|10165_calculated"
## [5] "GATAD1|57798_calculated"   "GSTK1|373156_calculated"  
## [7] "STYXL1|51657_calculated"   "SUMF2|25870_calculated"   
## [9] "CASP2|835_calculated"
```


Now get a sense of how many genes are shared between these 4 different sets:

```r
fts.all <- list(risk = common.fts.risk, trisomy8 = common.fts.trisomy8, del5 = common.fts.del5, 
    del7 = common.fts.del7)
plot.new()
venn.plot <- venn.diagram(fts.all, filename = NULL, fill = c("red", "yellow", 
    "blue", "green"))
grid.draw(venn.plot)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11.png) 



## 4) Train RF to predict cytogenetic risk, use correlations for feature selection

Run a 5-fold cross-validation for the data:

```r
cv.risk.res <- rf.cv(rfDat, rf.labels, rf.levels, K = 5, fs.method = "corr")
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-121.png) ![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-122.png) ![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-123.png) ![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-124.png) ![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-125.png) 

```r
cv.risk.res[1:3]
```

```
## $acc
## [1] 0.8239
## 
## $sens
## [1] 0.381
## 
## $spec
## [1] 0.9627
```

```r
all.results[["corr.poor"]] <- cv.risk.res[1:3]
fts <- cv.risk.res[[4]]

plot.new()
venn.plot <- venn.diagram(fts, filename = NULL, fill = c("red", "blue", "green", 
    "yellow", "purple"), margin = 0.1)
grid.draw(venn.plot)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-126.png) 

```r

(common.fts.risk <- intersect(fts[[1]], intersect(fts[[2]], intersect(fts[[3]], 
    intersect(fts[[4]], fts[[5]])))))
```

```
## [1] "PDE4DIP|9659_calculated"  "PHKA1|5255_calculated"   
## [3] "SDPR|8436_calculated"     "RECK|8434_calculated"    
## [5] "IL7|3574_calculated"      "STARD10|10809_calculated"
```



## 5) Predict the different cytogenetic mutations, use correlations for feature selection

First look at trisomy 8:

```r
cv.trisomy8.res <- rf.cv(rDat, factor(as.numeric(rDes$trisomy_8)), factor(as.numeric(rDes$trisomy_8)), 
    K = 5, fs.method = "corr")
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-131.png) ![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-132.png) ![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-133.png) ![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-134.png) ![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-135.png) 

```r
cv.trisomy8.res[1:3]
```

```
## $acc
## [1] 0.9609
## 
## $sens
## [1] 0.6842
## 
## $spec
## [1] 0.9938
```

```r
all.results[["cor.trisomy8"]] <- cv.trisomy8.res[1:3]
fts <- cv.trisomy8.res[[4]]

plot.new()
venn.plot <- venn.diagram(fts, filename = NULL, fill = c("red", "blue", "green", 
    "yellow", "purple"), margin = 0.1)
grid.draw(venn.plot)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-136.png) 

```r

(common.fts.trisomy8 <- intersect(fts[[1]], intersect(fts[[2]], intersect(fts[[3]], 
    intersect(fts[[4]], fts[[5]])))))
```

```
##  [1] "NEIL2|252969_calculated"   "KIAA1967|57805_calculated"
##  [3] "R3HCC1|203069_calculated"  "DOCK5|80005_calculated"   
##  [5] "BIN3|55909_calculated"     "POLR3D|661_calculated"    
##  [7] "PPP2R2A|5520_calculated"   "TSNARE1|203062_calculated"
##  [9] "ZNF7|7553_calculated"      "HEATR7A|727957_calculated"
## [11] "COMMD5|28991_calculated"   "IKBKB|3551_calculated"
```


Next look at deletions of chromosome 5:

```r
cv.del5.res <- rf.cv(rDat, factor(as.numeric(rDes$del_5)), factor(as.numeric(rDes$del_5)), 
    K = 5, fs.method = "corr")
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-141.png) ![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-142.png) ![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-143.png) ![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-144.png) ![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-145.png) 

```r
cv.del5.res[1:3]
```

```
## $acc
## [1] 0.9385
## 
## $sens
## [1] 0.4375
## 
## $spec
## [1] 0.9877
```

```r
all.results[["cor.del5"]] <- cv.del5.res[1:3]
fts <- cv.del5.res[[4]]

plot.new()
venn.plot <- venn.diagram(fts, filename = NULL, fill = c("red", "blue", "green", 
    "yellow", "purple"), margin = 0.1)
grid.draw(venn.plot)
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-146.png) 

```r

(common.fts.del5 <- intersect(fts[[1]], intersect(fts[[2]], intersect(fts[[3]], 
    intersect(fts[[4]], fts[[5]])))))
```

```
## [1] "DSCAM|1826_calculated"
```


Next look at deletions of chromosome 7:

```r
cv.del7.res <- rf.cv(rDat, factor(as.numeric(rDes$del_7)), factor(as.numeric(rDes$del_7)), 
    K = 5, fs.method = "corr")
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-151.png) ![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-152.png) ![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-153.png) ![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-154.png) ![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-155.png) 

```r
cv.del7.res[1:3]
```

```
## $acc
## [1] 0.9497
## 
## $sens
## [1] 0.7143
## 
## $spec
## [1] 0.981
```

```r
all.results[["cor.del7"]] <- cv.del7.res[1:3]
fts <- cv.del7.res[[4]]

plot.new()
venn.plot <- venn.diagram(fts, filename = NULL, fill = c("red", "blue", "green", 
    "yellow", "purple"), margin = 0.1)
grid.draw(venn.plot)
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-156.png) 

```r

(common.fts.del7 <- intersect(fts[[1]], intersect(fts[[2]], intersect(fts[[3]], 
    intersect(fts[[4]], fts[[5]])))))
```

```
## [1] "MKRN1|23608_calculated"  "LUC7L2|51631_calculated"
## [3] "GSTK1|373156_calculated" "PDAP1|11333_calculated" 
## [5] "FSCN3|29999_calculated"  "CASP2|835_calculated"   
## [7] "PRKAG2|51422_calculated" "CPSF4|10898_calculated" 
## [9] "ARF5|381_calculated"
```


Now get a sense of how many genes are shared between these 4 different sets:

```r
fts.all <- list(risk = common.fts.risk, trisomy8 = common.fts.trisomy8, del5 = common.fts.del5, 
    del7 = common.fts.del7)
plot.new()
venn.plot <- venn.diagram(fts.all, filename = NULL, fill = c("red", "yellow", 
    "blue", "green"))
grid.draw(venn.plot)
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.png) 



## 6) Train RF to predict "Good" cytogenetic risk, use linear models for feature selection

Set up the data. Remove samples where the cytogenetic risk category is not determined, and set categories to good and not good:

```r
rfDes <- rDes[rDes$Cytogenetic_risk != "N.D.", ]
rf.labels <- mapvalues(rfDes$Cytogenetic_risk, c("Good", "Intermediate", "Poor"), 
    c(1, 0, 0), warn_missing = TRUE)
rf.labels <- factor(rf.labels)
rf.levels <- mapvalues(rfDes$Cytogenetic_risk, c("Good", "Intermediate", "Poor"), 
    c(3, 2, 1), warn_missing = TRUE)
rf.levels <- factor(rf.levels)
rfDat <- rDat[, rownames(rfDes)]
```


Run a 5-fold cross-validation for the data:

```r
cv.risk.res <- rf.cv(rfDat, rf.labels, rf.levels, K = 5, fs.method = "lm")
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-181.png) ![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-182.png) ![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-183.png) ![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-184.png) ![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-185.png) 

```r
cv.risk.res[1:3]
```

```
## $acc
## [1] 0.9659
## 
## $sens
## [1] 0.8485
## 
## $spec
## [1] 0.993
```

```r
all.results[["lm.good"]] <- cv.risk.res[1:3]
fts <- cv.risk.res[[4]]

plot.new()
venn.plot <- venn.diagram(fts, filename = NULL, fill = c("red", "blue", "green", 
    "yellow", "purple"), margin = 0.1)
grid.draw(venn.plot)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-186.png) 

```r

(common.fts.risk <- intersect(fts[[1]], intersect(fts[[2]], intersect(fts[[3]], 
    intersect(fts[[4]], fts[[5]])))))
```

```
##  [1] "CPNE8|144402_calculated"  "HOXA7|3204_calculated"   
##  [3] "HOXA6|3203_calculated"    "HOXA5|3202_calculated"   
##  [5] "HOXA3|3200_calculated"    "HOXA4|3201_calculated"   
##  [7] "HOXA9|3205_calculated"    "HOXA10|3206_calculated"  
##  [9] "HOXA2|3199_calculated"    "FGFR1|2260_calculated"   
## [11] "CYP7B1|9420_calculated"   "PDE4DIP|9659_calculated" 
## [13] "HOXB5|3215_calculated"    "NKX2-3|159296_calculated"
## [15] "LPO|4025_calculated"      "RMND5B|64777_calculated" 
## [17] "PRDM16|63976_calculated"  "HOXB6|3216_calculated"
```



## 7) Train RF to predict "Good" cytogenetic risk, use correlations for feature selection

```r
cv.risk.res <- rf.cv(rfDat, rf.labels, rf.levels, K = 5, fs.method = "corr")
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-191.png) ![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-192.png) ![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-193.png) ![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-194.png) ![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-195.png) 

```r
cv.risk.res[1:3]
```

```
## $acc
## [1] 0.9602
## 
## $sens
## [1] 0.8485
## 
## $spec
## [1] 0.986
```

```r
all.results[["corr.good"]] <- cv.risk.res[1:3]
fts <- cv.risk.res[[4]]

plot.new()
venn.plot <- venn.diagram(fts, filename = NULL, fill = c("red", "blue", "green", 
    "yellow", "purple"), margin = 0.1)
grid.draw(venn.plot)
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-196.png) 

```r

(common.fts.risk <- intersect(fts[[1]], intersect(fts[[2]], intersect(fts[[3]], 
    intersect(fts[[4]], fts[[5]])))))
```

```
## [1] "PDE4DIP|9659_calculated"  "PHKA1|5255_calculated"   
## [3] "SDPR|8436_calculated"     "RECK|8434_calculated"    
## [5] "IL7|3574_calculated"      "STARD10|10809_calculated"
```



## 8) Train RF to predict "Intermediate" cytogenetic risk, use linear models for feature selection

Set up the data. Remove samples where the cytogenetic risk category is not determined, and set categories to intermediate and not intermediate:

```r
rfDes <- rDes[rDes$Cytogenetic_risk != "N.D.", ]
rf.labels <- mapvalues(rfDes$Cytogenetic_risk, c("Good", "Intermediate", "Poor"), 
    c(0, 1, 0), warn_missing = TRUE)
rf.labels <- factor(rf.labels)
rf.levels <- mapvalues(rfDes$Cytogenetic_risk, c("Good", "Intermediate", "Poor"), 
    c(3, 2, 1), warn_missing = TRUE)
rf.levels <- factor(rf.levels)
rfDat <- rDat[, rownames(rfDes)]
```


Run a 5-fold cross-validation for the data:

```r
cv.risk.res <- rf.cv(rfDat, rf.labels, rf.levels, K = 5, fs.method = "lm")
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-211.png) ![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-212.png) ![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-213.png) ![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-214.png) ![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-215.png) 

```r
cv.risk.res[1:3]
```

```
## $acc
## [1] 0.8636
## 
## $sens
## [1] 0.901
## 
## $spec
## [1] 0.8133
```

```r
all.results[["lm.intermediate"]] <- cv.risk.res[1:3]
fts <- cv.risk.res[[4]]

plot.new()
venn.plot <- venn.diagram(fts, filename = NULL, fill = c("red", "blue", "green", 
    "yellow", "purple"), margin = 0.1)
grid.draw(venn.plot)
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-216.png) 

```r

(common.fts.risk <- intersect(fts[[1]], intersect(fts[[2]], intersect(fts[[3]], 
    intersect(fts[[4]], fts[[5]])))))
```

```
##  [1] "NAV1|89796_calculated"    "SLC18A2|6571_calculated" 
##  [3] "LASS4|79603_calculated"   "PBX3|5090_calculated"    
##  [5] "HOXB6|3216_calculated"    "HOXB5|3215_calculated"   
##  [7] "NKX2-3|159296_calculated" "IQCE|23288_calculated"   
##  [9] "EVPL|2125_calculated"     "C7orf50|84310_calculated"
```



## 9) Train RF to predict "Intermediate" cytogenetic risk, use correlations for feature selection

```r
cv.risk.res <- rf.cv(rfDat, rf.labels, rf.levels, K = 5, fs.method = "corr")
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-221.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-222.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-223.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-224.png) ![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-225.png) 

```r
cv.risk.res[1:3]
```

```
## $acc
## [1] 0.7784
## 
## $sens
## [1] 0.901
## 
## $spec
## [1] 0.6133
```

```r
all.results[["corr.intermediate"]] <- cv.risk.res[1:3]
fts <- cv.risk.res[[4]]

plot.new()
venn.plot <- venn.diagram(fts, filename = NULL, fill = c("red", "blue", "green", 
    "yellow", "purple"), margin = 0.1)
grid.draw(venn.plot)
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-226.png) 

```r

(common.fts.risk <- intersect(fts[[1]], intersect(fts[[2]], intersect(fts[[3]], 
    intersect(fts[[4]], fts[[5]])))))
```

```
## [1] "PDE4DIP|9659_calculated"  "PHKA1|5255_calculated"   
## [3] "SDPR|8436_calculated"     "RECK|8434_calculated"    
## [5] "IL7|3574_calculated"      "STARD10|10809_calculated"
```



## 10) Summarize results

```r
show(all.results)
```

```
## $lm.poor
## $lm.poor$acc
## [1] 0.8523
## 
## $lm.poor$sens
## [1] 0.5714
## 
## $lm.poor$spec
## [1] 0.9403
## 
## 
## $lm.trisomy8
## $lm.trisomy8$acc
## [1] 0.9553
## 
## $lm.trisomy8$sens
## [1] 0.6842
## 
## $lm.trisomy8$spec
## [1] 0.9875
## 
## 
## $lm.del5
## $lm.del5$acc
## [1] 0.9721
## 
## $lm.del5$sens
## [1] 0.8125
## 
## $lm.del5$spec
## [1] 0.9877
## 
## 
## $lm.del7
## $lm.del7$acc
## [1] 0.9497
## 
## $lm.del7$sens
## [1] 0.7143
## 
## $lm.del7$spec
## [1] 0.981
## 
## 
## $corr.poor
## $corr.poor$acc
## [1] 0.8239
## 
## $corr.poor$sens
## [1] 0.381
## 
## $corr.poor$spec
## [1] 0.9627
## 
## 
## $cor.trisomy8
## $cor.trisomy8$acc
## [1] 0.9609
## 
## $cor.trisomy8$sens
## [1] 0.6842
## 
## $cor.trisomy8$spec
## [1] 0.9938
## 
## 
## $cor.del5
## $cor.del5$acc
## [1] 0.9385
## 
## $cor.del5$sens
## [1] 0.4375
## 
## $cor.del5$spec
## [1] 0.9877
## 
## 
## $cor.del7
## $cor.del7$acc
## [1] 0.9497
## 
## $cor.del7$sens
## [1] 0.7143
## 
## $cor.del7$spec
## [1] 0.981
## 
## 
## $lm.good
## $lm.good$acc
## [1] 0.9659
## 
## $lm.good$sens
## [1] 0.8485
## 
## $lm.good$spec
## [1] 0.993
## 
## 
## $corr.good
## $corr.good$acc
## [1] 0.9602
## 
## $corr.good$sens
## [1] 0.8485
## 
## $corr.good$spec
## [1] 0.986
## 
## 
## $lm.intermediate
## $lm.intermediate$acc
## [1] 0.8636
## 
## $lm.intermediate$sens
## [1] 0.901
## 
## $lm.intermediate$spec
## [1] 0.8133
## 
## 
## $corr.intermediate
## $corr.intermediate$acc
## [1] 0.7784
## 
## $corr.intermediate$sens
## [1] 0.901
## 
## $corr.intermediate$spec
## [1] 0.6133
```

