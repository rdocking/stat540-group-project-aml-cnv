#Differential Expression analysis 1

#Preliminaries and loading data
library(lattice) # if you don't already have this loaded ...
library(ggplot2) # we'll make figures with both
library(reshape2) # for the function melt
library(edgeR)

cleaned_data <- readRDS("cleaned_data.rds") #might need to fix pathname
RNA_seq_dat <- read.table("laml.rnaseq.179_v1.0_gaf2.0_rpkm_matrix.txt.tcgaID.txt",
                          sep = "\t", header = TRUE, row.names = 1) #might need to fix pathname

#making mini design matrix with just T8, D5 D7. Using as. numeric so it can be used as a design matrix for limma/voom--------
str(miniDes <- cleaned_data[,cbind("TCGA_patient_id","trisomy_8", "del_5", "del_7")])
miniDes$del_7 <- as.numeric(miniDes$del_7, levels = c("TRUE", "FALSE"))
miniDes$del_5 <- as.numeric(miniDes$del_5, levels = c("TRUE", "FALSE"))
miniDes$trisomy_8 <- as.numeric(miniDes$trisomy_8, levels = c("TRUE", "FALSE"))

#sample IDs are ordered differently in design and data matricies, so they will be merged by row name------
#first had to modify TCGA_patient_id in the design matrix, adding "TCGA.AB." to patient ID so it matches the data matrix row names.
des_temp <- miniDes
des_temp$bcr_patient_barcode <- cbind(paste0("TCGA.AB.", des_temp$TCGA_patient_id)) #added "TCGA.AB." to patient ID
merged <- merge(t(RNA_seq_dat), des_temp, by.x = "row.names",
                by.y = "bcr_patient_barcode") #merged the two matries by row name
head(merged[1:10]) #Design matrix got merged to rear of data matrix

#putting the design matrix portion in front of the data matrix portion in the merged matrix
ordMerged <- cbind(merged[1], merged[20444:20447], merged[2:20443])

#Reordering rows of the merged matrix (ordMerged), grouping by T8 D5 D7.
reordMerged <- ordMerged[order(ordMerged[,"trisomy_8"], ordMerged[,"del_5"], ordMerged[,"del_7"]), ]
rownames(reordMerged) <- c(seq_along(rownames(reordMerged)))
reordMerged[1:5] #there is some overlap in the groups. Therefore, they can't be reduced to levels of a single factor.

# if we made each combination of conditions its own group, there would be 7 groups.
# will try to avoid that for now... though, maybe that's what I should be doing? and then combine them with contrasts later?

#remaking separate design and data matrices from merged matrix for processing with limma and voom.---------
(prDes <- reordMerged[1:5])
prDat <- ordMerged[6:ncol(reordMerged)]

#now we have our design and data, and we can walk through seminar 7.--------------
#can we even do glm analysis if we can't make 'group'?------------------

dge.glm <- DGEList(counts = prDat, group = design3) #looks like not... what if there's no sample ovelap?
dge.glm <- DGEList(counts = prDat, group = design3[1:2]) #OK: DGE groups seem to be for 1 factor with multiples levels. Not our case.

#testing out VOOM.------------
#discovering genes with 0 counts prevents normalization
norm.factor <- calcNormFactors(prDat[1:16]) #Error message. so prDat(16) has some kind of missing value.
str(prDat[1:16]) #seems unable to deal with a column full of zeros. go figure. So, removing zeros....
prDat[15]

#RPKM filtering, based off of Mac's script (Mac does RNA-seq in the bohlmann lab)--------------
y <- DGEList(counts=t(prDat)) #had to transpose from original disordered form...
nrow(y)

# Keep only genes with at least 1 count-per-million reads (cpm) in at least 4 samples
z <- y[(rowSums(cpm(y) > 1) >= 4), ]

# Reset depth
z$samples$lib.size <- colSums(z$counts)

#Normalize by Depth (back to seminar) ------
norm.factor <- calcNormFactors(z$counts)

#preparing design matrix with reference + treatment effects parameterization
design2 <- prDes[3:5]
rownames(design2) <- prDes$Row.names
design2$Intercept <- c(rep(1, 179))
design3 <- cbind(design2[4], design2[1:3])
design3

#Voom and limma analysis
dat.voomed <- voom(z$counts, design3, plot = TRUE, lib.size = colSums(z$counts) * norm.factor)

fit <- lmFit(dat.voomed, design3)
fit <- eBayes(fit)
head(ttfit <- topTable(fit, number = Inf)) #oh, hurray!

#checking out F-value distribution
densityplot(~F, ttfit,
            plot.points=TRUE,
            main='Density plot of F-values',
            bw =0.00002, n = 1000,
            xlab='F-value')

#what is adj. P value even saying? This is what we used in homework 2, 
#but all genes have a huge adj. P value, so I don't trust what's going on.
densityplot(~adj.P.Val, ttfit,
            plot.points=TRUE,
            main='Density plot of adj.P.Values',
            n = 100000,
            xlab='adj.P.value')

#My Questions: 
#how to interpret adj. P value here (because it's ridiculously low for everything)
#Whether I should have made 7 groups, and then used contrast matrices to do toptable comparisons between them.

#optional: Making a Tall matrix------------
prMDat <- melt(ordMerged,
               id.vars=c('Row.names', 'TCGA_patient_id', 'trisomy_8', 'del_5', 'del_7'),
               variable.name='transcript', value.name='counts')
str(prMDat)

