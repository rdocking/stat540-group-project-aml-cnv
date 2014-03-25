# Read in RNA Seq data
RNA_seq_dat <- read.table("laml.rnaseq.179_v1.0_gaf2.0_rpkm_matrix.txt.tcgaID.txt",sep="\t",header=TRUE,row.names=1)

# Prepare RNA Seq data for merging with design data
RNA_seq_dat_t <- as.data.frame(t(RNA_seq_dat))
rownames(RNA_seq_dat_t) <- regmatches(rownames(RNA_seq_dat_t),regexpr("(?<=AB.).*",rownames(RNA_seq_dat_t),perl=TRUE))

# Merge RNA Seq and design
joined_dat <- merge(cleaned_data,RNA_seq_dat_t,by.x="TCGA_patient_id",by.y="row.names")