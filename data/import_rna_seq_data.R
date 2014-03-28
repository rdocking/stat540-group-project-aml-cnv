# Import RNA Seq Data

design <- read.table("clinical_patient_laml-1.tsv", sep = "\t", header = TRUE)
RNA_seq_dat <- read.table("laml.rnaseq.179_v1.0_gaf2.0_rpkm_matrix.txt.tcgaID.txt",
                          sep = "\t", header = TRUE, row.names = 1)
array_dat <- read.table("genome.wustl.edu_LAML.Genome_Wide_SNP_6.1.sdrf.txt",
                        sep = "\t", header = TRUE)

# Some inspection
str(design)
dim(design)
str(RNA_seq_dat)
dim(RNA_seq_dat)

# Let's make sure all the samples in `RNA_seq_dat` are present in `design`
des_temp <- design
des_temp$bcr_patient_barcode <- as.character(gsub("-",".",
                                                  des_temp$bcr_patient_barcode))
merged <- merge(t(RNA_seq_dat), des_temp, by.x = "row.names",
                by.y = "bcr_patient_barcode")
dim(merged) # Good, they are

# There are 179 samples (complete with clinical data) for 20442 genes
# The array data is weird... don't know how to interpret this...
str(array_dat)
head(array_dat)