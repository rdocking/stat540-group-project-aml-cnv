# NOTE: This was for the purposes of the proposal - don't bother using it now

# Import RNA Seq Data

design <- read.table("data/clinical_patient_laml.tsv", sep = "\t", 
                     header = TRUE)
RNA_seq_dat <- read.table("data/laml.rnaseq.179_v1.0_gaf2.0_rpkm_matrix.txt.tcgaID.txt",
                          sep = "\t", header = TRUE, row.names = 1)

# Some inspection
str(design)
dim(design)
str(RNA_seq_dat)
dim(RNA_seq_dat)
# There are 20442 transcripts (rows) for 179 patients (columns)

# Let's make sure all the samples in `RNA_seq_dat` are present in `design`
des_temp <- design
des_temp$bcr_patient_barcode <- as.character(gsub("-",".",
                                                  des_temp$bcr_patient_barcode))
merged <- merge(t(RNA_seq_dat), des_temp, by.x = "row.names",
                by.y = "bcr_patient_barcode")
dim(merged) # Good, they are