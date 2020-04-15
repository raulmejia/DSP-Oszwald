###################################
#### This script is the master program in https://github.com/raulmejia/DSP-Oszwald
#### Author: Raúl Mejía
#### The aim is to coordinate DSP data analysis 
###################################
###################################
#### Required libraries
###################################
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("pheatmap")) {
  BiocManager::install("pheatmap", ask =FALSE)
  library("pheatmap")
}
if (!require("ggplot2")) {
  BiocManager::install("ggplot2", ask =FALSE)
  library("ggplot2")
}
if (!require("argparse")) {
  install.packages("argparse", ask =FALSE)
  library("argparse")
}
if (!require("ggfortify")) {
  install.packages("ggfortify", ask =FALSE)
  library("ggfortify")
}
if (!require("sva")) {
  BiocManager::install("sva", ask =FALSE)
  library("sva")
}

###################################
#### Data given by the user
###################################
myargs <- commandArgs(trailingOnly = TRUE)
path_to_your_QC_file <-myargs[1] # Data must be separated by tab (it doesn't matter if the file has the .csv extension or other one)
path_to_your_QC_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/QC_All_Data_Human_IO_RNA.csv"

path_to_your_HK_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/HKNorm_All_Data_Human_IO_RNA.csv"
path_to_your_AreaNorm_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/AreaNorm_All_Data_Human_IO_RNA.csv"
path_to_your_SNRNorm_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/SNR_norm_All_Data_Human_IO_RNA.csv"
path_to_your_NucleiNorm_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/NucleiForm_All_Data_Human_IO_RNA.csv"

path_to_your_annotation_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Annotation/Annotation_file_Symplified.csv"
###################################
#### Reading the data
###################################
QC_table <- read.table( path_to_your_QC_file , sep = "\t", header = TRUE)
HK_table <- read.table( path_to_your_HK_file , sep = "\t", header = TRUE)
AreaNorm_table <- read.table( path_to_your_AreaNorm_file , sep = "\t", header = TRUE)
SNRNorm_table <- read.table( path_to_your_SNRNorm_file , sep = "\t", header = TRUE)
NucleiNorm_table <- read.table( path_to_your_NucleiNorm_file , sep = "\t", header = TRUE)

annot <- read.table( path_to_your_annotation_file , sep = "\t", header = TRUE)

QC_matrix <- QC_table[ ,10:dim(QC_table)[2]]
HK_matrix <- HK_table[ ,10:dim(HK_table)[2]]
AreaNorm_matrix <- AreaNorm_table[ ,10:dim(AreaNorm_table)[2]]
SNRNorm_matrix <- SNRNorm_table[ ,10:dim( SNRNorm_table )[2]]
NucleiNorm_matrix <- NucleiNorm_table[ ,10:dim( NucleiNorm_table )[2]]

mode(QC_table[,'ROI_ID']) <- "character"
mode(HK_table[,'ROI_ID']) <- "character"
mode(AreaNorm_table[,'ROI_ID']) <- "character"
mode(SNRNorm_table[,'ROI_ID']) <- "character"
mode( NucleiNorm_table[,'ROI_ID']) <- "character"

## Looking for batch effect
autoplot( prcomp( QC_matrix ), data = QC_table, colour= 'Scan_ID') +
  ggtitle("QC data")

autoplot( prcomp( HK_matrix ), data = HK_table, colour= 'Scan_ID') +
  ggtitle("Normalized trough HK")

autoplot( prcomp( AreaNorm_matrix ), data = HK_table, colour= 'Scan_ID') +
  ggtitle("Area Normalized")

autoplot( prcomp( SNRNorm_matrix ), data = HK_table, colour= 'Scan_ID') +
  ggtitle("SNR Normalized")

autoplot( prcomp( NucleiNorm_matrix ), data = HK_table, colour= 'Scan_ID') +
  ggtitle("Nuclei Normalized")

# We can observe a batch effect
# Plot the house keeping genes data

###
#############
## Looking for batch effect
#############
mymatrix <- as.matrix(SNRNorm_matrix)
mymatrix <- as.matrix(QC_matrix)
mymatrix <- as.matrix( NucleiNorm_matrix)
mymatrix <- as.matrix( HK_matrix)

pheno <- annot[, c("Morph.cat..Andre.","Histology.number")]
#pheno <- annot[, c("Morph.cat..Andre.","Scan_ID")]
head(annot)
colnames(pheno) <- c("subgroups","batch")
## Fixing colnames 
rownames(pheno) <- annot[,"Unique_ID"] ; rownames( mymatrix ) <- annot[,"Unique_ID"]
#
#
pheno$batch <- as.factor(pheno$batch)
batch<-pheno$batch
modcombat<-model.matrix(~1, data=pheno)
combat_mydata= ComBat(dat = t(mymatrix) , batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

modcombat<-model.matrix(~subgroups, data=pheno)
combat_mydata<-ComBat(dat= t(mymatrix), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

autoplot( prcomp( t(combat_mydata) ), data = HK_table, colour= 'Scan_ID') +
  ggtitle("SNR Normalized")
autoplot( prcomp( t(combat_mydata) ), data = pheno, colour= 'subgroups' , label = TRUE, label.size = 3) +
  ggtitle("SNR Normalized")

autoplot(prcomp( NucleiNorm_matrix ), data = pheno, colour= 'subgroups') +
  ggtitle("Nuclei Normalized")
autoplot( prcomp( NucleiNorm_matrix ), data = HK_table, colour= 'Scan_ID') +
  ggtitle("Nuclei Normalized")

str(combat_mydata)
head(combat_mydata)

arab <- readRDS("/media/rmejia/mountme88/Projects/DSP/Data/toy/arabidopsis.RDS") 
str(arab)
class()

class(mymatrix)
mydf[1:4,1:10]
gsub("12e RNA",12,)
mydf[,1]
str(pheno)
  