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

###################################
#### Data given by the user
###################################
myargs <- commandArgs(trailingOnly = TRUE)
path_to_your_QC_file <-myargs[1] # Data must be separated by tab (it doesn't matter if the file has the .csv extension or other one)
path_to_your_QC_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/QC_All_Data_Human_IO_RNA.csv"

path_to_your_HK_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/HKNorm_All_Data_Human_IO_RNA.csv"
path_to_your_AreaNorm_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/AreaNorm_All_Data_Human_IO_RNA.csv"
path_to_your_SNRNorm_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/SNR_norm_All_Data_Human_IO_RNA.csv"

###################################
#### Reading the data
###################################
QC_table <- read.table( path_to_your_QC_file , sep = "\t", header = TRUE)
HK_table <- read.table( path_to_your_HK_file , sep = "\t", header = TRUE)
AreaNorm_table <- read.table( path_to_your_AreaNorm_file , sep = "\t", header = TRUE)
SNRNorm_table <- read.table( path_to_your_SNRNorm_file , sep = "\t", header = TRUE)

dim(QC_table)[2]
QC_matrix <- QC_table[ ,10:dim(QC_table)[2]]
HK_matrix <- HK_table[ ,10:dim(HK_table)[2]]
AreaNorm_matrix <- AreaNorm_table[ ,10:dim(AreaNorm_table)[2]]
SNRNorm_matrix <- SNRNorm_table[ ,10:dim( SNRNorm_table )[2]]

mode(QC_table[,'ROI_ID']) <- "character"
mode(HK_table[,'ROI_ID']) <- "character"
mode(AreaNorm_table[,'ROI_ID']) <- "character"
mode(SNRNorm_table[,'ROI_ID']) <- "character"

## Looking for batch effect
autoplot( prcomp( QC_matrix ), data = QC_table, colour= 'Scan_ID') +
  ggtitle("QC data")

autoplot( prcomp( HK_matrix ), data = HK_table, colour= 'Scan_ID') +
  ggtitle("Normalized trough HK")

autoplot( prcomp( AreaNorm_matrix ), data = HK_table, colour= 'Scan_ID') +
  ggtitle("Area Normalized")

autoplot( prcomp( SNRNorm_matrix ), data = HK_table, colour= 'Scan_ID') +
  ggtitle("SNR Normalized")


# We can observe a batch effect
# Plot the house keeping genes data


###
#############
## Looking for batch effect
#############



str(iris[,1:4])

data(iris)
head(iris)
dim(iris)
extra_column <-c(rep("green",30),rep("yellow",30),rep("black",30),rep("purple",30),rep("brown",30))
iris2 <- cbind(iris, extra_column)
head(iris2)

autoplot( prcomp(iris2[,1:4]), data = iris2, colour= 'Species' )
autoplot( prcomp(iris2[,1:4]), data = iris2, colour= 'Species' , shape= 'extra_column', label =TRUE)

