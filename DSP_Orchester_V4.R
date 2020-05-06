###################################
#### This script is the master program in https://github.com/raulmejia/DSP-Oszwald
#### Author: Raúl Mejía
#### The aim is to coordinate DSP data analysis 
###################################
###################################
#### 0) loading and/or installing required libraries
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
if (!require("limma")) {
  BiocManager::install("limma", ask =FALSE)
  library("limma")
}
if (!require("smooth")) {
  install.packages("smooth", ask =FALSE)
  library("smooth")
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", ask =FALSE)
  library("RColorBrewer")
}
if (!require("plotrix")) {
  install.packages("plotrix", ask =FALSE)
  library("plotrix")
}
if (!require("datasets")) {
  install.packages("datasets", ask =FALSE)
  library("datasets")
}
if (!require("reshape2")) {
  install.packages("reshape2", ask =FALSE)
  library("reshape2")
}
if (!require("ggridges")) {
  install.packages("ggridges", ask =FALSE)
  library("ggridges")
}
if (!require("cowplot")) {
  install.packages("cowplot", ask =FALSE)
  library("cowplot")
}
if (!require("preprocessCore")) {
  BiocManager::install("preprocessCore", ask =FALSE)
  library("preprocessCore")
}
if (!require("affy")) {
  BiocManager::install("affy", ask =FALSE)
  library("affy")
}

###################################
#### Data given by the user
###################################
#myargs <- commandArgs(trailingOnly = TRUE)
#path_to_your_QC_file <-myargs[1] # Data must be separated by tab (it doesn't matter if the file has the .csv extension or other one)
#path_to_your_QC_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/QC_All_Data_Human_IO_RNA.csv"
#path_to_your_HK_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/HKNorm_All_Data_Human_IO_RNA.csv"
#path_to_your_AreaNorm_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/AreaNorm_All_Data_Human_IO_RNA.csv"
#path_to_your_SNRNorm_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/SNR_norm_All_Data_Human_IO_RNA.csv"
#path_to_your_NucleiNorm_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/NucleiForm_All_Data_Human_IO_RNA.csv"

path_to_your_annotation_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Annotation/Annotation_file_Symplified_Corrected_colnames_NOspaces.csv"
#path_to_your_annotation_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Annotation/DSP_ROI_annotations_outcome_v2_RM_cols.csv"
#path_to_your_table_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/HKNorm_All_Data_Human_IO_RNA.csv"
#data_label<- "HK"

path_to_your_table_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/NucleiForm_All_Data_Human_IO_RNA.csv"

Code_path <- "/media/rmejia/mountme88/Projects/DSP/Code/DSP-Oszwald/"  
# Path where are the rest of your scripts

path_Results_directory <-"/media/rmejia/mountme88/Projects/DSP/Results/Andre/Nuclei"

data_label<- "NucleiNorm"

colname_4_intra_batch_normalization <- "Scan_ID" # the name of your column to correct
  # please don't use spaces or parenthesis/brackets in the names of your columns
  # The files should have the folowing columns
  # annotation file Morphological_Categories
  
###################################
#### Normalize your paths
###################################
Code_path<-normalizePath(Code_path)
path_Results_directory  <- normalizePath( path_Results_directory  )

###################################
#### Creating your result folders if they doesn't exist yet
###################################
dir.create(path_Results_directory , recursive = TRUE)

###################################
#### Reading the annotation table and the table that cointains the expression data
###################################
table <- read.table( path_to_your_table_file , sep = "\t", header = TRUE )
annot <- read.table( path_to_your_annotation_file , sep = "\t", header = TRUE )

###################################
#### Extracting the Raw expression matrix
###################################
colpositions_withgenes <- grep("OAZ1",colnames(table)):dim(table)[2]
Rawmymatrix <- table[ ,colpositions_withgenes]
Rawmymatrix <- as.matrix(Rawmymatrix)
mode( table[,'ROI_ID']) <- "character" 
rownames( Rawmymatrix ) <- annot[,"Unique_ID"]

#####################
# Annotation object for plotting pcas
####################
annot_4_plotting_pca <- cbind( annot[, c("Morph_cat_Andre","Histology_number","Scan_ID","Biopsy_year","Morphological_Categories")])
annot_4_plotting_pca[,"Histology_number"] <- as.factor(annot_4_plotting_pca[,"Histology_number"])
annot_4_plotting_pca[,"Biopsy_year"] <- as.factor(annot_4_plotting_pca[,"Biopsy_year"])
rownames( annot_4_plotting_pca ) <- rownames( Rawmymatrix )
annot_4_plotting_pca$Morphological_Categories <- as.factor(annot_4_plotting_pca$Morphological_Categories)

# loading the function to melt (reshape) the data to preparation for ggplot2 functions
source( paste0( Code_path,"/","matrix_N_annotdf_2_melteddf.R") )
meltedrawdata <- matrix_N_annotdf_2_melteddf( Rawmymatrix , annot )
head( meltedrawdata )

########################################
########
########    0. Exploratory
########
########################################
#########
### Visualize your the Raw data
########
source(paste0( Code_path,"/","PCA_box_density_plots.R") )
PCA_box_density_plots(  paste0( path_Results_directory,"/Exploratory" )  ,
               Rawmymatrix,  annot_4_plotting_pca , meltedrawdata , paste( data_label, "Data as given" ))

table(annot$Morphological_Categories) # How many categories do you have
########################################
########
########    1. Preprocessing
########
########################################
##################
## log 2 transformation
##################
mymatrix <- log2(Rawmymatrix)+1
meltedlog2data <- matrix_N_annotdf_2_melteddf(mymatrix , annot )

# Visualize the data log2 transformed data
PCA_box_density_plots(  paste0( path_Results_directory,"/Preprocessing" )  ,
                        mymatrix ,  annot_4_plotting_pca , meltedlog2data , paste( data_label, "Log2" ))


###################################
#### Normalization intra batch
###################################
source(paste0( Code_path,"/","Matrix_2_list_of_sub_matrices.R") )
list_of_submatrices <- Matrix_2_list_of_sub_matrices( table$Scan_ID , mymatrix)


################################
###### quantile normalization   , normalizeQuantiles
################################
list_splitd_qnorm <- lapply( list_of_submatrices  ,  normalizeQuantiles)
qnormmat <- do.call(cbind, list_splitd_qnorm)
qnormmat_t <- t(qnormmat)

# All in once q norm
qnorm_all_of_once<- normalizeQuantiles(t(mymatrix))

#### visualize your normalized data

melteqnormmat_t <- matrix_N_annotdf_2_melteddf( qnormmat_t , annot )
PCA_box_density_plots(  paste0( path_Results_directory,"/Preprocessing" )  ,
                        qnormmat_t ,  annot_4_plotting_pca , melteqnormmat_t , paste( data_label, "Log2_Qnorm_intra" ))
melteqqnorm_all_of_once <- matrix_N_annotdf_2_melteddf( t(qnorm_all_of_once) , annot )
PCA_box_density_plots(  paste0( path_Results_directory,"/Preprocessing" )  ,
                        t(qnorm_all_of_once) ,  annot_4_plotting_pca , melteqqnorm_all_of_once , paste( data_label, "Log2_Qnorm_all_inonce" ))


#############
#############

#############
## Looking for batch effect
#############
## Pheno object for combat
pheno <- annot[, c("Morph_cat_Andre","Scan_ID")]
pheno <- annot[, c("Morphological_Categories","Scan_ID")]
colnames(pheno) <- c("subgroups","batch")
rownames(pheno) <- annot[,"Unique_ID"] 

# Batch effect with combat
pheno$batch <- as.factor(pheno$batch)
batch<-pheno$batch
modcombat<-model.matrix(~1, data=pheno)
combat_qnormBsep = ComBat(dat = t(qnormmat_t) , batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
combat_qnormAllinOne = ComBat(dat = qnorm_all_of_once , batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
#modcombat<-model.matrix(~subgroups, data=pheno)
#combat_mydata<-ComBat(dat= t(mymatrix), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

meltedcombat_qsepareted <- matrix_N_annotdf_2_melteddf( t(combat_qnormBsep) , annot )
PCA_box_density_plots(  paste0( path_Results_directory,"/Preprocessing" )  ,
                        t(combat_qnormBsep) ,  annot_4_plotting_pca , meltedcombat_qsepareted , paste( data_label, "Log2_QnormBatchSep_combat" ))

meltedcombat_qnormAllinOne <- matrix_N_annotdf_2_melteddf( t(combat_qnormAllinOne) , annot )
PCA_box_density_plots(  paste0( path_Results_directory,"/Preprocessing" )  ,
                        t(combat_qnormAllinOne) ,  annot_4_plotting_pca , meltedcombat_qnormAllinOne , paste( data_label, "Log2_QnormInBatch_combatAcross" ))



