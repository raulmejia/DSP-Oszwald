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

path_Results_directory <-"/media/rmejia/mountme88/Projects/DSP/Results/prezoom/Nuclei"

data_label<- "NucleiNorm"

colname_4_intra_batch_normalization <- "Scan_ID" # the name of your column to correct
  # please don't use spaces or parenthesis/brackets in the names of your columns
  # The files should have the folowing columns
  # annotation file Morphological_Categories
  


###################################
#### Creating your result folders if they doesn't exist yet
###################################
dir.create(path_Results_directory , recursive = TRUE)

###################################
#### Reading the annotation table and the table that cointains the expression data
###################################
table <- read.table( path_to_your_table_file , sep = "\t", header = TRUE)
annot <- read.table( path_to_your_annotation_file , sep = "\t", header = TRUE)

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
source( paste0( Code_path,"matrix_N_annotdf_2_melteddf.R") )
meltedrawdata <- matrix_N_annotdf_2_melteddf( Rawmymatrix , annot )
head( meltedrawdata )

#########
### Generating plots to describe the Raw data
########


dir.create(paste0(path_Results_directory,"/Exploratory"), recursive = TRUE)
pdf( file=paste0(path_Results_directory,"/Exploratory","/" ,data_label,"_Exploring_Data_as_given.pdf"),
     width = 10, height = 7)
autoplot( prcomp( Rawmymatrix ), data = annot_4_plotting_pca, colour= 'Scan_ID') +
  ggtitle(paste(data_label,"Exploring Data as given"))
autoplot( prcomp( Rawmymatrix ), data = annot_4_plotting_pca, colour= 'Histology_number', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Exploring Data as given"))
autoplot( prcomp( Rawmymatrix ), data = annot_4_plotting_pca, colour= 'Biopsy_year', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Exploring Data as given"))
autoplot( prcomp( Rawmymatrix ), data = annot_4_plotting_pca, colour= 'Morphological_Categories', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Exploring Data as given"))
q <- ggplot( meltedrawdata, aes(Unique_ID, value, fill=Scan_ID))
q + geom_boxplot( )+ ggtitle(paste(data_label,"Data as given"))
q <- ggplot( meltedrawdata, aes(Unique_ID, value, fill=Histology_number))
q + geom_boxplot( )+ ggtitle(paste(data_label,"Data as given"))
q <- ggplot( meltedrawdata, aes(Unique_ID, value, fill=Biopsy_year))
q + geom_boxplot( )+ ggtitle(paste(data_label,"Data as given"))
q <- ggplot( meltedrawdata, aes(Unique_ID, value, fill= Morphological_Categories))
q + geom_boxplot( )+ ggtitle(paste(data_label,"Data as given"))
dev.off()

table(annot$Morphological_Categories)

##################
## log 2 transformation
##################
mymatrix <- log2(Rawmymatrix)+1
meltedrawdata <- matrix_N_annotdf_2_melteddf(mymatrix,annot )
# Exploring Raw Data with log2 trandformation
dir.create(paste0(path_Results_directory,"/Exploratory"), recursive = TRUE)
pdf( file=paste0(path_Results_directory,"/Exploratory","/" ,data_label,"Raw_Data_log2_tranformed.pdf"),
     width = 10, height = 7)
autoplot( prcomp( mymatrix ), data = annot_4_plotting_pca, colour= 'Scan_ID') +
  ggtitle(paste(data_label,"Raw Data log 2 tranformed"))
autoplot( prcomp( mymatrix ), data = annot_4_plotting_pca, colour= 'Histology_number', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Raw Data log 2 tranformed"))
autoplot( prcomp( mymatrix ), data = annot_4_plotting_pca, colour= 'Biopsy_year', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Raw Data log 2 tranformed"))
autoplot( prcomp( mymatrix ), data = annot_4_plotting_pca, colour= 'Morphological_Categories', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Raw Data log 2 tranformed"))
q <- ggplot( meltedrawdata, aes(Unique_ID, value, fill=Scan_ID))
q + geom_boxplot( )+ ggtitle(paste(data_label,"Raw Data log 2 tranformed"))
q <- ggplot( meltedrawdata, aes(Unique_ID, value, fill=Histology_number))
q + geom_boxplot( )+ ggtitle(paste(data_label,"Raw Data log 2 tranformed"))
q <- ggplot( meltedrawdata, aes(Unique_ID, value, fill=Biopsy_year))
q + geom_boxplot( )+ ggtitle(paste(data_label,"Raw Data log 2 tranformed"))
q <- ggplot( meltedrawdata, aes(Unique_ID, value, fill= Morphological_Categories))
q + geom_boxplot( )+ ggtitle(paste(data_label,"Raw Data log 2 tranformed"))
dev.off()

###################################
#### Normalization between batch
###################################
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

#test the function
df <- data.frame(one=c(5,2,3,4),
                 two=c(4,1,4,2),
                 three=c(3,4,6,8)
)
rownames(df) <- toupper(letters[1:4])
#uniqbatch_colors <- brewer.pal(n = length(unique(batches)), name = "Dark2")
#uniqbatch_colors <- brewer.pal(n = length(unique(batches)), name = "Set3")
annot[ , colname_4_intra_batch_normalization]
uniqbatch_colors <- brewer.pal( n = length( unique( annot[ , colname_4_intra_batch_normalization] )), name = "Set1")
sapply( uniqbatch_colors , color.id) ; humanreadablecolors<- as.character(sapply( uniqbatch_colors , color.id))
colbatchcolors <- rep( uniqbatch_colors  , times = as.vector(table( annot[ , colname_4_intra_batch_normalization])))

plotDensities(Rawmymatrix, col=colbatchcolors, main="Scan ID", legend= FALSE)
plotDensities(  mymatrix, col=colbatchcolors, main="Scan ID", legend= FALSE)



ScanIDEdited <- as.character(table$Scan_ID) 
ScanIDEdited <- gsub("e RNA","_",ScanIDEdited)
ROI_ScanID <- paste0(ScanIDEdited,as.character(table$ROI_ID))
table_with_uniqID <- cbind(ROI_ScanID,table)
colpositions_withgenes <- grep("OAZ1",colnames(table_with_uniqID)):dim(table_with_uniqID)[2]
table_melted_MtxAndScanID <-melt(data=table_with_uniqID, id.vars="ROI_ScanID",measure.vars =  colpositions_withgenes)
#table_melted_MtxAndScanID <-melt(data=table_with_uniqID, id.vars="Scan_ID",measure.vars =  colpositions_withgenes)
colnames(table_melted_MtxAndScanID)[2] <- "gene"
ScanID4melted <- as.character( rep(table$Scan_ID, length(colpositions_withgenes)) )
ROI4melted <- as.character( rep(table$ROI_ID, length(colpositions_withgenes)) )
#PathoDescrip4melted <- as.character( rep(table$, length(colpositions_withgenes)) )
table_melted_MtxAndScanID <- cbind(table_melted_MtxAndScanID , ScanID4melted)

head(table)
head(table_melted_MtxAndScanID, 100)

plot3 <- ggplot(table_melted_MtxAndScanID, aes(x=value, fill = ScanID4melted, y = ROI_ScanID)) + 
  geom_density_ridges() + 
  scale_x_continuous(trans = "log10") +
  ggtitle("Density Scaled to 1")
plot3

plot4 <- ggplot(table_melted_MtxAndScanID, aes(x=value, color = ScanID4melted)) + 
  geom_density() + 
  scale_x_continuous(trans = "log2")+
  ggtitle("Height Scaled to X")

plot_grid(plot3, plot4, nrow = 1)
plot_grid(plot3, plot4, plot4, plot3, nrow = 2)

################################
###### quantile normalization 
################################
matrix4intrabatchnorm <- cbind(as.character(table$Scan_ID),as.data.frame(mymatrix))
colnames(matrix4intrabatchnorm) <- "batch"
matrix4intrabatchnorm[,"batch"] <- as.factor( matrix4intrabatchnorm[,"batch"] )
matrix4intrabatchnorm[1:4,1:4]

list_splited_table <- split(matrix4intrabatchnorm, matrix4intrabatchnorm$batch)
bye1stcol <- function(df){
  print(dim(df))
  df<- df[,-1]
  print(dim(df))
  return(df)
} 
list_splited_table_no1stcol <- lapply(list_splited_table, bye1stcol)
putcolnmaes <- function(df){
  colnames(df) <- colnames(mymatrix)
  return(df)
}
list_splited_table_propercolnames <- lapply(list_splited_table_no1stcol, putcolnmaes)
list_splited_table_propercolnames[[1]][1:4,1:4]

list_splited_table_t <- lapply(list_splited_table_propercolnames, t )
makeme_num_mat <-function(df){
  df <- as.matrix(df)
  mode(df) <- "numeric"
  return(df)
}
list_splited_nummat <- lapply(list_splited_table_t, makeme_num_mat)
list_splitd_qnorm <- lapply(list_splited_nummat, quantile_normalisation)
qnormmat <- do.call(cbind, list_splitd_qnorm)
dim(qnormmat)
qnormmat_t <- t(qnormmat)
dim(qnormmat_t)
qnormmat_t[12:15, 1:5]
mymatrix[12:15,1:5]

### Achtung!! all normalization
#qnormmat_t<-quantile_normalisation(t(mymatrix))
qnormmat_t[1:4,1:4]
t(mymatrix)[1:4,1:4]
plotDensities(t(qnormmat_t), legend = FALSE)
plotDensities(t(qnormmat_t))

#############
#### visualize your normalized data
##############
qnormmat_t[1:4,1:4]
head(table)

colpositions_withoutgenes <- 1:(grep("OAZ1",colnames(table))-1)
table_qnorm <- cbind(table[,colpositions_withoutgenes], qnormmat_t)

ScanIDEdited_qnorm <- as.character(table_qnorm$Scan_ID) ; ScanIDEdited_qnorm <- gsub("e RNA","_",ScanIDEdited_qnorm)
ROI_ScanID_qnorm <- paste0(ScanIDEdited_qnorm,as.character(table_qnorm$ROI_ID))
table_with_uniqID_qnorm <- cbind(ROI_ScanID_qnorm,table_qnorm)
colpositions_withgenes_qnorm <- grep("OAZ1",colnames(table_with_uniqID_qnorm)):dim(table_with_uniqID_qnorm)[2]
table_melted_MtxAndScanID_qnorm <-melt(data=table_with_uniqID_qnorm, id.vars="ROI_ScanID_qnorm",measure.vars =  colpositions_withgenes_qnorm)
#table_melted_MtxAndScanID <-melt(data=table_with_uniqID, id.vars="Scan_ID",measure.vars =  colpositions_withgenes)
colnames(table_melted_MtxAndScanID_qnorm)[2] <- "gene"
ScanID4melted_qnorm <- as.character( rep(table_qnorm$Scan_ID, length(colpositions_withgenes)) )
ROI4melted_qnorm <- as.character( rep(table_qnorm$ROI_ID, length(colpositions_withgenes)) )
#PathoDescrip4melted <- as.character( rep(table$, length(colpositions_withgenes)) )
table_melted_MtxAndScanID_qnorm <- cbind(table_melted_MtxAndScanID_qnorm , ScanID4melted_qnorm)

head(table_qnorm)
head(table_melted_MtxAndScanID_qnorm, 100)

plot3 <- ggplot(table_melted_MtxAndScanID_qnorm, aes(x=value, fill = ScanID4melted_qnorm, y = ROI_ScanID_qnorm)) + 
  geom_density_ridges() 
  #scale_x_continuous(trans = "log10") +
  #ggtitle("Density Scaled to 1")
plot3

plot4 <- ggplot(table_melted_MtxAndScanID_qnorm, aes(x=value, color = ScanID4melted_qnorm)) + 
  geom_density() 
  #scale_x_continuous(trans = "log2")+
  #ggtitle("Height Scaled to X")

plot_grid(plot3, plot4, nrow = 1)


##### normalize all togheter qnorm



#############
## Looking for batch effect
#############
## Pheno object for combat
#pheno <- annot[, c("Morphological_Categories","Histology_number")]
pheno <- annot[, c("Morph_cat_Andre","Scan_ID")]
pheno <- annot[, c("Morphological_Categories","Scan_ID")]
colnames(pheno) <- c("subgroups","batch")
rownames(pheno) <- annot[,"Unique_ID"] 

# Batch effect with combat
pheno$batch <- as.factor(pheno$batch)
batch<-pheno$batch
modcombat<-model.matrix(~1, data=pheno)
combat_mydata= ComBat(dat = t(qnormmat_t) , batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
t(qnormmat_t)
#modcombat<-model.matrix(~subgroups, data=pheno)
#combat_mydata<-ComBat(dat= t(mymatrix), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)


# After combat
pdf( file=paste0(path_Results_directory,"/Exploratory","/" ,data_label,"_Postcombat.pdf"),
     width = 10, height = 7)
autoplot( prcomp( t(combat_mydata) ), data = annot_4_plotting_pca, colour= 'Scan_ID', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Post Combat"))
autoplot( prcomp( t(combat_mydata) ), data = annot_4_plotting_pca, colour= 'Histology_number', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Post Combat"))
autoplot( prcomp( t(combat_mydata) ), data = annot_4_plotting_pca, colour= 'Biopsy_year', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Post Combat"))
autoplot( prcomp( t(combat_mydata) ), data = annot_4_plotting_pca, colour= 'Morphological_Categories', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Post Combat"))
boxplot(combat_mydata, col=annot_4_plotting_pca$Scan_ID, main="Scan ID")
boxplot(combat_mydata , col=annot_4_plotting_pca$Histology_number, main="Histology number")
boxplot( combat_mydata , col=annot_4_plotting_pca$Biopsy_year, main="Biopsy year")
boxplot( combat_mydata, col=annot_4_plotting_pca$Morphological_Categories , main="Morph Cat André"  )
dev.off()

pheatmap(combat_mydata, annotation= annot)
pheatmap(combat_mydata, annotation= annot_4_plotting_pca)

?pheatmap
#############
## Looking for batch effect
#############
## Pheno object for combat
#pheno <- annot[, c("Morphological_Categories","Histology_number")]
pheno <- annot[, c("Morph_cat_Andre","Scan_ID")]
pheno <- annot[, c("Morphological_Categories","Scan_ID")]
colnames(pheno) <- c("subgroups","batch")
rownames(pheno) <- annot[,"Unique_ID"] 

# Batch effect with combat
pheno$batch <- as.factor(pheno$batch)
batch<-pheno$batch
modcombat<-model.matrix(~1, data=pheno)
combat_mydata= ComBat(dat = t(mymatrix) , batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

#modcombat<-model.matrix(~subgroups, data=pheno)
#combat_mydata<-ComBat(dat= t(mymatrix), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)


# After combat
pdf( file=paste0(path_Results_directory,"/Exploratory","/" ,data_label,"_Postcombat.pdf"),
     width = 10, height = 7)
autoplot( prcomp( t(combat_mydata)  ), data = annot_4_plotting_pca, colour= 'Scan_ID', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Post Combat"))
autoplot( prcomp( t(combat_mydata) ), data = annot_4_plotting_pca, colour= 'Histology_number', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Post Combat"))
autoplot( prcomp( t(combat_mydata) ), data = annot_4_plotting_pca, colour= 'Biopsy.year', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Post Combat"))
autoplot( prcomp( t(combat_mydata) ), data = annot_4_plotting_pca, colour= 'Morphological_Categories', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Post Combat"))
boxplot(combat_mydata, col=annot_4_plotting_pca$Scan_ID, main="Scan ID")
boxplot(combat_mydata , col=annot_4_plotting_pca$Histology.number, main="Histology number")
boxplot( combat_mydata , col=annot_4_plotting_pca$Biopsy.year, main="Biopsy year")
boxplot( combat_mydata, col=annot_4_plotting_pca$Morph.cat..Andre. , main="Morph Cat André"  )
dev.off()


