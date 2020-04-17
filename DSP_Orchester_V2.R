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
if (!require("limma")) {
  BiocManager::install("limma", ask =FALSE)
  library("limma")
}
# ???
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

path_to_your_annotation_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Annotation/Annotation_file_Symplified_Corrected.csv"
#path_to_your_table_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/HKNorm_All_Data_Human_IO_RNA.csv"
#data_label<- "HK"

path_to_your_table_file <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_in_CSV_format/NucleiForm_All_Data_Human_IO_RNA.csv"
data_label<- "NucleiNorm"

path_Results_directory <-"/media/rmejia/mountme88/Projects/DSP/Results/prezoom/Nuclei"

intra_batch_normalization <- # the name of your column to correct
  # please don't use spaces in the names of your columns
  # The files should have the folowing columns
  # annotation file Morphological_Categories
###################################
#### Creating your files
###################################
dir.create(path_Results_directory , recursive = TRUE)

###################################
#### Reading the data
###################################
table <- read.table( path_to_your_table_file , sep = "\t", header = TRUE)
annot <- read.table( path_to_your_annotation_file , sep = "\t", header = TRUE)
mymatrix <- table[ ,10:dim(table)[2]]
mymatrix <- as.matrix(mymatrix)
mode( table[,'ROI_ID']) <- "character" # ? 
rownames( mymatrix ) <- annot[,"Unique_ID"]

#####################
# Annotation object for plotting
####################

# annot_4_pca <- cbind( annot[, c("Morph.cat..Andre.","Histology.number","Scan_ID","Biopsy.year")], pheno )
annot_4_pca <- cbind( annot[, c("Morph.cat..Andre.","Histology.number","Scan_ID","Biopsy.year")])
annot_4_pca[,"Histology.number"] <- as.factor(annot_4_pca[,"Histology.number"])
annot_4_pca[,"Biopsy.year"] <- as.factor(annot_4_pca[,"Biopsy.year"])
rownames( annot_4_pca ) <- rownames(mymatrix)

##   ggplot2

# meltme 4  ggplot2
annot_4_pca
ScanIDEdited <- as.character(table$Scan_ID) 
ScanIDEdited <- gsub("e RNA","_",ScanIDEdited)
ROI_ScanID <- paste0(ScanIDEdited,as.character(table$ROI_ID))
table_with_uniqID <- cbind(ROI_ScanID,table)
colpositions_withgenes <- grep("OAZ1",colnames(table_with_uniqID)):dim(table_with_uniqID)[2]
table_melted_MtxAndScanID <-melt(data=table_with_uniqID, id.vars="ROI_ScanID",measure.vars =  colpositions_withgenes)
table_melted_MtxAndScanID <-melt(data=table_with_uniqID, id.vars="Scan_ID",measure.vars =  colpositions_withgenes)
colnames(table_melted_MtxAndScanID)[2] <- "gene"
ScanID4melted <- as.character( rep(table$Scan_ID, length(colpositions_withgenes)) )
ROI4melted <- as.character( rep(table$ROI_ID, length(colpositions_withgenes)) )
#PathoDescrip4melted <- as.character( rep(table$, length(colpositions_withgenes)) )
table_melted_MtxAndScanID <- cbind(table_melted_MtxAndScanID , ScanID4melted)


head(table)
head(table_melted_MtxAndScanID, 100)
q <- ggplot(table_melted_MtxAndScanID, aes(ROI_ScanID, value))
q + geom_boxplot()
annot

#########
###
########



dir.create(paste0(path_Results_directory,"/Exploratory"), recursive = TRUE)
pdf( file=paste0(path_Results_directory,"/Exploratory","/" ,data_label,"_Exploring_Data_as_given.pdf"),
     width = 10, height = 7)
autoplot( prcomp( mymatrix ), data = annot_4_pca, colour= 'Scan_ID') +
  ggtitle(paste(data_label,"Exploring Data as given"))
autoplot( prcomp( mymatrix ), data = annot_4_pca, colour= 'Histology.number', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Exploring Data as given"))
autoplot( prcomp( mymatrix ), data = annot_4_pca, colour= 'Biopsy.year', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Exploring Data as given"))
autoplot( prcomp( mymatrix ), data = annot_4_pca, colour= 'Morph.cat..Andre.', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Exploring Data as given"))
boxplot(t(mymatrix), col=annot_4_pca$Scan_ID, main="Scan ID")
boxplot(t(mymatrix) , col=annot_4_pca$Histology.number, main="Histology number")
boxplot( t(mymatrix) , col=annot_4_pca$Biopsy.year, main="Biopsy year")
boxplot( t(mymatrix), col=annot_4_pca$Morph.cat..Andre. , main="Morph Cat André"  )
dev.off()
dim(mymatrix)
table(annot_4_pca$Morph.cat..Andre.)
##################
## log 2 transformation
##################
mymatrix <- log2(mymatrix)+1

# Exploring Raw Data with log2 trandformation
dir.create(paste0(path_Results_directory,"/Exploratory"), recursive = TRUE)
pdf( file=paste0(path_Results_directory,"/Exploratory","/" ,data_label,"Raw_Data_log2_tranformed.pdf"),
     width = 10, height = 7)
autoplot( prcomp( mymatrix ), data = annot_4_pca, colour= 'Scan_ID') +
  ggtitle(paste(data_label,"Raw Data log 2 tranformed"))
autoplot( prcomp( mymatrix ), data = annot_4_pca, colour= 'Histology.number', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Raw Data log 2 tranformed"))
autoplot( prcomp( mymatrix ), data = annot_4_pca, colour= 'Biopsy.year', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Raw Data log 2 tranformed"))
autoplot( prcomp( mymatrix ), data = annot_4_pca, colour= 'Morph.cat..Andre.', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Raw Data log 2 tranformed"))
boxplot(t(mymatrix), col=annot_4_pca$Scan_ID, main="Scan ID")
boxplot(t(mymatrix) , col=annot_4_pca$Histology.number, main="Histology number")
boxplot( t(mymatrix) , col=annot_4_pca$Biopsy.year, main="Biopsy year")
boxplot( t(mymatrix), col=annot_4_pca$Morph.cat..Andre. , main="Morph Cat André"  )
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
uniqbatch_colors <- brewer.pal(n = length(unique(batches)), name = "Set1")
sapply( uniqbatch_colors , color.id)
humanreadablecolors<- as.character(sapply( uniqbatch_colors , color.id))
colbatchcolors <- rep( uniqbatch_colors  , times = as.vector(table(batches)))

plotDensities(data, col=colbatchcolors, main="Scan ID", legend= FALSE)
plotDensities(data, col=rep( uniqbatch_colors  , times = as.vector(table(batches))), main="Scan ID")
boxplot(data, col=colbatchcolors, main="Scan ID")


ScanIDEdited <- as.character(table$Scan_ID) 
ScanIDEdited <- gsub("e RNA","_",ScanIDEdited)
ROI_ScanID <- paste0(ScanIDEdited,as.character(table$ROI_ID))
table_with_uniqID <- cbind(ROI_ScanID,table)
colpositions_withgenes <- grep("OAZ1",colnames(table_with_uniqID)):dim(table_with_uniqID)[2]
table_melted_MtxAndScanID <-melt(data=table_with_uniqID, id.vars="ROI_ScanID",measure.vars =  colpositions_withgenes)
table_melted_MtxAndScanID <-melt(data=table_with_uniqID, id.vars="Scan_ID",measure.vars =  colpositions_withgenes)
colnames(table_melted_MtxAndScanID)[2] <- "gene"
ScanID4melted <- as.character( rep(table$Scan_ID, length(colpositions_withgenes)) )
ROI4melted <- as.character( rep(table$ROI_ID, length(colpositions_withgenes)) )
#PathoDescrip4melted <- as.character( rep(table$, length(colpositions_withgenes)) )
table_melted_MtxAndScanID <- cbind(table_melted_MtxAndScanID , ScanID4melted)


head(table)
head(table_melted_MtxAndScanID, 100)
q <- ggplot(table_melted_MtxAndScanID, aes(ROI_ScanID, value))
q + geom_boxplot()

head(mpg, n= 30)
p <- ggplot(mpg, aes(class, hwy))
p + geom_boxplot()





head(table)

head(mpg, n= 20)
p <- ggplot(mpg, aes(class, hwy))
p + geom_boxplot()
table[1:5,1:5]
1:5; 
letters[1:5]
melt()
data(mpg)
class(mpg)
dim(mpg)
head(mpg, n= 20)
p <- ggplot(mpg, aes(class, hwy))
p + geom_boxplot()

plotDensities(quantile_normalisation(df), col=c("red","red","green"),legend=FALSE)
plotDensities(df, col=c("red","red","green"),legend=FALSE)

matrix_list <- split(mymatrix,annot_4_pca$Scan_ID )

matrix_list
boxplot(matrix_list[[1]])
class(matrix_list[[1]])
length(matrix_list[[1]])
matrix_list[[1]] <- matrix(matrix_list[[1]] , nrow=72)
?as.matrix
dim(matrix_list[[1]])

# Normalization Between slides ~ Micoarrays (It is pertinent?) 
data <- t(mymatrix)
plotDensities(t(data))
batches <- annot_4_pca$Scan_ID



plotDensities(t(data), col=rainbow(length(annot_4_pca$Scan_ID)),legend=FALSE)
#boxplot(log2(data)+1)
boxplot(data, col=annot_4_pca$Scan_ID, main="Scan ID")
boxplot(data, col=annot_4_pca$Histology.number, main="Histology number")
boxplot(data, col=annot_4_pca$Biopsy.year, main="Biopsy year")
boxplot(data, col=annot_4_pca$Morph.cat..Andre. , main="Morph Cat André"  )
geom_boxplot(aes=data)
data(mpg)
p <- ggplot(mpg, aes(class, hwy))
p + geom_boxplot()


plotDensities(data[1:84,1:71])

rownames(data)
dim(data)
str(data)
boxplot(log2(data)+1)
plotDensities( t(log2(data)+1))
plotDensities(matrix(rnorm(100,0,1), ncol=10))
?plotDensities
head(data)
str(data)
rownames(data) <- data[,1]
data_mat <- data.matrix(data[,-1]) 
head(data_mat)
data_norm <- normalize.quantiles(data_mat, copy = TRUE)


#############
## Looking for batch effect
#############
## Pheno object for combat
pheno <- annot[, c("Morph.cat..Andre.","Histology.number")]
pheno <- annot[, c("Morph.cat..Andre.","Scan_ID")]
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
autoplot( prcomp( t(combat_mydata) ), data = annot_4_pca, colour= 'Scan_ID', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Post Combat"))
autoplot( prcomp( t(combat_mydata) ), data = annot_4_pca, colour= 'Histology.number', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Post Combat"))
autoplot( prcomp( t(combat_mydata) ), data = annot_4_pca, colour= 'Biopsy.year', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Post Combat"))
autoplot( prcomp( t(combat_mydata) ), data = annot_4_pca, colour= 'Morph.cat..Andre.', label = TRUE, label.size = 3) +
  ggtitle(paste(data_label,"Post Combat"))
boxplot(combat_mydata, col=annot_4_pca$Scan_ID, main="Scan ID")
boxplot(combat_mydata , col=annot_4_pca$Histology.number, main="Histology number")
boxplot( combat_mydata , col=annot_4_pca$Biopsy.year, main="Biopsy year")
boxplot( combat_mydata, col=annot_4_pca$Morph.cat..Andre. , main="Morph Cat André"  )
dev.off()




library(limma)
library(sma)
help.start()
data(MouseArray)
# just bkg
corr
MA.n <- MA.RG(mouse.data)                               # no
normalization

MA.q <- normalizeBetweenArrays(MA.p, method = "q")      # quantile
norm

G.q <- normalizeBetweenArraysRG.MA(MA.n)$G,method="q") # only green
sc's
                                                        # takes a
matrix

tmp<-cbindRG.MA(MA.n)$R,RG.MA(MA.n)$G)[,c(1,3,8,9,12)] # select sc's
tmp.q <- normalizeBetweenArrays(tmp, method = "q")      # takes a
matrix

MA.p <- normalizeWithinArrays(MA.n, mouse.setup)        # default
p-loess
MA.pq <- normalizeBetweenArrays(MA.p, method = "q")     # pq norm

MA.MpAq <- normalizeBetweenArrays(MA.p, method = "Aq")  # MpAq norm
# Yang &
Thorne 03


plotDensities(MA.n)                                     # default

plotDensities(MA.n,arrays=c(1:6),                       # same as
              default
              groups=c(rep(1,6),rep(2,6)),col=c("red","green"))

plotDensities(MA.n,arrays=NULL,groups=NULL,             # diff cols
              col=c("blue","purple"))

