###################################
#### This script receives a expression matrix, an annotation table, and a feature from your exp mat which you want to see in a box plot.
####    The program retrieves a box plot of the values of the feature according each category from your AnnotTable
####     
####    https://github.com/raulmejia/DSP-Oszwald
####      
####     Input description:
####        expression matrix: tab separated, rows = featrures, cols= sources/samples, there should not be colname over the first column(for the rownames)
####        annotdf: rows = cols from the expression matrix, you should provide the name of column that contain your labels
####                Try to use character labels and if you have a "reference class" you can give that name as well for example "Control" 
####                against this category will be arranged
####     Output:
####        Boxplot of your feaute's data over the categories
####
####    Author of the script: Raúl Mejía
#### 
###################################
#### 0) loading and/or installing required libraries
###################################
if (!require("argparse")) {
  install.packages("argparse", ask =FALSE)
  library("argparse")
}
if (!require("gtools")) {
  install.packages("gtools", ask =FALSE)
  library("gtools")
}
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("limma")) {
  BiocManager::install("limma", ask =FALSE)
  library("limma")
}
if (!require("tidyverse")) {
  install.packages("tidyverse", ask =FALSE)
  library("tidyverse")
}
########################################
#### data given by the user
#########################################
myargs <- commandArgs(trailingOnly = TRUE)

# path_expression_matrix <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_sent_by_DrAO_2020_12_11/Protein/Prunned/ExpMat_from_the_fitting_of--NucleiLog2_QnormBatchSep_combat_Prot_titleclean--n--DSP_protein_annotation_Raul_20201211_odd_deleted_characters_instead_numbers.tsv"
path_expression_matrix <- myargs[1]

# path_annotation_table <-"/media/rmejia/mountme88/Projects/DSP/Data/Data_sent_by_DrAO_2020_12_11/Protein/Prunned/Annot_table_from_the_fitting_of--NucleiLog2_QnormBatchSep_combat_Prot_titleclean--n--DSP_protein_annotation_Raul_20201211_odd_deleted_characters_instead_numbers.tsv"
path_annotation_table <- myargs[2]

# column_of_labels_in_the_annotdf <- "label" # name of the column in your annotation data frame that contains the labels
column_of_labels_in_the_annotdf <- myargs[3]

# your_interesting_feature <- "CD44"
your_interesting_feature <- myargs[4]

# Folder_to_save_results <- "/media/rmejia/mountme88/Projects/DSP/Results/DEG/Protein/"
Folder_to_save_results <- myargs[5]

# label_for_your_results <- "Box_plot_-ExpMat-Annot-Feature-.pdf"
# label_for_your_results <- "Box_plot--DSPGeoMx_Protein---Annot--CD44"
label_for_your_results <- myargs[6]

#your_title <- "Protein Digital Spatial Profiling  GEOmx"
your_title <- myargs[7]

#your_x_label <- "Diagnoses"
your_x_label <- myargs[8]

############################
###### Body
#############################
Folder_to_save_results <- normalizePath(Folder_to_save_results)
# reading your matrices
myexpmat <- read.table( file= path_expression_matrix , sep="\t", header = TRUE, row.names = 1, check.names = FALSE)
annotdf <- read.table( file= path_annotation_table, sep="\t", header = TRUE, row.names = 1)

if( length(which( colnames(myexpmat)  %in% rownames(annotdf) )) != length(colnames(myexpmat))  ){
  print("The number columns in your expmat is different than the rows in your Annotation Table")
  quit()
}
if( length(which( colnames(myexpmat)  %in% rownames(annotdf) )) != length( rownames(annotdf) )  ){
  print("The number columns in your expmat is different than the rows in your Annotation Table")
  quit()
}

obj_4_gplot <- as.data.frame( t( myexpmat[your_interesting_feature ,] ) )
obj_4_gplot$categories <- annotdf[  , column_of_labels_in_the_annotdf ]
obj_4_gplot$categories <- as.factor( obj_4_gplot$categories )
categories1 <- "categories"

path_2pdf <-  paste0(Folder_to_save_results,"/",label_for_your_results,".pdf" )
pdf(file=path_2pdf, width=8, height = 8 )
    q <- ggplot( obj_4_gplot, aes_string(x= categories1 , y = your_interesting_feature, color = categories1 ) ) +  geom_boxplot()
    q + labs(title= your_title ,x= your_x_label , y =your_interesting_feature )+ geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3) +   theme_classic()
dev.off()

