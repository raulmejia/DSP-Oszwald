###################################
#### This script performs the sample´s intersection between an expression matrix and its annotation table. 
####    Inputs: an expression matrix and an annotation table. both in a tsv format (tab separated).
####         Input description:
####             expression matrix: tab separated, rows = featrures, cols= sources/samples, there should not be colname over the first column(for the rownames)
####                Please try to avoid X in or number as the first character of your columns names, that doens´t like to R              
####              annotdf: rows = cols from the expression matrix
####    Output: The program retrieves a subexpression matrix that only contains the common samples (in the columns)
####            likewise for samples (rows) in the the sub annotation matrix. the program shows the dimensions of 
####            the input and output data on the screen and save it in a file
####    Example: Rscript /pathA/Pruning_to_Match_ExpMat_n_AnnotTable.R  /pathB/Input-myExpMAt /pathC/Input-AnnotTable /pathD/path_to_store_results /pathE/2Output-ExpMat /pathF/2Output-AnnotTable /pathG/2OutputStatistics
####
####    Repository: https://github.com/raulmejia/DSP-Oszwald
####      
####    Author of the script: Raúl Mejía
####  Annot_table_
####  path_annotation_table <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_sent_by_DrAO_2020_12_11/RNA/Final_annot_RNA_20201212_odd_deleted.tsv"
####      Maybe you can produce the outpus under this formar = "ExpMat_from_the_fitting_of--",basename?InputExpMat,","--n--",basename?InputDF,".tsv"
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
if (!require("tidyverse")) {
  install.packages("tidyverse", ask =FALSE)
  library("tidyverse")
}
########################################
#### data given by the user
########################################
myargs <- commandArgs(trailingOnly = TRUE)

# path_expression_matrix <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_sent_by_DrAO_2020_12_11/Protein/NucleiLog2_QnormBatchSep_combat_Prot_titleclean.txt"
# path_expression_matrix <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_sent_by_DrAO_2020_12_11/RNA/NucleiLog2_QnormBatchSep_combat_2.txt"
path_expression_matrix <- myargs[1]

# path_annotation_table <-"/media/rmejia/mountme88/Projects/DSP/Data/Data_sent_by_DrAO_2020_12_11/Protein/DSP_protein_annotation_Raul_20201211_odd_deleted.tsv"
# path_annotation_table <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_sent_by_DrAO_2020_12_11/RNA/Final_annot_RNA_20201212_odd_deleted_characters_instead_numbers.tsv"
path_annotation_table <- myargs[2]

# Folder_to_save_results <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_sent_by_DrAO_2020_12_11/Protein/"
# Folder_to_save_results <- "/media/rmejia/mountme88/Projects/DSP/Data/Data_sent_by_DrAO_2020_12_11/RNA/Prunned"
Folder_to_save_results <- myargs[3]

# label_for_the_result_expMat <- "ExpMat_from_the_fitting_of_NucleiLog2_QnormBatchSep_combat_Prot_titleclean--n--protein_annotation_Raul_20201211_odd_deleted.tsv"
# label_for_the_result_expMat <- "ExpMat_from_the_fitting_of--RNANucleiLog2_QnormBatchSep_combat_Prot_titleclean--n--Final_annot_RNA_20201212_odd_deleted_characters_instead_numbers.tsv"
# Usually pasting the name of the expmat and annotTable is a good idea or adding the suffix "intersected"/"prunned to"/"fitted"
label_for_the_result_expMat <- myargs[4]

# label_for_the_result_AnnotTable <- "Annot_table_from_the_fitting_of_NucleiLog2_QnormBatchSep_combat_Prot_titleclean--n--fitted2-DSP_protein_annotation_Raul_20201211_odd_deleted.tsv"
# label_for_the_result_AnnotTable <-  "Annot_table_from_the_fitting_of--RNANucleiLog2_QnormBatchSep_combat_Prot_titleclean--n--Final_annot_RNA_20201212_odd_deleted_characters_instead_numbers.tsv"
# Usually pasting the name of the expmat and annotTable is a good idea or adding the suffix "intersected"/"prunned to"/"fitted"
label_for_the_result_AnnotTable <- myargs[5]

# label_for_the_output_statistics <- "Output_statistics_from_the_fitting_of_NucleiLog2_QnormBatchSep_combat_Prot_titleclean--n--fitted2-DSP_protein_annotation_Raul_20201211_odd_deleted.txt"
# label_for_the_output_statistics <- "Output_statistics_from_the_fitting_of_RNANucleiLog2_QnormBatchSep_combat_Prot_titleclean--n--fitted2-Final_annot_RNA_20201212_odd_deleted_characters_instead_numbers.txt"
label_for_the_output_statistics <- myargs[6]

############################
###### Basic Processing of data given by the user for example Paths
############################
dir.create( Folder_to_save_results , recursive = TRUE )
Folder_to_save_results <- normalizePath( Folder_to_save_results ) 

############################
###### Body
#############################
# reading your matrices

input_expmat <- read.table( file = path_expression_matrix , sep="\t", header = TRUE , row.names = 1, check.names = FALSE )
input_annotdf <- read.table( file = path_annotation_table , sep="\t", header = TRUE , row.names = 1 )

colnames_expmat_nomarlized <- colnames( input_expmat )
#colnames_expmat_nomarlized <- str_remove( colnames( input_expmat ) , "^X" ) # Deleting the X at the beginnig of the colnames ( If
# a column´s name start with a number R append automatically an X at the beginning of the name )

positions_inputexpmat_in_annot <- which( colnames_expmat_nomarlized %in% rownames( input_annotdf ) )
colnames( input_expmat[ , positions_inputexpmat_in_annot ] )
expmat_intersected_columns <- input_expmat[ , positions_inputexpmat_in_annot ] 

positions_input_annotdf_in_expmat <- which( rownames(input_annotdf ) %in% colnames_expmat_nomarlized )
anntdf_intersected_rows <- input_annotdf[ positions_input_annotdf_in_expmat , ]

################################## 
## writing down the results
##################################
# Output expression file
out_expmat_path <- paste0(  Folder_to_save_results, "/", label_for_the_result_expMat)
write.table(  expmat_intersected_columns , file = out_expmat_path ,  sep= "\t", quote=FALSE)

# Output Annotation file
out_AnnotTable_path <- paste0(  Folder_to_save_results, "/", label_for_the_result_AnnotTable)
write.table(  anntdf_intersected_rows , file = out_AnnotTable_path ,  sep= "\t", quote=FALSE)

# Output Statistcs:  showing the original and final sizes as well as the differences
OutStatistics_path <- paste0(  Folder_to_save_results, "/", label_for_the_output_statistics)
sink( OutStatistics_path  )
cat( "dimension of your input expression matrix: " )
cat( dim( input_expmat ) )
cat( "\n" )
cat( "dimension of your output expression matrix: " )
cat( dim( expmat_intersected_columns ) )
cat( "\n" )
cat( c("This are the column names included in the input expression matrix but not in the input annotation table:") )
cat( "\n" )
cat( setdiff( colnames_expmat_nomarlized ,  rownames(  input_annotdf )) )
cat( "\n" )
cat( "\n" )
cat( "dimension of your input annotation Table: " )
cat( dim( input_annotdf ) )
cat( "\n" )
cat( "dimension of your output annotation Table: " )
cat( dim( anntdf_intersected_rows  ) )
cat( "\n" )
cat( "The following are the rows mentioned in the input annotation table but not in the columns of the input expression matrix:" )
cat( "\n" )
cat( setdiff( rownames(  input_annotdf ) , colnames_expmat_nomarlized  ) )
cat( "\n" )
sink()
