##################################################################
### This program receives a design matrix as input and retrieves a contrasts' matrix 
### with all the pairwise possible contrasts
### The design matrix should be created like the output of stats::model.matrix
######  Note: Maybe you can improve it asking for a label called "Control" 
######  If one of the levels/labels of the design matrix has the name reference/control,
######  first all the comparisons will be made against it and then the other pairwise comparisons
######  will follow an alphabetical order.
##################################################################


########################################
##### required loading required packages 
########################################
if (!require("gtools")) {
  install.packages("gtools", ask =FALSE)
  library("gtools")
}

##########################
##### Body
##########################
mymakeContrasts <- function( design_input ){ 
  n = length( colnames(design_input) )
  ncols_contrast_mat = (n * (n-1)) /2  # total number of pairwise comparisons
  
  # creating the scaffold of the contrasts' matrix  
  required_zero_entries_for_contrast_mat <-  rep(0, ncols_contrast_mat * dim(design_input)[2]  )
  result_Contrasts_matrix <- matrix( required_zero_entries_for_contrast_mat , nrow = dim(design_input)[2] )  
  rownames(result_Contrasts_matrix) <- colnames(design_input)
  
  
    comb <-combinations(n,2,v=colnames(design_input) ) 
    colnames( result_Contrasts_matrix ) <- paste0(comb[,1]," - ",comb[,2])
  
    for(k in 1:ncols_contrast_mat){
      jo1<- grep( comb[k,1] , rownames(result_Contrasts_matrix) )
      result_Contrasts_matrix[jo1,k] <- 1
      jo_1 <- grep( comb[k,2] , rownames( result_Contrasts_matrix ) )
      result_Contrasts_matrix[jo_1,k] <- -1

   #if(length(which(colnames(design_input) %in% "Control"))!=0 ){  
   # print("There is not a label named Control") if / else 
  }
  names(attributes( result_Contrasts_matrix )$dimnames) <- c( "Levels" , "Contrasts" )
  return(result_Contrasts_matrix)
  } 
