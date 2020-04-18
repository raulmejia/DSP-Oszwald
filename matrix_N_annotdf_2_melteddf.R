
matrix_N_annotdf_2_melteddf <- function( onematrix, oneannotdf){
  table_with_uniqID <- cbind( as.character(oneannotdf$Unique_ID) , as.data.frame(onematrix) )
  
  colnames(table_with_uniqID) <- c("Unique_ID", colnames(table[ ,colpositions_withgenes]))
  table_melted_MtxAndUniqueID <-melt(data=table_with_uniqID, id.vars="Unique_ID",measure.vars = colnames(table_with_uniqID)[-1]   )
  
  colnames(table_melted_MtxAndUniqueID)[2] <- "gene"
  
  Scan_ID <-                  as.character( rep( oneannotdf$Scan_ID, length(colpositions_withgenes)) )
  Biopsy_year <-              as.character( rep( oneannotdf$Biopsy_year , length(colpositions_withgenes)) )
  Histology_number <-         as.character( rep( oneannotdf$Histology_number , length(colpositions_withgenes)) )
  Morphological_Categories <- as.character( rep( oneannotdf$Morphological_Categories , length(colpositions_withgenes)) )
  ROI_ID <-                   as.character( rep( oneannotdf$ROI_ID, length(colpositions_withgenes)) )
  
  result_table_melted_MtxAndScanID <- cbind(table_melted_MtxAndUniqueID , Scan_ID
                                            ,Biopsy_year,Histology_number,ROI_ID,Morphological_Categories )
  
  return(result_table_melted_MtxAndScanID)
}
