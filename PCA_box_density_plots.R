

PCA_box_density_plots <- function(result_dir, exp_matrix, annotdf, melteddf, label4title  ){
  dir.create( result_dir, recursive = TRUE )
  pdf( file=paste0(result_dir,"/" , label4title ,".pdf"),
       width = 10, height = 7)
  print(autoplot( prcomp( exp_matrix ), data = annotdf, colour= 'Scan_ID') +
          ggtitle(paste(label4title )))
  print(autoplot( prcomp( exp_matrix ), data = annotdf , colour= 'Histology_number', label = TRUE, label.size = 3) +
          ggtitle(paste( label4title  )))
  print(autoplot( prcomp( exp_matrix ), data = annotdf , colour= 'Biopsy_year', label = TRUE, label.size = 3) +
          ggtitle(paste( label4title  ) ) )
  print(  autoplot( prcomp( exp_matrix ), data = annotdf , colour= 'Morphological_Categories', label = TRUE, label.size = 3) +
            ggtitle(paste( label4title )) )
  q <- ggplot( melteddf , aes(Unique_ID, value, fill=Scan_ID))
  print( q + geom_boxplot( )+ ggtitle(paste( label4title )) )
  q <- ggplot( melteddf , aes(Unique_ID, value, fill=Histology_number))
  print( q + geom_boxplot( )+ ggtitle(paste( label4title )) )
  q <- ggplot( melteddf , aes(Unique_ID, value, fill=Biopsy_year))
  print( q + geom_boxplot( )+ ggtitle(paste( label4title )) )
  q <- ggplot( melteddf , aes(Unique_ID, value, fill= Morphological_Categories))
  print( q + geom_boxplot( )+ ggtitle(paste( label4title )) )
  
  plot1 <- ggplot(melteddf, aes(x=value, fill = Scan_ID, y = Unique_ID)) + 
    geom_density_ridges() + 
    ggtitle(label4title)
  print(plot1)
  plot2 <- ggplot(melteddf, aes(x=value, color = Scan_ID)) + 
    geom_density() + 
    ggtitle(label4title)
  print(plot2)
  plot3 <- ggplot(melteddf, aes(x=value, fill = Morphological_Categories  , y = Unique_ID)) + 
    geom_density_ridges() + 
    ggtitle(label4title)
  print(plot3)
  plot4 <- ggplot(melteddf, aes(x=value, color = Morphological_Categories )) + 
    geom_density() + 
    ggtitle(label4title)
  print(plot4)
  print(plot_grid(plot1, plot2, plot3, plot4, nrow = 2))
  pheatmap( t(exp_matrix) ,cutree_cols = 4,  annotation_col  = annotdf[,c("Histology_number","Morphological_Categories","Scan_ID" )] ) 
  
  dev.off()  
}
