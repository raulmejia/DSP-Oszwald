# result_dir <- paste0( path_Results_directory,"/Exploratory" )
# exp_matrix <- Raw_expmat
# annotdf <- annot_4_plotting_pca 
# melteddf <- meltedrawdata 
# label4title <- paste( data_label, "Data as given" )



PCA_box_density_plots <- function(result_dir, exp_matrix, annotdf, melteddf, label4title  ){
  exp_matrix_T <- t(exp_matrix)
  dir.create( result_dir, recursive = TRUE )
  pdf( file=paste0(result_dir,"/" , label4title ,".pdf"),
       width = 10, height = 7)
  print(autoplot( prcomp( exp_matrix_T ), data = annotdf, colour= 'Scan_ID') +
          ggtitle(paste(label4title )))
  print(autoplot( prcomp( exp_matrix_T ), data = annotdf , colour= 'Histology_number', label = TRUE, label.size = 3) +
          ggtitle(paste( label4title  )))
  print(autoplot( prcomp( exp_matrix_T ), data = annotdf , colour= 'Biopsy_year', label = TRUE, label.size = 3) +
          ggtitle(paste( label4title  ) ) )
  print(  autoplot( prcomp( exp_matrix_T ), data = annotdf , colour= 'Morphological_Categories', label = TRUE, label.size = 3) +
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
  #pheatmap( exp_matrix ,cutree_cols = 5,  annotation_col  = annotdf[,c("Histology_number","Morphological_Categories","Scan_ID" )], fontsize = 4) 
  pheatmap( exp_matrix ,cutree_cols = 5, col = brewer.pal(  length(table(annotdf[,"Morphological_Categories"])) ,"Set3"),  annotation_col  = annotdf[,c("Histology_number","Morphological_Categories","Scan_ID" )], fontsize = 4) 

  tsne_model_1 = Rtsne(  t(exp_matrix)  , check_duplicates=FALSE, pca=TRUE, perplexity=5, theta=0.5, dims=2)
  d_tsne_1 = as.data.frame( tsne_model_1$Y )
  
  d_tsne_1_simplecols <- cbind(d_tsne_1, annotdf$Scan_ID)
  mytsneplot_Nocolors <- ggplot(d_tsne_1_simplecols, aes(x=V1, y=V2 )) +
    geom_point(size=2.5) +
    guides(colour=guide_legend(override.aes=list(size=6))) +
    xlab("") + ylab("") +
    ggtitle( paste("t-SNE",label4title) ) +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) +
    scale_colour_brewer(palette = "Set3")
  print(mytsneplot_Nocolors)
  
  d_tsne_1_simplecols <- cbind(d_tsne_1, annotdf$Scan_ID)
  colnames(d_tsne_1_simplecols)[3] <- "Scan_ID"
  mytsneplot_colScanID <- ggplot(d_tsne_1_simplecols, aes(x=V1, y=V2 , col=Scan_ID )) +
    geom_point(size=2.5) +
    guides(colour=guide_legend(override.aes=list(size=6))) +
    xlab("") + ylab("") +
    ggtitle( paste("t-SNE",label4title) ) +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) +
    scale_colour_brewer(palette = "Set3")
  print(mytsneplot_colScanID)
  
  d_tsne_1_simplecols <- cbind(d_tsne_1, annotdf$Morphological_Categories  )
  colnames(d_tsne_1_simplecols)[3] <- "Morphological_Categories"
  mytsneplot_colMorphological_Categories   <- ggplot(d_tsne_1_simplecols, aes(x=V1, y=V2 , col=Morphological_Categories   )) +
    geom_point(size=2.5) +
    guides(colour=guide_legend(override.aes=list(size=6))) +
    xlab("") + ylab("") +
    ggtitle( paste("t-SNE",label4title) ) +
    theme_light(base_size=20) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank()) +
    scale_colour_brewer(palette = "Set3")
  print( mytsneplot_colMorphological_Categories )
  
  
  
  
  
  
  dev.off()  
}
