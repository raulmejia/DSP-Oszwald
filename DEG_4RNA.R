###################################
#### This script is the differential expression 4 RNA program in https://github.com/raulmejia/DSP-Oszwald
#### 
#### Author: Raúl Mejía
#### It should receive an expression matrix, an annotation file and retrieve
####  the differential expressed genes and their metrics in a table  
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
if (!require("cowplot")) {
  install.packages("cowplot", ask =FALSE)
  library("cowplot")
}
if (!require("affy")) {
  BiocManager::install("affy", ask =FALSE)
  library("affy")
}


## Data giveb by the user
## Provide the vector of desired contrasts

# RNA degradation with "AffyRNAdeg"




mydesign <- model.matrix( ~ 0 + annot$Morphological_Categories)
colnames(mydesign) <- gsub("annot\\$Morphological_Categories","",colnames(mydesign))
contrast.matrix <- makeContrasts(paste0("Normal","-",setdiff(colnames(mydesign), "Normal")[1]),
                                 paste0("Normal","-",setdiff(colnames(mydesign), "Normal")[2]),
                                 paste0("Normal","-",setdiff(colnames(mydesign), "Normal")[3]),
                                 paste0("Normal","-",setdiff(colnames(mydesign), "Normal")[4]),
                                 paste0("Normal","-",setdiff(colnames(mydesign), "Normal")[5]),levels=mydesign)
fit <- lmFit( combat_qnormBsep , mydesign )
fita1c <-contrasts.fit(fit , contrast.matrix)
fita1c<-eBayes(fita1c)



fit<-lmFit(a1c[,1:200],design)
contrast.matrix<-makeContrasts(Normal-Medium, Normal-High, Medium-High, levels=design)
fita1c<-contrasts.fit(fit,contrast.matrix)
fita1c<-eBayes(fit2)
fita1c<-eBayes(fita1c)



table( annot$Morphological_Categories )


dim(combat_cyclic_loess_AO)


combat_qnormBsep

# Create design matrix for leukemia study
testchr<- as.character(annot$Morphological_Categories)
testchr[ testchr == "Normal"  ] <- "Normal"
testchr[ testchr != "Normal"  ] <- "disease"
test2df   <- data.frame(testchr)
colnames(test2df) <- "Mycat"

design <- model.matrix(~Mycat, data = test2df)

# Count the number of samples modeled by each coefficient
colSums(design)

# Load package
library(limma)

# Fit the model
fit <- lmFit(combat_qnormBsep, design)

# Calculate the t-statistics
fit <- eBayes(fit)

# Summarize results
results <- decideTests(fit[, "Mycatnormal"])
results <- decideTests(fit[, "MycatNormal"])

#options(digits=2) 

genes<- topTable(fit, coef=2, n=60, adjust="BH") 
summary(results)


#---------------------------------
Morph_cat_chr <- as.character(annot$Morphological_Categories)
log_Norm_necr_w_cell <- (grepl("Normal",Morph_cat_chr ) | grepl("necrosis_w_cell_cresc",Morph_cat_chr ) )
Normal_nec_w_cresc_df <- data.frame(Morph_cat_chr [log_Norm_necr_w_cell] )
colnames(Normal_nec_w_cresc_df) <- "Mycat"

design_crec_w <- model.matrix(~Mycat, data = Normal_nec_w_cresc_df)

# Fit the model
fit_crec_w <- lmFit( combat_qnormBsep[,log_Norm_necr_w_cell], design_crec_w )

# Calculate the t-statistics
fit_crec_w <- eBayes(fit_crec_w)

# Summarize results
results <- decideTests(fit_crec_w[, "MycatNormal"])

#options(digits=2) 

genes <- topTable(fit, coef=2, n=60, adjust="BH") 
genes[,c("logFC","adj.P.Val")]
summary(results)


#--------------------------------
Morph_cat_chr <- as.character(annot$Morphological_Categories)
Normal_nec_w_cresc[ Normal_nec_w_cresc == "Normal"  ] <- "Normal"
log_fibrocell_Crescent <- (grepl("Normal",Morph_cat_chr ) | grepl("fibrocell_Crescent",Morph_cat_chr ) )
Normal_fibrocell_Crescent_df <- data.frame(Morph_cat_chr[log_fibrocell_Crescent] )
colnames(Normal_fibrocell_Crescent_df) <- "Mycat"

design_ibrocell_Crescent <- model.matrix(~Mycat, data = Normal_fibrocell_Crescent_df)

# Fit the model
fit_fibrocell_Crescent <- lmFit( combat_qnormBsep[,log_fibrocell_Crescent] , design_ibrocell_Crescent )

# Calculate the t-statistics
fit_fibrocell_Crescent <- eBayes(fit_fibrocell_Crescent)

# Summarize results
results_fibrocell_Crescent <- decideTests(fit_fibrocell_Crescent[, "MycatNormal"])

#options(digits=2) 

genes_fibrocell_Crescent <- topTable(fit_fibrocell_Crescent , coef=2, n=60, adjust="BH") 
genes_fibrocell_Crescent[,c("logFC","adj.P.Val")]
summary(results)


#-----------------------------

Morph_cat_chr <- as.character(annot$Morphological_Categories)
log_necrosis_only <- (grepl("Normal",Morph_cat_chr ) | grepl("necrosis_only",Morph_cat_chr ) )
Normal_necrosis_only_df <- data.frame(Morph_cat_chr[log_necrosis_only] )
colnames(Normal_necrosis_only_df) <- "Mycat"

design_necrosis_only <- model.matrix(~Mycat, data = Normal_necrosis_only_df)

# Fit the model
fit_necrosis_only <- lmFit( combat_qnormBsep[,log_necrosis_only] , design_necrosis_only )

# Calculate the t-statistics
fit_necrosis_only <- eBayes(fit_necrosis_only)

# Summarize results
results_necrosis_only <- decideTests(fit_necrosis_only[, "MycatNormal"])

#options(digits=2) 

genes_necrosis_only <- topTable(fit_necrosis_only , coef=2, n=60, adjust="BH") 

genes_necrosis_only[,c("logFC","adj.P.Val")]
summary( results_necrosis_only )

