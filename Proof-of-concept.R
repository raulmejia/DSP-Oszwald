rm(list = ls())
library("limma")
library("tidyverse")

# the name of your groups could be the non duplicated elements of
# annotations dataframe

sd <- 0.3*sqrt(4/rchisq(100,df=4))
y <- matrix(rnorm(100*6,sd=sd),60,10)
rownames(y) <- paste("Gene",1:60)
colnames(y) <- paste("Sample" ,1:10 )
y[1:10,1:3] <- y[1:10,1:3] + 2

annotdf <- data.frame(c(paste("Sample" ,1:10 )), 
                      c("ein","zwei","ein","zwei","ein","zwei","drei","vier","drei","vier"))
colnames(annotdf) <- c("samples","labels")

# ---

designmat <- model.matrix( ~0+  as.factor( annotdf$labels) )  # design matrix
    colnames(designmat) <- levels(as.factor( annotdf$labels))
´`´`++**´``´´´´++´stats::model.matrix
    
myfit <- lmFit(y, designmat) # fitting the linear model

mycontMat <- mymakeContrasts(designmat) # make your matrix of contrasts

myefit <- eBayes(myfit) # moderate standard errors of the estimated log-fold changes

topTable( myefit , coef=1, adjust="BH", number = 20)
topTable( myefit , coef=2, adjust="BH", number = 10)

# write.table

#results <- decideTests(myefit)
#vennDiagram(results)






