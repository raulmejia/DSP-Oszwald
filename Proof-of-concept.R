rm(list = ls())
library("limma")

# the name of your groups could be the non duplicated elements of
# annotations dataframe

sd <- 0.3*sqrt(4/rchisq(100,df=4))
y <- matrix(rnorm(100*6,sd=sd),100,6)
rownames(y) <- paste("Gene",1:100)
y[1:2,4:6] <- y[1:2,4:6] + 2

# An appropriate desing matrix can be created and linear model
# fitted using:
design <- model.matrix(~ 0+factor(c(1,1,2,2,3,3)))
colnames(design) <- c("group1", "group2", "group3")
fit <- lmFit(y, design)
efit <- eBayes(fit)

# To do all the possible pair-wise comparisons the appropriate
# contrast matrix can be created by:
contrast.matrix <- makeContrasts(group2-group1, group3-group2, group3-group1, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef=1, adjust="BH", number = 10)
topTable(fit2, coef=2, adjust="BH")
topTable(fit2, coef=3, adjust="BH")

results <- decideTests(fit2)
?decideTests

vennDiagram(results)
