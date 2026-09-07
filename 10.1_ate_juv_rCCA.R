setwd("/home/jc828813/RNAseq/results/aten_juveniles/mixomics")

library(limma)
library(Glimma)
library(edgeR)
library(tidyverse)
library(mixOmics)
set.seed(5249)

#load sample info
coldata<- read.csv("aten_juvenile_read_info2.csv", row.names="ID") |> filter(Symbiont == "D1")
coldata <- coldata[,c("Cross","Temp","Symbiont","group")]
#convert columns to factors
coldata$Cross <- factor(coldata$Cross)
coldata$group <- factor(coldata$group)
coldata$Temp <- factor(coldata$Temp)
coldata$Symbiont <- factor(coldata$Symbiont)

#make sample lists for each symbiont treatment
coldata2<- coldata |> rownames_to_column(var = "Samples")
D1samples<- coldata2 |> dplyr::filter(Symbiont == "D1")
D1samples<- D1samples$Samples

#load data
x1<- read.csv("x_D1_host.csv", row.names = 1) |> data.matrix()
x2<- read.csv("x_D1_symb.csv", row.names = 1) |> data.matrix()

#define groups
group <- as.factor(c("AmbientD1", "AmbientD1", "AmbientD1", "AmbientD1", "AmbientD1", "HotD1", "AmbientD1", "AmbientD1", "HotD1", "AmbientD1", "AmbientD1", "HotD1", "HotD1", "AmbientD1", "HotD1", "HotD1", "AmbientD1", "HotD1", "AmbientD1", "HotD1", "HotD1", "AmbientD1", "HotD1", "HotD1", "HotD1", "AmbientD1", "HotD1", "HotD1", "HotD1", "AmbientD1", "AmbientD1", "AmbientD1", "HotD1", "AmbientD1", "HotD1", "HotD1"))
y<-group

#Initial Analysis
#Exploration of feature correlation
# produce a heat map of the cross correlation matrix
#pdf("imgCor.pdf", width = 10, height = 10)
#imgCor(x1, x2, sideColors = c("purple", "green"))
#dev.off()

# set grid search values for each regularisation parameter
grid1 <- seq(0.001, 0.2, length = 5) 
grid2 <- seq(0.001, 0.2, length = 5)

# optimise the regularisation parameter values
pdf("cv_tune_rcc_nutrimouse.pdf", width = 10, height = 10)
cv.tune.rcc.nutrimouse <- tune.rcc(x1, x2, grid1 = grid1, grid2 = grid2, 
                                   validation = "loo")
dev.off()

cv.tune.rcc.nutrimouse # examine the results of CV tuning

opt.l1 <- cv.tune.rcc.nutrimouse$opt.lambda1 # extract the optimal lambda values
opt.l2 <- cv.tune.rcc.nutrimouse$opt.lambda2

# formed optimised CV rCCA
CV.rcc.nutrimouse <- rcc(x1, x2, method = "ridge", 
                         lambda1 = opt.l1, lambda2 = opt.l2) 

# run the rCCA method using shrinkage
shrink.rcc.nutrimouse <- rcc(x1,x2, method = 'shrinkage') 
# examine the optimal lambda values after shrinkage 
shrink.rcc.nutrimouse$lambda 

# barplot of cross validation method rCCA canonical correlations
pdf("CV_rcc_nutrimouse_bar.pdf", width = 10, height = 10)
plot(CV.rcc.nutrimouse, type = "barplot", main = "Cross Validation") 
dev.off()

# barplot of shrinkage method rCCA canonical correlations
pdf("shrink_rcc_nutrimouse_bar.pdf", width = 10, height = 10)
plot(shrink.rcc.nutrimouse, type = "barplot", main = "Shrinkage") 
dev.off()

#Final Model
#Sample plots
# plot the projection of samples for CV rCCA data
pdf("CV_rcc_nutrimouse_plot.pdf", width = 10, height = 10)
plotIndiv(CV.rcc.nutrimouse, comp = 1:2, 
          ind.names = nutrimouse$genotype,
          group = nutrimouse$diet, rep.space = "XY-variate", 
          legend = TRUE, title = '(a) rCCA CV XY-space')
dev.off()

# plot the projection of samples for shrinkage rCCA data
pdf("shrink_rcc_nutrimouse_plot.pdf", width = 10, height = 10)
plotIndiv(shrink.rcc.nutrimouse, comp = 1:2, 
          ind.names = nutrimouse$genotype,
          group = nutrimouse$diet, rep.space = "XY-variate", 
          legend = TRUE, title = '(b) rCCA shrinkage XY-space')
dev.off()

# plot the arrow plot of samples for CV rCCA data
pdf("CV_rcc_nutrimouse_arrow.pdf", width = 10, height = 10)
plotArrow(CV.rcc.nutrimouse, group = nutrimouse$diet, 
          col.per.group = color.mixo(1:5),
          title = '(a) CV method')
dev.off()

# plot the arrow plot of samples for shrinkage rCCA data
pdf("shrink_rcc_nutrimouse_arrow.pdf", width = 10, height = 10)
plotArrow(shrink.rcc.nutrimouse, group = nutrimouse$diet, 
          col.per.group = color.mixo(1:5),
          title = '(b) shrinkage method')
dev.off()

#circle plots
pdf("CV_rcc_nutrimouse_circle.pdf", width = 10, height = 10)
plotVar(CV.rcc.nutrimouse, var.names = c(TRUE, TRUE),
        cex = c(4, 4), cutoff = 0.5,
        title = '(a) rCCA CV comp 1 - 2')
dev.off()

pdf("shrink_rcc_nutrimouse_circle.pdf", width = 10, height = 10)
plotVar(shrink.rcc.nutrimouse, var.names = c(TRUE, TRUE),
        cex = c(4, 4), cutoff = 0.5,
        title = '(b) rCCA shrinkage comp 1 - 2')
dev.off()

pdf("CV_rcc_nutrimouse_network.pdf", width = 10, height = 10)
network(CV.rcc.nutrimouse, comp = 1:2, interactive = FALSE,
        lwd.edge = 2,
        cutoff = 0.5)
dev.off()

pdf("CV_rcc_nutrimouse_heat.pdf", width = 10, height = 10)
cim(CV.rcc.nutrimouse, comp = 1:2, xlab = "genes", ylab = "lipids")
dev.off()
