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

###############
#Preliminary Analysis with PCA
pca.host <- pca(x1, ncomp = 10, center = TRUE, scale = TRUE)
pca.sym <- pca(x2, ncomp = 10, center = TRUE, scale = TRUE)

pdf("pca_host.pdf", width = 10, height = 10)
plot(pca.host)
dev.off()

pdf("pca_sym.pdf", width = 10, height = 10)
plot(pca.sym)
dev.off()

pdf("pca_host_indiv.pdf", width = 10, height = 10)
plotIndiv(pca.host, comp = c(1, 2), 
          group = y, 
          ind.names = y, 
          legend = TRUE, 
          title = 'Host PCA comp 1 - 2')
dev.off()

pdf("pca_sym_indiv.pdf", width = 10, height = 10)
plotIndiv(pca.sym, comp = c(1, 2), 
          group = y, 
          ind.names = y, 
          legend = TRUE, title = 'Symbiont PCA comp 1 - 2')
dev.off()

#Initial sPLS model
spls.aten <- spls(X = x1, Y = x2, ncomp = 6, mode = 'regression')
saveRDS(spls.aten, file = "10.2_spls.aten.rds")

#Tuning sPLS
# repeated CV tuning of component count
perf.spls.aten <- perf(spls.aten, validation = 'Mfold',
                        folds = 4, nrepeat = 50,
                        )
saveRDS(perf.spls.aten, file = "10.2_perf.spls.aten.rds")

pdf("perf_spls_aten.pdf", width = 10, height = 10)
plot(perf.spls.aten, criterion = 'Q2.total')
dev.off()

#Selecting the number of variables
# set range of test values for number of variables to use from X dataframe
list.keepX <- c(seq(20, 50, 5))
# set range of test values for number of variables to use from Y dataframe
list.keepY <- c(3:10)


tune.spls.aten <- tune.spls(x1, x2, ncomp = 2,
                             test.keepX = list.keepX,
                             test.keepY = list.keepY,
                             nrepeat = 10, folds = 4, # use 4 folds
                             mode = 'regression', measure = 'cor') 
saveRDS(tune.spls.aten, file = "10.2_tune.spls.aten.rds")

pdf("tune_spls_aten.pdf", width = 10, height = 10)
plot(tune.spls.aten)         # use the correlation measure for tuning
dev.off()

tune.spls.aten$choice.keepX

tune.spls.aten$choice.keepY

# extract optimal number of variables for X dataframe
optimal.keepX <- tune.spls.aten$choice.keepX 

# extract optimal number of variables for Y datafram
optimal.keepY <- tune.spls.aten$choice.keepY

optimal.ncomp <-  length(optimal.keepX) # extract optimal number of components


#Final Model
# use all tuned values from above
final.spls.aten <- spls(x1, x2, ncomp = optimal.ncomp, 
                         keepX = optimal.keepX,
                         keepY = optimal.keepY,
                         mode = "regression") # explanitory approach being used, 
# hence use regression mode
saveRDS(final.spls.aten, file = "10.2_final.spls.aten.rds")

#Plots
pdf("final_spls_aten_x1.pdf", width = 10, height = 10)
plotIndiv(final.spls.aten, ind.names = FALSE, 
          rep.space = "X-variate", # plot in X-variate subspace
          group = aten.toxicity$treatment$Time.Group, # colour by time group
          pch = as.factor(aten.toxicity$treatment$Dose.Group), 
          col.per.group = color.mixo(1:4), 
          legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')
dev.off()

pdf("final_spls_aten_x2.pdf", width = 10, height = 10)
plotIndiv(final.spls.aten, ind.names = FALSE,
          rep.space = "Y-variate", # plot in Y-variate subspace
          group = aten.toxicity$treatment$Time.Group, # colour by time group
          pch = as.factor(aten.toxicity$treatment$Dose.Group), 
          col.per.group = color.mixo(1:4), 
          legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')
dev.off()

pdf("final_spls_aten_XY.pdf", width = 10, height = 10)
plotIndiv(final.spls.aten, ind.names = FALSE, 
          rep.space = "XY-variate", # plot in averaged subspace
          group = aten.toxicity$treatment$Time.Group, # colour by time group
          pch = as.factor(aten.toxicity$treatment$Dose.Group), # select symbol
          col.per.group = color.mixo(1:4),                      # by dose group
          legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')
dev.off()

pdf("final_spls_aten_arrow.pdf", width = 10, height = 10)
plotArrow(final.spls.aten, ind.names = FALSE,
          group = aten.toxicity$treatment$Time.Group, # colour by time group
          col.per.group = color.mixo(1:4),
          legend.title = 'Time.Group')
dev.off()

# form new perf() object which utilises the final model
perf2.spls.aten <- perf(final.spls.aten, 
                        folds = 4, nrepeat = 50, # use repeated cross-validation
                        validation = "Mfold", 
                        dist = "max.dist",  # use max.dist measure
                        progressBar = FALSE)
saveRDS(perf2.spls.aten, file = "10.2_perf2.spls.aten.rds")

# plot the stability of each feature for the first two components, 
# 'h' type refers to histogram
pdf("final_spls_aten_histogram.pdf", width = 10, height = 10)
par(mfrow=c(1,2)) 
plot(perf.spls.aten$features$stability.X[[1]], type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(a) Comp 1', las =2,
     xlim = c(0, 150))
plot(perf.spls.aten$features$stability.X$comp2, type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(b) Comp 2', las =2,
     xlim = c(0, 300))
dev.off()

pdf("final_spls_aten_circle.pdf", width = 10, height = 10)
plotVar(final.spls.aten, cex = c(3,4), var.names = c(FALSE, TRUE))
dev.off()

color.edge <- color.GreenRed(50)  # set the colours of the connecting lines

# X11() # To open a new window for Rstudio
pdf("final_spls_aten_network.pdf", width = 10, height = 10)
network(final.spls.aten, comp = 1:2,
        cutoff = 0.7, # only show connections with a correlation above 0.7
        shape.node = c("rectangle", "circle"),
        color.node = c("cyan", "pink"),
        color.edge = color.edge,
        save = 'png', # save as a png to the current working directory
        name.save = 'sPLS Network Plot')
dev.off()
