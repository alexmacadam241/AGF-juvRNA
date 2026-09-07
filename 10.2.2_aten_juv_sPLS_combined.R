setwd("/Users/alexmacadam/Dropbox/PhD/Chapter4_RNA/Scripts/Codes for RNAseq analysis")

library(limma)
library(Glimma)
library(edgeR)
library(tidyverse)
library(mixOmics)
set.seed(5249)

#load sample info
coldata<- read.csv("/Users/alexmacadam/Dropbox/PhD/Chapter4_RNA/Data/Gene_counts/read_info/aten_juvenile_read_info2.csv", row.names="ID") |> filter(Symbiont == "D1")
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
x1<- read.csv("//Users/alexmacadam/Dropbox/PhD/Chapter4_RNA/Data/Gene_counts/aten_juveniles/mixomics_combined/x_D1_host.csv", row.names = 1) |> data.matrix()
x2<- read.csv("/Users/alexmacadam/Dropbox/PhD/Chapter4_RNA/Data/Gene_counts/aten_juveniles/mixomics_combined/x_D1_symb.csv", row.names = 1) |> data.matrix()

#define groups
group <- as.factor(c("AmbientD1", "AmbientD1", "AmbientD1", "AmbientD1", "AmbientD1", "HotD1", "AmbientD1", "AmbientD1", "HotD1", "AmbientD1", "AmbientD1", "HotD1", "HotD1", "AmbientD1", "HotD1", "HotD1", "AmbientD1", "HotD1", "AmbientD1", "HotD1", "HotD1", "AmbientD1", "HotD1", "HotD1", "HotD1", "AmbientD1", "HotD1", "HotD1", "HotD1", "AmbientD1", "AmbientD1", "AmbientD1", "HotD1", "AmbientD1", "HotD1", "HotD1"))
y<-group

###############
#choose componenets X
dfX <- data.frame(comp1 = 34, comp2 = 35)

# Assign it to tune.spls.liver$choice.keepX
tune.spls.liverX <- list()  # Initialize an empty list if it's not already created
tune.spls.liverX$choice.keepX <- dfX

# Check the result
tune.spls.liverX$choice.keepX

#choose componenets Y
dfY <- data.frame(comp1 = 3, comp2 = 6)

# Assign it to tune.spls.liver$choice.keepX
tune.spls.liverY <- list()  # Initialize an empty list if it's not already created
tune.spls.liverY$choice.keepY <- dfY

# Assign it to tune.spls.liver$choice.keepY
tune.spls.liverY <- list()  # Initialize an empty list if it's not already created
tune.spls.liverY$choice.keepY <- dfY
# Check the result
tune.spls.liverY$choice.keepY

# extract optimal number of variables for X dataframe
optimal.keepX <- tune.spls.liverX$choice.keepX

# extract optimal number of variables for Y datafram
optimal.keepY <- tune.spls.liverY$choice.keepY

optimal.ncomp <-  length(optimal.keepX) # extract optimal number of components

#Final Model
# use all tuned values from above
#final.spls.liver <- spls(x1, x2, ncomp = optimal, 
#                         keepX = optimal.keepX,
#                         keepY = optimal.keepY,
#                         mode = "regression") # explanitory approach being used, 
# hence use regression mode

final.spls.liver <- spls(x1, x2, ncomp = 2, 
                         keepX = c(34,35),
                         keepY = c(3,6),
                         mode = "regression") # explanitory approach being used, 
# hence use regression mode
#saveRDS(final.spls.liver, file = "/Users/alexmacadam/Dropbox/PhD/Chapter4_RNA/Scripts/Codes for RNAseq analysis/10.2.2_aten_host_mixomics_D1_final_spls.rds")
final.spls.liver<- readRDS("/Users/alexmacadam/Dropbox/PhD/Chapter4_RNA/Scripts/Codes for RNAseq analysis/10.2.2_aten_host_mixomics_D1_final_spls.rds")

#Plots
pdf("final_spls_liver_X.pdf", width = 10, height = 10)
plotIndiv(final.spls.liver, ind.names = FALSE, 
          rep.space = "X-variate", # plot in X-variate subspace
          group = coldata$group, # colour by time group
          pch = as.factor(coldata$group), 
          legend = TRUE, legend.title = 'Temp', legend.title.pch = 'Temp')
dev.off()

pdf("final_spls_liver_Y.pdf", width = 10, height = 10)
plotIndiv(final.spls.liver, ind.names = FALSE,
          rep.space = "Y-variate", # plot in Y-variate subspace
          group = coldata$group, # colour by time group
          pch = as.factor(coldata$group), 
          #col.per.group = color.mixo(1:4), 
          legend = TRUE, legend.title = 'Temp', legend.title.pch = 'Temp')
dev.off()

pdf("final_spls_liver_XY.pdf", width = 10, height = 10)
plotIndiv(final.spls.liver, ind.names = FALSE, 
          rep.space = "XY-variate", # plot in averaged subspace
          group = coldata$group, # colour by time group
          pch = as.factor(coldata$group), # select symbol
          #col.per.group = color.mixo(1:4),                      # by dose group
          legend = TRUE, legend.title = 'Temp', legend.title.pch = 'Temp')
dev.off()

pdf("final_spls_liver_arrow.pdf", width = 10, height = 10)
plotArrow(final.spls.liver, ind.names = FALSE,
          group = coldata$group, # colour by time group
          #col.per.group = color.mixo(1:4),
          legend.title = 'Temp')
dev.off()

# form new perf() object which utilises the final model
perf.spls.liver <- perf(final.spls.liver, 
                        folds = 5, nrepeat = 50, # use repeated cross-validation
                        validation = "Mfold", 
                        dist = "max.dist",  # use max.dist measure
                        progressBar = TRUE)

#saveRDS(perf.spls.liver, file = "/Users/alexmacadam/Dropbox/PhD/Chapter4_RNA/Scripts/Codes for RNAseq analysis/10.2.2_aten_host_mixomics_D1_perf_splsda.rds")
perf.spls.liver<- readRDS("/Users/alexmacadam/Dropbox/PhD/Chapter4_RNA/Scripts/Codes for RNAseq analysis/10.2.2_aten_host_mixomics_D1_perf_splsda.rds")
# plot the stability of each feature for the first two components, 
# 'h' type refers to histogram
pdf("final_spls_liver_histogram.pdf", width = 10, height = 10)
par(mfrow=c(1,2)) 
plot(perf.spls.liver$features$stability.X[[1]], type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(a) Comp 1', las =2,
     xlim = c(0, 150))
plot(perf.spls.liver$features$stability.X$comp2, type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = '(b) Comp 2', las =2,
     xlim = c(0, 300))
dev.off()

pdf("final_spls_liver_circle.pdf", width = 10, height = 10)
plotVar(final.spls.liver, cex = c(3,4), var.names = c(FALSE, TRUE))
dev.off()


color.edge <- color.GreenRed(50)  # set the colours of the connecting lines

# X11() # To open a new window for Rstudio
network(final.spls.liver, comp = 1:2,
        cutoff = 0.7, # only show connections with a correlation above 0.7
        shape.node = c("rectangle", "circle"),
        color.node = c("cyan", "pink"),
        color.edge = color.edge,
        save = 'png', # save as a png to the current working directory
        name.save = 'sPLS Network Plot')

