#set wd
setwd("/home/jc828813/RNAseq/results/aten_juveniles/mixomics")
#load packages
library(limma)
library(Glimma)
library(edgeR)
library(tidyverse)
library(mixOmics)
set.seed(5249)

#load sample info
coldata<- read.csv("aten_juvenile_read_info2.csv", row.names="ID") |> filter(Symbiont == "SS")
coldata <- coldata[,c("Cross","Temp","Symbiont","group")]
#convert columns to factors
coldata$Cross <- factor(coldata$Cross)
coldata$group <- factor(coldata$group)
coldata$Temp <- factor(coldata$Temp)
coldata$Symbiont <- factor(coldata$Symbiont)
print("info loaded")

#make sample lists for each symbiont treatment
coldata2<- coldata |> rownames_to_column(var = "Samples")
SSsamples<- coldata2 |> dplyr::filter(Symbiont == "SS")
SSsamples<- SSsamples$Samples

#load count data
x1<- read.csv("x_SS_host.csv", row.names = 1) |> data.matrix()
x2<- read.csv("x_SS_symb.csv", row.names = 1) |> data.matrix()
print("counts loaded")

#define groups
group <- as.factor(c("AmbientSS","HotSS","HotSS","AmbientSS","AmbientSS","HotSS","AmbientSS","AmbientSS","AmbientSS","AmbientSS","HotSS","AmbientSS","HotSS","AmbientSS","AmbientSS","HotSS","HotSS","HotSS","AmbientSS","AmbientSS","HotSS","AmbientSS","HotSS","AmbientSS","HotSS","HotSS","AmbientSS","AmbientSS","HotSS","HotSS","AmbientSS","AmbientSS"))
y<-group

#initital spls
spls.liver <- spls(X = x1, Y = x2, ncomp = 5, mode = 'regression')

# repeated CV tuning of component count
perf.spls.liver <- perf(spls.liver, validation = 'Mfold',
                        folds = 5, nrepeat = 50,
                        progressBar = TRUE)

pdf("perf.spls.liver_SS.pdf", width = 10, height = 10)
plot(perf.spls.liver, criterion = 'Q2.total')
dev.off()

print("perf done")

# set range of test values for number of variables to use from X dataframe
list.keepX <- c(seq(20, 50, 5))
# set range of test values for number of variables to use from Y dataframe
list.keepY <- c(3:10) 

tune.spls.liver <- tune.spls(x1, x2, ncomp = 2,
                             test.keepX = list.keepX,
                             test.keepY = list.keepY,
                             nrepeat = 5, folds = 5,
                             mode = 'regression', measure = 'cor',
                             progressBar = TRUE) 


pdf("tune.spls.liver_SS.pdf", width = 10, height = 10)
plot(tune.spls.liver)         # use the correlation measure for tuning
dev.off()
print("spls done")