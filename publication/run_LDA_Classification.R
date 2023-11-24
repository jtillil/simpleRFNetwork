library(parallel)
library(pracma)
library(tictoc)
library(Rfast)
library(matrixcalc)
library(simpleRFNetwork)

# set scenarios
n_networks = c(2)
n_genes = c(1000, 3000)
n_samples = c(1000)
n_disease_modules = c(1, 2)
main_disease_gene = c(F)
# main_disease_gene = c(T, F)
average_beta = c(1)
# average_beta = c(0.5, 1, 2)


# load data
datroot = paste0(
  "./data/ndclassif",
  "_nn", scenario$n_networks,
  "_ng", scenario$n_genes,
  "_ns", scenario$n_samples,
  "_ndm", scenario$n_disease_modules,
  "_mdg", scenario$main_disease_gene,
  "_ab", scenario$average_beta,
  ".Rdata"
)
setwd(getSrcDirectory(function(){})[1])
load(datroot)

rf = boruta(dat[[1]], "LDA", 50, 1)

