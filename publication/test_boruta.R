setwd(getSrcDirectory(function(){})[1])
source("./source_files.R")

#### Load ####

library(simpleRFNetwork)
library(SeqNet)
library(parallel)
library(pracma)
library(tictoc)
library(Rfast)
library(matrixcalc)

# Generate Network Data

dat <- genGeneNetworkData(
  n_networks = 1,
  n_genes = 1000,
  n_samples = 500,
  n_disease_modules = 2,
  main_disease_gene = F,
  prop_disease_genes = 0.5,
  average_beta = 2,
  num_threads = 4
)

# set scenarios
n_networks = c(100)
n_genes = c(1000)
n_samples = c(1000)
n_disease_modules = c(1, 2)
main_disease_gene = c(F)
# main_disease_gene = c(T, F)
# average_beta = c(1)
average_beta = c(0.5, 1, 2)

method = "LDA"
importance = "Gini"
n_iterations = 20

saveroot = "./file.Rdata"

print(paste0("This is test run"))
borutares = list()
save(borutares, file = saveroot)

# run boruta
borutares = mclapply(
  1:length(dat),
  function(i) {
    tic()
    print(paste("Boruta for Network Nr", i, "started."))
    res = boruta(dat[[i]], i, method, importance, 50, 4, n_iterations, i, saveroot)
    print(paste("Boruta for Network Nr", i, "finished!"))
    toc()
    return(res)
  },
  mc.cores = 1
)
