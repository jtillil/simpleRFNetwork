library(parallel)
library(pracma)
library(tictoc)
library(Rfast)
library(matrixcalc)
library(simpleRFNetwork)

setwd(getSrcDirectory(function(){})[1])

# set scenarios
n_networks = c(2)
n_genes = c(1000)
n_samples = c(1000)
n_disease_modules = c(1, 2)
main_disease_gene = c(F)
# main_disease_gene = c(T, F)
# average_beta = c(1)
average_beta = c(0.5, 1, 2)

scenarios = expand.grid(
  n_networks = n_networks,
  n_genes = n_genes,
  n_samples = n_samples,
  n_disease_modules = n_disease_modules,
  main_disease_gene = main_disease_gene,
  average_beta = average_beta
)
scenarios = rbind(scenarios, c(100, 1000, 1000, 0, F, 0))
# scenarios = rbind(scenarios, c(100, 3000, 1000, 0, F, 0))

scenarios = scenarios[1,]

# generate networks
for (i in 1:nrow(scenarios)) {
  # read scenario
  scenario = scenarios[i,]
  
  # load data
  datroot = paste0(
    "./data/ndclassif",
    "_nn", n_networks,
    "_ng", n_genes,
    "_ns", n_samples,
    "_ndm", n_disease_modules,
    "_mdg", main_disease_gene,
    "_ab", average_beta,
    ".Rdata"
  )
  load(datroot)
  
  # run boruta
  borutares = boruta(dat[[1]], "LDA", 50, 1)
  
  # save results
  setwd(getSrcDirectory(function(){})[1])
  save(dat, file = saveroot)
}


