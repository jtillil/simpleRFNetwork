library(SeqNet)
library(parallel)

setwd(getSrcDirectory(function(){})[1])
source("../R/genGeneNetworkDataCLassification.R")

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
scenarios = rbind(scenarios, c(100, 3000, 1000, 0, F, 0))

scenarios = scenarios[1,]
seed = 1

# generate networks
for (i in 1:nrow(scenarios)) {
  # read scenario
  scenario = scenarios[i,]
  
  # generate networkdat
  dat = genGeneNetworkDataClassification(
    n_networks = scenario$n_networks,
    n_genes = scenario$n_genes,
    n_samples = scenario$n_samples,
    n_disease_modules = scenario$n_disease_modules,
    main_disease_gene = scenario$main_disease_gene,
    average_beta = scenario$average_beta
  )
  
  # build saveroot
  saveroot = paste0(
    "./data/ndclassif",
    "_nn", scenario$n_networks,
    "_ng", scenario$n_genes,
    "_ns", scenario$n_samples,
    "_ndm", scenario$n_disease_modules,
    "_mdg", scenario$main_disease_gene,
    "_ab", scenario$average_beta,
    ".Rdata"
  )
  
  # save networkdat
  save(dat, file = saveroot)
}
