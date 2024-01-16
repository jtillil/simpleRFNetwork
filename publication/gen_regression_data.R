library(SeqNet)
library(parallel)
library(pracma)
library(tictoc)
library(simpleRFNetwork)

setwd(getSrcDirectory(function(){})[1])

# set scenarios
n_networks = c(2)
n_genes = c(1000)
n_samples = c(1000)
disease_modules = c(T, F)
prop_disease_genes = c(0.5)

scenarios = expand.grid(
  n_networks = n_networks,
  n_genes = n_genes,
  n_samples = n_samples,
  disease_modules = disease_modules,
  prop_disease_genes = prop_disease_genes
)

scenarios = scenarios[1,]

tic()
# generate networks
for (i in 1:nrow(scenarios)) {
  print(i)
  
  # read scenario
  scenario = scenarios[i,]
  
  # generate networkdat
  dat = genGeneNetworkDataRegression(
    n_networks = scenario$n_networks,
    n_genes = scenario$n_genes,
    n_samples = scenario$n_samples,
    disease_modules = scenario$disease_modules,
    prop_disease_genes = scenario$prop_disease_genes,
    num_threads = 2
  )
  
  # build saveroot
  saveroot = paste0(
    "./data/ndregress",
    "_nn", scenario$n_networks,
    "_ng", scenario$n_genes,
    "_ns", scenario$n_samples,
    "_dm", scenario$disease_modules,
    "_pdg", scenario$prop_disease_genes,
    ".Rdata"
  )
  
  # save networkdat
  save(dat, file = saveroot)
}
toc()

