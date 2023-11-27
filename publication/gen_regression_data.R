library(SeqNet)
library(parallel)
library(pracma)
library(tictoc)
library(simpleRFNetwork)

setwd(getSrcDirectory(function(){})[1])

# set scenarios
n_networks = c(100)
n_genes = c(1000)
n_samples = c(1000)
disease_modules = c(T, F)

scenarios = expand.grid(
  n_networks = n_networks,
  n_genes = n_genes,
  n_samples = n_samples,
  disease_modules = disease_modules
)

# scenarios = scenarios[1,]

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
    num_threads = 60
  )
  
  # build saveroot
  saveroot = paste0(
    "./data/ndregress",
    "_nn", scenario$n_networks,
    "_ng", scenario$n_genes,
    "_ns", scenario$n_samples,
    "_dm", scenario$disease_modules,
    ".Rdata"
  )
  
  # save networkdat
  save(dat, file = saveroot)
}
toc()

