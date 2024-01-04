setwd(getSrcDirectory(function(){})[1])
source("./source_files.R")

# set scenarios
n_networks = c(100)
n_genes = c(1000)
n_samples = c(1000)
n_disease_modules = c(1, 2)
main_disease_gene = c(F)
# main_disease_gene = c(T, F)
# average_beta = c(1)
average_beta = c(0.5, 1, 2)
prop_disease_genes = c(0.5)

scenarios = expand.grid(
  n_networks = n_networks,
  n_genes = n_genes,
  n_samples = n_samples,
  n_disease_modules = n_disease_modules,
  main_disease_gene = main_disease_gene,
  average_beta = average_beta,
  prop_disease_genes = prop_disease_genes
)
scenarios = rbind(scenarios, c(100, 1000, 1000, 0, F, 0, 0.5))
# scenarios = rbind(scenarios, c(100, 3000, 1000, 0, F, 0, 1))

method = "LDA"
importance = "permutation"
n_iterations = 20

print(paste0("This is run ", method, " Classification"))
# generate networks
for (i in 1:nrow(scenarios)) {
  print(i)
  
  # read scenario
  scenario = scenarios[i,]
  
  # load data
  datroot = paste0(
    "./data/ndclassif",
    # "_datasc", data_scenario,
    "_nn", scenario$n_networks,
    "_ng", scenario$n_genes,
    "_ns", scenario$n_samples,
    "_ndm", scenario$n_disease_modules,
    "_mdg", scenario$main_disease_gene,
    "_pdg", scenario$prop_disease_genes,
    "_ab", scenario$average_beta,
    ".Rdata"
  )
  print(getwd())
  print(datroot)
  load(datroot)
  
  # save root
  saveroot = paste0(
    "./results/resclassif",
    # "_datasc", data_scenario,
    "_", method,
    "_", importance,
    "_ni", n_iterations,
    "_nn", scenario$n_networks,
    "_ng", scenario$n_genes,
    "_ns", scenario$n_samples,
    "_ndm", scenario$n_disease_modules,
    "_mdg", scenario$main_disease_gene,
    "_pdg", scenario$prop_disease_genes,
    "_ab", scenario$average_beta,
    ".Rdata"
  )
  borutares = list()
  save(borutares, file = saveroot)
  
  # run boruta
  borutares = mclapply(
    1:length(dat),
    function(i) {
      tic()
      print(paste("Boruta for Network Nr", i, "started."))
      res = boruta(dat[[i]], i, method, importance, 500, 1, n_iterations, i, saveroot)
      print(paste("Boruta for Network Nr", i, "finished!"))
      toc()
      return(res)
    },
    mc.cores = 40
  )
  
  # clear dat
  rm(dat)
}
