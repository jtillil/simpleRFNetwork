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

fraction_of_positives = c()
for (i in 1:nrow(scenarios)) {
  # read scenario
  scenario = scenarios[i,]
  print(scenario)
  
  datroot = paste0(
    # "./resclassif",
    "./data/ndclassif",
    "_nn", 100,
    "_ng", 1000,
    "_ns", 1000,
    "_ndm", scenario$n_disease_modules,
    "_mdg", 0,
    "_pdg", 0.5,
    "_ab", scenario$average_beta,
    ".Rdata"
  )
  load(datroot)
  
  fractions_in_scenario = c()
  for (j in 1:100) {
    fractions_in_scenario = c(fractions_in_scenario, sum(dat[[j]]$data[,1] == 1) / 1000)
  }
  
  fraction_of_positives = c(fraction_of_positives, mean(fractions_in_scenario))
}
