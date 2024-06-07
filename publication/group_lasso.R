setwd(getSrcDirectory(function(){})[1])
source("./source_files.R")

library(MLGL)

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

# go through scenarios
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
    "./results/grplasso",
    "_nn", scenario$n_networks,
    "_ng", scenario$n_genes,
    "_ns", scenario$n_samples,
    "_ndm", scenario$n_disease_modules,
    "_mdg", scenario$main_disease_gene,
    "_pdg", scenario$prop_disease_genes,
    "_ab", scenario$average_beta,
    ".Rdata"
  )
  grplassores = list()
  # save(grplassores, file = saveroot)
  
  # run grp lasso
  grplassores = lapply(
    1:length(dat),
    function(i) {
      tic()
      print(paste("GRP lasso for Network Nr", i, "started."))
      
      netdat = dat[[i]]$data[, -1]
      
      var_gglasso = c()
      group_gglasso = c()
      for (nmodule in 1:length(dat[[i]]$modules)) {
        for (variable in dat[[i]]$modules[[nmodule]]) {
          var_gglasso = c(var_gglasso, variable)
          group_gglasso = c(group_gglasso, nmodule)
        }
      }
      
      gr = overlapgglasso(X = as.matrix(dat[[i]]$data[, -1]),
            y = 2*as.numeric(dat[[i]]$data[, 1]) - 3,
            var = var_gglasso,
            group = group_gglasso,
            # lambda = 1,
            loss="logit",
            intercept = F)
      
      print(paste("GRP lasso for Network Nr", i, "finished!"))
      toc()
      
      return(list(
        selected = unique(gr$group$s20),
        causal = dat[[i]]$causal_modules,
        number_correctly_selected = sum(dat[[i]]$causal_modules %in% unique(gr$group$s50))
      ))
    }
  )
  
  # save grp lasso
  save(grplassores, file = saveroot)
  
  # clear dat
  rm(dat)
}
