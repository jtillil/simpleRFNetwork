setwd(getSrcDirectory(function(){})[1])
source("./source_files.R")

library(simpleRFNetwork)
library(SeqNet)
library(parallel)
library(pracma)
library(tictoc)
library(Rfast)
library(matrixcalc)

library(ridge)
library(glmnet)

library(ggplot2)

library(ranger)

#### setup

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
importance = "permutation"

#### prediction error RF

pred_err = function(network) {
  # print(ID_network)
  #
  # network = dat[[ID_network]]
  rfdat = network$data[1:500,]
  preddat = as.matrix(network$data[501:1000,-1])
  predlabel = network$data[501:1000,1]
  modules = network$modules
  
  rf = simpleRFNetwork(
    pheno ~ .,
    data = rfdat,
    num_trees=500,
    num_threads=1,
    splitobject="module",
    splitmethod=splitmethod,
    varselection="none",
    mtry="root",
    varclusters = modules,
    seed = 1L
  )
  
  print(length(rf$trees[[1]]$split_clusterIDs))
  
  res = list()
  res$method = splitmethod
  # err = rf$predictionErrorForestAndTrees()
  # err = rf$predictionErrorForest()
  pred = rf$predict(preddat)
  res$prederr = sum(pred != predlabel) / 500
  
  split_counts = numeric(length(modules))
  for (treeID in 1:length(rf$trees)) {
    clusterIDs = rf$trees[[treeID]]$split_clusterIDs
    for (cluster in clusterIDs) {
      if (!is.na(cluster)) {
        split_counts[cluster] = split_counts[cluster] + 1
      }
    }
  }
  res$split_group_counts = split_counts
  
  res$modules = modules
  res$causal_modules = network$causal_modules
  
  return(res)
}

# for (i in 1:1) {
for (i in 1:nrow(scenarios)) {
  # read scenario
  scenario = scenarios[i,]
  print(scenario)
  
  for (method in c("LDA", "PCA", "logridge1")) {
    print(method)
    # for (method in c("LDA")) {
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
    
    splitmethod = method
    
    # prederr_res = lapply(dat, pred_err)
    prederr_res = mclapply(dat, pred_err, mc.cores = 60)
    
    saveroot = paste0(
      "./results/prederrres_",
      splitmethod,
      "_ndm", scenario$n_disease_modules,
      "_ab", scenario$average_beta,
      ".Rdata"
    )
    save(prederr_res, file = saveroot)
  }
}
