setwd(getSrcDirectory(function(){})[1])
source("./source_files.R")

library(igraph)

library(simpleRFNetwork)
library(SeqNet)
library(parallel)
library(pracma)
library(tictoc)
library(Rfast)
library(matrixcalc)

library(ridge)
library(glmnet)

# ## load and extract tcga data
# datroot = "./data/tcga_breast_pr.rdata"
# load(datroot)
# 
# pheno = tcga_breast_pr$microarray$pheno
# pheno[pheno == "Positive"] = 1
# pheno[pheno == "Negative"] = 0
# pheno = as.numeric(pheno)
# microarray = data.frame(pheno = as.factor(pheno))
# microarray = cbind(microarray, tcga_breast_pr$microarray$geno)
# 
# pheno = tcga_breast_pr$rna_seq$pheno
# pheno[pheno == "Positive"] = 1
# pheno[pheno == "Negative"] = 0
# pheno = as.numeric(pheno)
# rna_seq = data.frame(pheno = as.factor(pheno))
# rna_seq = cbind(rna_seq, tcga_breast_pr$rna_seq$geno)
# 
# #### select genes
# tcganames = colnames(microarray[, -1])
# pvals_micro = numeric(length(tcganames))
# pvals_rna = numeric(length(tcganames))
# 
# for (i in 1:length(tcganames)) {
#   print(i)
#   pvals_micro[i] = summary(glm(pheno ~ ., data = microarray[, c(1, i+1)], family = "binomial"))$coefficients[2, "Pr(>|z|)"]
#   pvals_rna[i] = summary(glm(pheno ~ ., data = rna_seq[, c(1, i+1)], family = "binomial"))$coefficients[2, "Pr(>|z|)"]
# }
# 
# signif_micro = pvals_micro < 1e-5
# signif_rna = pvals_rna < 1e-5

# for (resolution in c(5, 10, 15)) {
for (resolution in c(5)) {
  # #### build modules
  # igraph_network = upgrade_graph(tcga_breast_pr$network)
  # igraph_network <- induced_subgraph(igraph_network, vids = (1:(length(tcganames)))[signif_rna])
  # set.seed(1)
  # igraph_modules = cluster_louvain(igraph_network, weights = NULL, resolution = resolution)
  # sizes(igraph_modules)
  # modules = list()
  # for (i in 1:max(igraph_modules$membership)) {
  #   modules[[i]] = numeric()
  # }
  # for (i in 1:length(igraph_modules$membership)) {
  #   membership = igraph_modules$membership[i]
  #   modules[[membership]] = c(modules[[membership]], i)
  # }
  # for (i in length(modules):1) {
  #   if (length(modules[[i]]) < 2 ) {
  #     modules = modules[-i]
  #   }
  # }
  # 
  # tcgadat_rnaseq = list()
  # tcgadat_rnaseq$data = rna_seq
  # tcgadat_rnaseq$modules = modules
  
  #### load modules
  load(paste0("./data/tcga_modules_rnaseq_res", resolution, ".Rdata"))
  
  #### start calculations
  method = "LDA"
  importance = "permutation"
  n_iterations = 20
  
  saveroot = paste0(
    "./results/tcga",
    "_", method,
    "_rnaseq",
    "_all_modulesgt1",
    "_res", resolution,
    ".Rdata"
  )
  borutares = list()
  save(borutares, file = saveroot)
  boruta_TCGA(tcgadat_rnaseq, 1, method, importance, 500, 30, n_iterations, 1, saveroot)
  
  # #### build modules
  # igraph_network = upgrade_graph(tcga_breast_pr$network)
  # igraph_network <- induced_subgraph(igraph_network, vids = (1:(length(tcganames)))[signif_micro])
  # set.seed(1)
  # igraph_modules = cluster_louvain(igraph_network, weights = NULL, resolution = resolution)
  # sizes(igraph_modules)
  # modules = list()
  # for (i in 1:max(igraph_modules$membership)) {
  #   modules[[i]] = numeric()
  # }
  # for (i in 1:length(igraph_modules$membership)) {
  #   membership = igraph_modules$membership[i]
  #   modules[[membership]] = c(modules[[membership]], i)
  # }
  # for (i in length(modules):1) {
  #   if (length(modules[[i]]) < 2 ) {
  #     modules = modules[-i]
  #   }
  # }
  # 
  # tcgadat_micro = list()
  # tcgadat_micro$data = microarray
  # tcgadat_micro$modules = modules
  
  #### load modules
  load(paste0("./data/tcga_modules_micro_res", resolution, ".Rdata"))
  
  saveroot = paste0(
    "./results/tcga",
    "_", method,
    "_micro",
    "_all_modulesgt1",
    "_res", resolution,
    ".Rdata"
  )
  borutares = list()
  save(borutares, file = saveroot)
  boruta_TCGA(tcgadat_micro, 1, method, importance, 500, 30, n_iterations, 1, saveroot)
}


