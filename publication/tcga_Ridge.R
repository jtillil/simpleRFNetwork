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

## load and extract tcga data
datroot = "./data/tcga_breast_pr.rdata"
load(datroot)

pheno = tcga_breast_pr$microarray$pheno
pheno[pheno == "Positive"] = 1
pheno[pheno == "Negative"] = 0
pheno = as.numeric(pheno)
microarray = data.frame(pheno = as.factor(pheno))
microarray = cbind(microarray, tcga_breast_pr$microarray$geno)

pheno = tcga_breast_pr$rna_seq$pheno
pheno[pheno == "Positive"] = 1
pheno[pheno == "Negative"] = 0
pheno = as.numeric(pheno)
rna_seq = data.frame(pheno = as.factor(pheno))
rna_seq = cbind(rna_seq, tcga_breast_pr$rna_seq$geno)

## build modules
igraph_network = tcga_breast_pr$network
set.seed(1)
igraph_modules = cluster_louvain(igraph_network, weights = NULL, resolution = 4)
sizes(igraph_modules)
# sum(sizes(igraph_modules) >= 10)
# mean(sizes(igraph_modules))
modules = list()
for (i in 1:max(igraph_modules$membership)) {
  modules[[i]] = numeric()
}
for (i in 1:length(igraph_modules$membership)) {
  membership = igraph_modules$membership[i]
  modules[[membership]] = c(modules[[membership]], i)
}
for (i in length(modules):1) {
  if (length(modules[[i]]) < 10 ) {
    modules = modules[-i]
  }
}
lengths(modules)

tcgadat_rnaseq = list()
tcgadat_rnaseq$dat = rna_seq
tcgadat_rnaseq$modules = modules

tcgadat_micro = list()
tcgadat_micro$dat = microarray
tcgadat_micro$modules = modules

#### start calculations
method = "logridge1"
importance = "permutation"
n_iterations = 20

saveroot = paste0(
  "./results/tcga",
  "_", method,
  "_rnaseq",
  ".Rdata"
)
borutares = list()
save(borutares, file = saveroot)
boruta_TCGA(tcgadat_rnaseq, 1, method, importance, 500, 60, n_iterations, 1, saveroot)

saveroot = paste0(
  "./results/tcga",
  "_", method,
  "_micro",
  ".Rdata"
)
borutares = list()
save(borutares, file = saveroot)
boruta_TCGA(tcgadat_micro, 1, method, importance, 500, 60, n_iterations, 1, saveroot)

