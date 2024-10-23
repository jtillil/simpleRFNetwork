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

## generate hunames
tcganames = colnames(rna_seq[, -1])
hunames = (1:14167)[tcganames %in% c(
  "PGR",
  "AR",
  "WDR19",
  "GATA3",
  "GREB1",
  "ESR1",
  "CA12",
  "SLC39A6",
  "SCUBE2",
  "C6ORF97",
  "DNALI1",
  "SERPINA11",
  "ZMYND10",
  "FGD3",
  "ABAT",
  "IL6ST",
  "PREX1",
  "THSD4",
  "B3GNT5",
  "PSAT1",
  "MAPT"
)]

for (resolution in c(5, 10, 15)) {
  ## build modules
  igraph_network = tcga_breast_pr$network
  set.seed(1)
  igraph_modules = cluster_louvain(igraph_network, weights = NULL, resolution = resolution)
  igraph_modules$membership[hunames] = max(igraph_modules$membership) + 1
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
  tcgadat_rnaseq$data = rna_seq
  tcgadat_rnaseq$modules = modules
  
  tcgadat_micro = list()
  tcgadat_micro$data = microarray
  tcgadat_micro$modules = modules
  
  #### start calculations
  method = "PCA"
  importance = "permutation"
  n_iterations = 20
  
  saveroot = paste0(
    "./results/tcga",
    "_", method,
    "_rnaseq",
    "_res", resolution,
    ".Rdata"
  )
  borutares = list()
  save(borutares, file = saveroot)
  boruta_TCGA(tcgadat_rnaseq, 1, method, importance, 500, 60, n_iterations, 1, saveroot)
  
  saveroot = paste0(
    "./results/tcga",
    "_", method,
    "_micro",
    "_res", resolution,
    ".Rdata"
  )
  borutares = list()
  save(borutares, file = saveroot)
  boruta_TCGA(tcgadat_micro, 1, method, importance, 500, 60, n_iterations, 1, saveroot)
}


