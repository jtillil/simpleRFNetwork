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

## run rf on rna-seq data
predlabel_micro = microarray$pheno
tcgares = list()

tcgares$prederr_micro_LDA = c()
for (i in 1:100) {
  print("LDA")
  print(i)
  rf = simpleRFNetwork(
    pheno ~ .,
    data = rna_seq,
    num_trees=500,
    splitobject="module",
    splitmethod="LDA",
    varselection="none",
    mtry="root",
    varclusters = modules,
    seed = as.integer(i),
    num_threads = 60
  )
  pred_micro_LDA = rf$predict(as.matrix(microarray[, -1]))
  prederr_micro_LDA = sum(pred_micro_LDA != predlabel_micro) / 283
  tcgares$prederr_micro_LDA = c(tcgares$prederr_micro_LDA, prederr_micro_LDA)
}

tcgares$prederr_micro_Ridge = c()
for (i in 1:100) {
  print("Ridge")
  print(i)
  rf = simpleRFNetwork(
    pheno ~ .,
    data = rna_seq,
    num_trees=500,
    splitobject="module",
    splitmethod="logridge1",
    varselection="none",
    mtry="root",
    varclusters = modules,
    seed = as.integer(i),
    num_threads = 60
  )
  pred_micro_Ridge = rf$predict(as.matrix(microarray[, -1]))
  prederr_micro_Ridge = sum(pred_micro_Ridge != predlabel_micro) / 283
  tcgares$prederr_micro_Ridge = c(tcgares$prederr_micro_Ridge, prederr_micro_Ridge)
}

tcgares$prederr_micro_PCA = c()
for (i in 1:100) {
  print("PCA")
  print(i)
  rf = simpleRFNetwork(
    pheno ~ .,
    data = rna_seq,
    num_trees=500,
    splitobject="module",
    splitmethod="PCA",
    varselection="none",
    mtry="root",
    varclusters = modules,
    seed = as.integer(i),
    num_threads = 60
  )
  pred_micro_PCA = rf$predict(as.matrix(microarray[, -1]))
  prederr_micro_PCA = sum(pred_micro_PCA != predlabel_micro) / 283
  tcgares$prederr_micro_PCA = c(tcgares$prederr_micro_PCA, prederr_micro_PCA)
}

saveroot = paste0(
  "./results/tcgares_micro.Rdata"
)
tcgares_micro = tcgares
save(tcgares_micro, file = saveroot)