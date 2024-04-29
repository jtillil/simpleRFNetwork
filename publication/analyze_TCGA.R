setwd(getSrcDirectory(function(){})[1])
source("./source_files.R")

library(igraph)

## load and extract tcga data
datroot = "./data/tcga_breast_pr.rdata"
load(datroot)

pheno = tcga_breast_pr$microarray$pheno
pheno[pheno == "Positive"] = 1
pheno[pheno == "Negative"] = 0
pheno = as.numeric(pheno)
microarray = data.frame(pheno = as.factor(pheno))
microarray = cbind(microarray, tcga_breast_pr$microarray$geno)

# rna_seq = data.frame(pheno = as.factor(tcga_breast_pr$rna_seq$pheno))
# rna_seq = cbind(rna_seq, tcga_breast_pr$rna_seq$geno)

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
lengths(modules)

## run rf
rf = simpleRFNetwork(
  pheno ~ .,
  data = microarray,
  num_trees=2,
  num_threads=1,
  splitobject="module",
  splitmethod="LDA",
  varselection="none",
  mtry="root",
  varclusters = modules,
  seed = 1L
)
