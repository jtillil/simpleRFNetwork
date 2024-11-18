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
library(SIS)

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

## get hunames
tcganames = colnames(microarray[, -1])
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

#### select genes
pvals = numeric(ncol(microarray) - 1)
# alpha = 0.01

for (i in 1:length(signif)) {
  print(i)
  pvals[i] = summary(glm(pheno ~ ., data = microarray[, c(1, i+1)], family = "binomial"))$coefficients[2, "Pr(>|z|)"]
}

signif = pvals < 5e-2
sum(signif)
sum(hunames %in% ((1:(ncol(microarray) - 1))[signif]))

sis_selected <- SIS(as.matrix(microarray[,-1]), as.numeric(microarray[,1])-1, family = "binomial")
sel = c(sis_selected$sis.ix0, sis_selected$ix)
sum(hunames %in% sel)
hunames[hunames %in% sel]

hunames[hunames %in% (1:14167)[signif]]


#### build modules
igraph_network = upgrade_graph(tcga_breast_pr$network)
igraph_network <- induced_subgraph(igraph_network, vids = (1:(ncol(microarray) - 1))[signif])
set.seed(1)
igraph_modules = cluster_louvain(igraph_network, weights = NULL, resolution = 5)
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
# for (i in length(modules):1) {
#   if (length(modules[[i]]) < 10 ) {
#     modules = modules[-i]
#   }
# }
sort(lengths(modules))
sum(lengths(modules))
length(modules)

# identify modules containing the hunames genes

containing_modules = c()
for (i in 1:length(modules)) {
  if (any(hunames %in% modules[[i]])) {
    containing_modules = c(containing_modules, i)
  }
}
containing_modules
lengths(modules)[containing_modules]
for (mod in containing_modules) {
  print(mod)
  print(sum(hunames %in% modules[[mod]]))
  print("")
}
