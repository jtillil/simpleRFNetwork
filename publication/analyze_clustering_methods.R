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
igraph_network = upgrade_graph(tcga_breast_pr$network)
set.seed(1)
igraph_modules = cluster_louvain(igraph_network, weights = NULL, resolution = 10)
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
sort(lengths(modules))

tcgadat_rnaseq = list()
tcgadat_rnaseq$dat = rna_seq
tcgadat_rnaseq$modules = modules

tcgadat_micro = list()
tcgadat_micro$dat = microarray
tcgadat_micro$modules = modules

################################################################################
################################################################################
################################################################################

ign = igraph_network
comp = components(igraph_network)
ign = largest_component(ign)
ign = simplify(ign)

igncov1 = cova(as.matrix(microarray[,-1]), large=TRUE)
igncov2 = cova(as.matrix(rna_seq[,-1]), large=TRUE)
igncov = 0.5*(igncov1 + igncov2)
igncovabs = abs(igncov)
rm(igncov1)
rm(igncov2)

################################################################################
################################################################################
################################################################################

## HU genes
# get detected gene positions
tcganames = colnames(rna_seq[, -1]) # [comp$membership==1]
hunames = (1:14157)[tcganames %in% c(
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

hucov = igncov[hunames, hunames]

clustering_methods = c(
  # "cluster_edge_betweenness" = cluster_edge_betweenness, # long
    # long calc times
    # no arguments
  # "cluster_fast_greedy" = cluster_fast_greedy, # no multi-edges
    # few large communities
    # no arguments
  # "cluster_fluid_communities" = cluster_fluid_communities # only simple graphs
    # high no.of.communities ==> similarly sized smaller modules
    # no weights possible
    # hunames distributed
  # "cluster_infomap" = cluster_infomap,
  # "cluster_label_prop" = cluster_label_prop,
    # can set initial labels
  # "cluster_leading_eigen" = cluster_leading_eigen,
  "cluster_leiden" = cluster_leiden
    # can set initial membership
  # "cluster_louvain" = cluster_louvain,
  # "cluster_optimal" = cluster_optimal, # too many columns
    # NP-complete, runs in exponential time
  # "cluster_spinglass" = cluster_spinglass, # has to be connected
  # "cluster_walktrap" = cluster_walktrap
    # few large communities
)

#### graph with weights ########################################################

edges = E(ign)
weights = c()
for (eID in 1:length(edges)) {
  print(eID)
  weights = c(weights, igncovabs[edges[eID][[1]]$tid, edges[eID][[1]]$hid])
}

#### LOUVAIN init ##############################################################

set.seed(1)
louvain_modules = cluster_louvain(ign, resolution = 10)
print(sort(sizes(louvain_modules)))
louvain_membership = louvain_modules$membership
louvain_membership[hunames] = max(louvain_membership) + 1

#### LEIDEN further optimization ###############################################

leiden_modules = cluster_leiden(ign, initial_membership = louvain_membership, resolution = 1e-2)
print(sort(sizes(leiden_modules)))
print(max(leiden_modules$membership))
# some very large modules, similarly many as LOUVAIN

#### LABEL_PROP further optimization ###########################################

init_hunames = 0:(sum(comp$membership==1) - 1)
init_hunames[hunames] = hunames[1] - 1
label_modules = cluster_label_prop(ign, initial = init_hunames, fixed = (1:sum(comp$membership==1) %in% hunames))
print(sort(sizes(label_modules)))
print(max(label_modules$membership))
# single module

################################################################################

for (i in 1:length(clustering_methods)) {
  print(names(clustering_methods)[i])
  # igraph_modules = clustering_methods[[i]](ign, no.of.communities = 50)
  igraph_modules = clustering_methods[[i]](ign, resolution = 0.01)
  # igraph_modules = clustering_methods[[i]](ign)
  print(sum(sizes(igraph_modules)))
  print(sizes(igraph_modules))
}

membership = louvain_membership
modules = list()
for (i in 1:max(membership)) {
  modules[[i]] = numeric()
}
for (i in 1:length(membership)) {
  moduleID = membership[i]
  modules[[moduleID]] = c(modules[[moduleID]], i)
}

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
  # print("")
}
print(sort(lengths(modules)))

