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
igraph_modules = cluster_louvain(igraph_network, weights = NULL, resolution = 15)
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
length(modules)

tcgadat_rnaseq = list()
tcgadat_rnaseq$dat = rna_seq
tcgadat_rnaseq$modules = modules

tcgadat_micro = list()
tcgadat_micro$dat = microarray
tcgadat_micro$modules = modules

## run rf on microarray data
# predlabel_rnaseq = rna_seq$pheno
# tcgares = list()
# 
# tcgares$prederr_rnaseq_LDA = c()
# for (i in 1:100) {
#   print("LDA")
#   print(i)
#   rf = simpleRFNetwork(
#     pheno ~ .,
#     data = microarray,
#     num_trees=500,
#     splitobject="module",
#     splitmethod="LDA",
#     varselection="none",
#     mtry="root",
#     varclusters = modules,
#     seed = as.integer(i),
#     num_threads = 2
#   )
#   pred_rnaseq_LDA = rf$predict(as.matrix(rna_seq[, -1]))
#   prederr_rnaseq_LDA = sum(pred_rnaseq_LDA != predlabel_rnaseq) / 284
#   tcgares$prederr_rnaseq_LDA = c(tcgares$prederr_rnaseq_LDA, prederr_rnaseq_LDA)
# }
# 
# tcgares$prederr_rnaseq_Ridge = c()
# for (i in 1:100) {
#   print("Ridge")
#   print(i)
#   rf = simpleRFNetwork(
#     pheno ~ .,
#     data = microarray,
#     num_trees=500,
#     splitobject="module",
#     splitmethod="logridge1",
#     varselection="none",
#     mtry="root",
#     varclusters = modules,
#     seed = as.integer(i),
#     num_threads = 2
#   )
#   pred_rnaseq_Ridge = rf$predict(as.matrix(rna_seq[, -1]))
#   prederr_rnaseq_Ridge = sum(pred_rnaseq_Ridge != predlabel_rnaseq) / 284
#   tcgares$prederr_rnaseq_Ridge = c(tcgares$prederr_rnaseq_Ridge, prederr_rnaseq_Ridge)
# }
# 
# tcgares$prederr_rnaseq_PCA = c()
# for (i in 1:100) {
#   print("PCA")
#   print(i)
#   rf = simpleRFNetwork(
#     pheno ~ .,
#     data = microarray,
#     num_trees=500,
#     splitobject="module",
#     splitmethod="PCA",
#     varselection="none",
#     mtry="root",
#     varclusters = modules,
#     seed = as.integer(i),
#     num_threads = 2
#   )
#   pred_rnaseq_PCA = rf$predict(as.matrix(rna_seq[, -1]))
#   prederr_rnaseq_PCA = sum(pred_rnaseq_PCA != predlabel_rnaseq) / 284
#   tcgares$prederr_rnaseq_PCA = c(tcgares$prederr_rnaseq_PCA, prederr_rnaseq_PCA)
# }
# 
# saveroot = paste0(
#   "./results/tcgares_rnaseq.Rdata"
# )
# tcgares_rnaseq = tcgares
# save(tcgares_rnaseq, file = saveroot)

# ## run rf on rna-seq data## run rf on rna-seq data
# predlabel_micro = microarray$pheno
# tcgares = list()
# 
# tcgares$prederr_micro_LDA = c()
# for (i in 1:100) {
#   print("LDA")
#   print(i)
#   rf = simpleRFNetwork(
#     pheno ~ .,
#     data = rna_seq,
#     num_trees=500,
#     splitobject="module",
#     splitmethod="LDA",
#     varselection="none",
#     mtry="root",
#     varclusters = modules,
#     seed = as.integer(i),
#     num_threads = 2
#   )
#   pred_micro_LDA = rf$predict(as.matrix(microarray[, -1]))
#   prederr_micro_LDA = sum(pred_micro_LDA != predlabel_micro) / 283
#   tcgares$prederr_micro_LDA = c(tcgares$prederr_micro_LDA, prederr_micro_LDA)
# }
# 
# tcgares$prederr_micro_Ridge = c()
# for (i in 1:100) {
#   print("Ridge")
#   print(i)
#   rf = simpleRFNetwork(
#     pheno ~ .,
#     data = rna_seq,
#     num_trees=500,
#     splitobject="module",
#     splitmethod="logridge1",
#     varselection="none",
#     mtry="root",
#     varclusters = modules,
#     seed = as.integer(i),
#     num_threads = 2
#   )
#   pred_micro_Ridge = rf$predict(as.matrix(microarray[, -1]))
#   prederr_micro_Ridge = sum(pred_micro_Ridge != predlabel_micro) / 283
#   tcgares$prederr_micro_Ridge = c(tcgares$prederr_micro_Ridge, prederr_micro_Ridge)
# }
# 
# tcgares$prederr_micro_PCA = c()
# for (i in 1:100) {
#   print("PCA")
#   print(i)
#   rf = simpleRFNetwork(
#     pheno ~ .,
#     data = rna_seq,
#     num_trees=500,
#     splitobject="module",
#     splitmethod="PCA",
#     varselection="none",
#     mtry="root",
#     varclusters = modules,
#     seed = as.integer(i),
#     num_threads = 2
#   )
#   pred_micro_PCA = rf$predict(as.matrix(microarray[, -1]))
#   prederr_micro_PCA = sum(pred_micro_PCA != predlabel_micro) / 283
#   tcgares$prederr_micro_PCA = c(tcgares$prederr_micro_PCA, prederr_micro_PCA)
# }
# 
# saveroot = paste0(
#   "./results/tcgares_micro.Rdata"
# )
# tcgares_micro = tcgares
# save(tcgares_micro, file = saveroot)

#### ranger

devtools::install_github("silkeszy/Pomona")

library(ranger)
library(Pomona)

rangerrf = ranger(
  dependent.variable.name = "pheno",
  data = rna_seq,
  num.trees = 500,
  seed = 10
)
rangerrf$predict(as.matrix(microarray[, -1]))

rangerpred = predict(rangerrf, microarray[,-1])$predictions

set.seed(1)
vita_rnaseq = var.sel.vita(
  x = rna_seq[, -1],
  y = rna_seq[, 1],
  type = "classification"
)
sum(vita_rnaseq$info$selected)
length(vita_rnaseq$info$selected)

tcganames = colnames(rna_seq[, -1])
brca1 = (1:14167)[tcganames == "BRCA1"]
brca2 = (1:14167)[tcganames == "BRCA2"]
tcganames[brca1]
vita_rnaseq$info$selected[brca1]
tcganames[brca2]
vita_rnaseq$info$selected[brca2]

set.seed(1)
vita_microarray = var.sel.vita(
  x = microarray[, -1],
  y = microarray[, 1],
  type = "classification"
)
sum(vita_microarray$info$selected)

tcganames[12442]
vita_microarray$info$selected[12442]
tcganames[9898]
vita_microarray$info$selected[9898]

prederr_micro_ranger = c()
for (i in 1:100) {
  print("ranger")
  print(i)
  
  rangerrf = ranger(
    dependent.variable.name = "pheno",
    data = rna_seq,
    num.trees = 500,
    seed = i
  )
  pred_micro_ranger = predict(rangerrf, microarray[, -1])$predictions
  prederr = sum(pred_micro_ranger != predlabel_micro) / 283
  prederr_micro_ranger = c(prederr_micro_ranger, prederr)
}

prederr_rnaseq_ranger = c()
for (i in 1:100) {
  print("ranger")
  print(i)
  
  rangerrf = ranger(
    dependent.variable.name = "pheno",
    data = microarray,
    num.trees = 500,
    seed = i
  )
  pred_rnaseq_ranger = predict(rangerrf, rna_seq[, -1])$predictions
  prederr = sum(pred_rnaseq_ranger != predlabel_rnaseq) / 284
  prederr_rnaseq_ranger = c(prederr_rnaseq_ranger, prederr)
}


#### plotting

library(ggplot2)

tcgares_dat = data.frame(
  evaldata = c(rep("Rna-Seq", 400), rep("Microarray", 400)),
  prederr = c(prederr_rnaseq_ranger, tcgares_rnaseq$prederr_rnaseq_LDA, tcgares_rnaseq$prederr_rnaseq_Ridge, tcgares_rnaseq$prederr_rnaseq_PCA, prederr_micro_ranger, tcgares_micro$prederr_micro_LDA, tcgares_micro$prederr_micro_Ridge, tcgares_micro$prederr_micro_PCA),
  method = factor(c(rep("RF", 100), rep("Group LDA", 100), rep("Group Ridge", 100), rep("Group PCA", 100), rep("RF", 100), rep("Group LDA", 100), rep("Group Ridge", 100), rep("Group PCA", 100)), levels = c("RF", "Group LDA", "Group Ridge", "Group PCA"), ordered = TRUE)
)
colnames(tcgares_dat) = c("evaldata", "prederr", "Method")

ggplot(tcgares_dat, aes(x = evaldata, y = prederr, fill = Method)) +
  geom_boxplot() +
  theme_bw() +
  ylab("Prediction error") +
  xlab("Evaluation dataset") +
  scale_fill_hue()

# scale_fill_discrete(name = "Splitmethod", labels = c("LDA", "Ridge", "PCA"))

ggsave("box_prederr_micro_rnaseq.pdf", width = 7, height = 4)

#### Comparison to Hu and Scymczak detected genes

# get detected gene positions
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
tcganames[hunames]
# hunames_network = (1:14167)[tcganames %in% c(
#   
# )]
# hunames_standard = (1:14167)[tcganames %in% c(
#   
# )]

# identify modules containing the detected genes

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

#### TCGA module selection results

containing_modules = c()
for (i in 1:length(modules)) {
  if (brca1 %in% modules[[i]] | brca2 %in% modules[[i]]) {
    containing_modules = c(containing_modules, i)
  }
}
containing_modules
lengths(modules)[containing_modules]

selected = borutares[[1]]$aggregated_classifications == 1
selected_modules = (1:length(modules))[selected]
selected_modules
lengths(modules)[selected_modules]

binomres = borutares[[1]]$first_binomresults
binomres
binomres[[37]]

vim = borutares[[1]]$first_vim
vimsum = colsums(vim)
vimsum[1:(length(vim)/2)][selected_modules]
vim[,1:(ncol(vim)/2)][,selected_modules]

vimsum[1:(length(vim)/2)][containing_modules]
vim[,1:(ncol(vim)/2)][,containing_modules]

## REDO EVERYTHING WITH ACTUALLY ALL MODULES

#### resolution = 15

## good module
# 37, length 134
## LDA micro
# 71
# b = 0
## LDA rnaseq
# -
# b = 0
## Ridge micro
# -
# b = 0
## Ridge rnaseq
# -
# b = 0
## PCA micro
# -
# b = 0
## PCA rnaseq
# 4 64 104 145 241
# b = 0

#### resolution = 10

## good module
# 39, length 189
## LDA micro
# -
# b = 0
## LDA rnaseq
# 22
# b = 0
## Ridge micro
# 47 59 61
# b = 0
## Ridge rnaseq
# 43
# b = 0
## PCA micro
# 53 55
# b = 0
## PCA rnaseq
# 12
# b = 0

#### resolution = 5

## good module
# 37, length 189
## LDA micro
# 30  34 117 121
# b = 0
## LDA rnaseq
# -
# b = 0
## Ridge micro
# -
# b = 0
## Ridge rnaseq
# -
# b = 0
## PCA micro
# 13  71 100 158 187
# b = 0
## PCA rnaseq
# 59 142  77  80  42  30  53  41  72  55  18  17  16
# b = 0

#### resolution = 4

## good module
# 28, length 479
## LDA micro
# -
# b = 0
## LDA rnaseq
# 22
# b = 0
## Ridge micro
# 47 59 61
# b = 0
## Ridge rnaseq
# 43
# b = 0
## PCA micro
# 53 55
# b = 0
## PCA rnaseq
# 12
# b = 0

