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

#### select genes
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
pvals_micro = numeric(length(tcganames))
pvals_rna = numeric(length(tcganames))

for (i in 1:length(tcganames)) {
  print(i)
  pvals_micro[i] = summary(glm(pheno ~ ., data = microarray[, c(1, i+1)], family = "binomial"))$coefficients[2, "Pr(>|z|)"]
  pvals_rna[i] = summary(glm(pheno ~ ., data = rna_seq[, c(1, i+1)], family = "binomial"))$coefficients[2, "Pr(>|z|)"]
}

signif_micro = pvals_micro < 1e-5
signif_rna = pvals_rna < 1e-5

#### build modules
resolution = 5

igraph_network = upgrade_graph(tcga_breast_pr$network)
igraph_network <- induced_subgraph(igraph_network, vids = (1:(length(tcganames)))[signif_rna])
set.seed(1)
igraph_modules = cluster_louvain(igraph_network, weights = NULL, resolution = resolution)
sizes(igraph_modules)
modules = list()
for (i in 1:max(igraph_modules$membership)) {
  modules[[i]] = numeric()
}
for (i in 1:length(igraph_modules$membership)) {
  membership = igraph_modules$membership[i]
  modules[[membership]] = c(modules[[membership]], i)
}
for (i in length(modules):1) {
  if (length(modules[[i]]) < 2 ) {
    modules = modules[-i]
  }
}
modules_rna = modules

## build modules
# igraph_network = tcga_breast_pr$network
# set.seed(1)
# igraph_modules = cluster_louvain(igraph_network, weights = NULL, resolution = 15)
# sizes(igraph_modules)
# # sum(sizes(igraph_modules) >= 10)
# # mean(sizes(igraph_modules))
# modules = list()
# for (i in 1:max(igraph_modules$membership)) {
#   modules[[i]] = numeric()
# }
# for (i in 1:length(igraph_modules$membership)) {
#   membership = igraph_modules$membership[i]
#   modules[[membership]] = c(modules[[membership]], i)
# }
# for (i in length(modules):1) {
#   if (length(modules[[i]]) < 10 ) {
#     modules = modules[-i]
#   }
# }
# sort(lengths(modules))
# length(modules)

tcgadat_rnaseq = list()
tcgadat_rnaseq$dat = rna_seq
tcgadat_rnaseq$modules = modules_rna

igraph_network = upgrade_graph(tcga_breast_pr$network)
igraph_network <- induced_subgraph(igraph_network, vids = (1:(length(tcganames)))[signif_micro])
set.seed(1)
igraph_modules = cluster_louvain(igraph_network, weights = NULL, resolution = resolution)
sizes(igraph_modules)
modules = list()
for (i in 1:max(igraph_modules$membership)) {
  modules[[i]] = numeric()
}
for (i in 1:length(igraph_modules$membership)) {
  membership = igraph_modules$membership[i]
  modules[[membership]] = c(modules[[membership]], i)
}
for (i in length(modules):1) {
  if (length(modules[[i]]) < 2 ) {
    modules = modules[-i]
  }
}
modules_micro = modules

tcgadat_micro = list()
tcgadat_micro$dat = microarray
tcgadat_micro$modules = modules_micro

## run rf on microarray data
predlabel_rnaseq = rna_seq$pheno
tcgares = list()

tcgares$prederr_rnaseq_LDA = c()
for (i in 1:100) {
  print("LDA")
  print(i)
  rf = simpleRFNetwork(
    pheno ~ .,
    data = microarray,
    num_trees=500,
    splitobject="module",
    splitmethod="LDA",
    varselection="none",
    mtry="root",
    varclusters = modules,
    seed = as.integer(i),
    num_threads = 2
  )
  pred_rnaseq_LDA = rf$predict(as.matrix(rna_seq[, -1]))
  prederr_rnaseq_LDA = sum(pred_rnaseq_LDA != predlabel_rnaseq) / 284
  tcgares$prederr_rnaseq_LDA = c(tcgares$prederr_rnaseq_LDA, prederr_rnaseq_LDA)
}

tcgares$prederr_rnaseq_Ridge = c()
for (i in 1:100) {
  print("Ridge")
  print(i)
  rf = simpleRFNetwork(
    pheno ~ .,
    data = microarray,
    num_trees=500,
    splitobject="module",
    splitmethod="logridge1",
    varselection="none",
    mtry="root",
    varclusters = modules,
    seed = as.integer(i),
    num_threads = 2
  )
  pred_rnaseq_Ridge = rf$predict(as.matrix(rna_seq[, -1]))
  prederr_rnaseq_Ridge = sum(pred_rnaseq_Ridge != predlabel_rnaseq) / 284
  tcgares$prederr_rnaseq_Ridge = c(tcgares$prederr_rnaseq_Ridge, prederr_rnaseq_Ridge)
}

tcgares$prederr_rnaseq_PCA = c()
for (i in 1:100) {
  print("PCA")
  print(i)
  rf = simpleRFNetwork(
    pheno ~ .,
    data = microarray,
    num_trees=500,
    splitobject="module",
    splitmethod="PCA",
    varselection="none",
    mtry="root",
    varclusters = modules,
    seed = as.integer(i),
    num_threads = 2
  )
  pred_rnaseq_PCA = rf$predict(as.matrix(rna_seq[, -1]))
  prederr_rnaseq_PCA = sum(pred_rnaseq_PCA != predlabel_rnaseq) / 284
  tcgares$prederr_rnaseq_PCA = c(tcgares$prederr_rnaseq_PCA, prederr_rnaseq_PCA)
}

saveroot = paste0(
  "./results/tcgares_rnaseq.Rdata"
)
tcgares_rnaseq = tcgares
save(tcgares_rnaseq, file = saveroot)

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
install.packages("ranger")

library(ranger)
library(Pomona)

# rangerrf = ranger(
#   dependent.variable.name = "pheno",
#   data = rna_seq,
#   num.trees = 500,
#   seed = 10
# )
# rangerrf$predict(as.matrix(microarray[, -1]))
# 
# rangerpred = predict(rangerrf, microarray[,-1])$predictions

## vita selection

set.seed(1)
vita_rnaseq = var.sel.vita(
  x = rna_seq[, -1],
  y = rna_seq[, 1],
  type = "classification"
)
sum(vita_rnaseq$info$selected)
length(vita_rnaseq$info$selected)
vita_rnaseq$info$selected[hunames]

set.seed(1)
vita_microarray = var.sel.vita(
  x = microarray[, -1],
  y = microarray[, 1],
  type = "classification"
)
sum(vita_microarray$info$selected)
vita_microarray$info$selected[hunames]

## network prediction

prederr_micro_ranger = c()
for (i in 1:100) {
  print("ranger")
  print(i)
  
  rangerrf = ranger(
    dependent.variable.name = "pheno",
    data = rna_seq[, c(T, signif_rna)],
    num.trees = 500,
    seed = i
  )
  pred_micro_ranger = predict(rangerrf, microarray[, -1][, signif_rna])$predictions
  prederr = sum(pred_micro_ranger != microarray[, 1]) / 283
  prederr_micro_ranger = c(prederr_micro_ranger, prederr)
}

prederr_rnaseq_ranger = c()
for (i in 1:100) {
  print("ranger")
  print(i)
  
  rangerrf = ranger(
    dependent.variable.name = "pheno",
    data = microarray[, c(T, signif_micro)],
    num.trees = 500,
    seed = i
  )
  pred_rnaseq_ranger = predict(rangerrf, rna_seq[, -1][, signif_micro])$predictions
  prederr = sum(pred_rnaseq_ranger != rna_seq[, 1]) / 284
  prederr_rnaseq_ranger = c(prederr_rnaseq_ranger, prederr)
}

#### group lasso

## network prediction

## RNAseq evaluation dataset

# read X and y
X = as.matrix(microarray[, -1][, signif_micro])
y = 2*as.numeric(microarray[, 1]) - 3

preddat = as.matrix(rna_seq[, -1][, signif_micro])
predlabels = 2*as.numeric(rna_seq[, 1]) - 3

# create var and group vectors
var_gglasso = c()
group_gglasso = c()
for (nmodule in 1:length(modules_micro)) {
  for (variable in modules_micro[[nmodule]]) {
    var_gglasso = c(var_gglasso, variable)
    group_gglasso = c(group_gglasso, nmodule)
  }
}

# order group (for gglasso)
ord <- order(group_gglasso)
groupord <- group_gglasso[ord]
# order var according to group
varord <- var_gglasso[ord]

# transform group to have consecutive numbers (for gglasso)
groupb <- cumsum(!duplicated(groupord))

# new data
Xb <- X[, varord]
preddatb = preddat[, varord]

# optimize lambda
gr_cv <- gglasso::cv.gglasso(x=Xb, y=y, group=groupb, 
                             pred.loss="loss", 
                             loss = "logit",
                             # intercept = F, 
                             nfolds=5)

# predict data
pred = predict(gr_cv$gglasso.fit, preddatb, s=gr_cv$lambda.min)
prederr_gglasso_rna = sum(pred != predlabels) / length(pred)

## Microarray evaluation dataset

# read X and y
X = as.matrix(rna_seq[, -1][, signif_rna])
y = 2*as.numeric(rna_seq[, 1]) - 3

preddat = as.matrix(microarray[, -1][, signif_rna])
predlabels = 2*as.numeric(microarray[, 1]) - 3

# create var and group vectors
var_gglasso = c()
group_gglasso = c()
for (nmodule in 1:length(modules_rna)) {
  for (variable in modules_rna[[nmodule]]) {
    var_gglasso = c(var_gglasso, variable)
    group_gglasso = c(group_gglasso, nmodule)
  }
}

# order group (for gglasso)
ord <- order(group_gglasso)
groupord <- group_gglasso[ord]
# order var according to group
varord <- var_gglasso[ord]

# transform group to have consecutive numbers (for gglasso)
groupb <- cumsum(!duplicated(groupord))

# new data
Xb <- X[, varord]
preddatb = preddat[, varord]

# optimize lambda
gr_cv <- gglasso::cv.gglasso(x=Xb, y=y, group=groupb, 
                             pred.loss="loss", 
                             loss = "logit",
                             # intercept = F, 
                             nfolds=5)

# predict data
pred = predict(gr_cv$gglasso.fit, preddatb, s=gr_cv$lambda.1se)
prederr_gglasso_micro = sum(pred != predlabels) / length(pred)

#### plotting

library(ggplot2)

tcgares_dat = data.frame(
  evaldata = c(rep("Rna-Seq", 402), rep("Microarray", 402)),
  prederr = c(prederr_gglasso_rna, prederr_rnaseq_ranger, tcgares_rnaseq$prederr_rnaseq_LDA, tcgares_rnaseq$prederr_rnaseq_Ridge, tcgares_rnaseq$prederr_rnaseq_PCA, prederr_gglasso_micro, prederr_micro_ranger, tcgares_micro$prederr_micro_LDA, tcgares_micro$prederr_micro_Ridge, tcgares_micro$prederr_micro_PCA),
  method = factor(c("Group Lasso", rep("RF", 100), rep("Group LDA", 100), rep("Group Ridge", 100), rep("Group PCA", 100), rep("RF", 100), rep("Group LDA", 100), rep("Group Ridge", 100), rep("Group PCA", 100)), levels = c("RF", "Group LDA", "Group Ridge", "Group PCA"), ordered = TRUE)
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

