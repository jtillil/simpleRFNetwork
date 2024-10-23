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

resolution = 5

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
method = "LDA"
importance = "permutation"
n_iterations = 20

saveroot = paste0(
  "./results/tcga",
  "_", method,
  "_micro",
  "_res", resolution,
  ".Rdata"
)
borutares = list()
save(borutares, file = saveroot)

#### BORUTA SECTION ####

dat = tcgadat_rnaseq
networkID = 1
splitmethod = method
# importance
num_trees = 500
num_threads = 2
num_iterations = n_iterations
seed = 1
# saveroot

# boruta_TCGA <- function(dat, networkID, splitmethod, importance, num_trees, num_threads, num_iterations, seed, saveroot) {

# seed
set.seed(seed)

# initiate
first_binomresults = rep(0, length(dat$modules))
first_bernresults = zeros(num_iterations, length(dat$modules))
first_vim = zeros(num_iterations, 2*length(dat$modules))
iteration = 0
run = 0

print(paste("This is", splitmethod, "Boruta"))
print(saveroot)
print("First run!")

# boruta loop
# while (iteration < num_iterations) {
# advance iteration
iteration = iteration + 1
print(paste("Iteration:", iteration))

# tic()

# add shadow modules
original_and_shadow_modules = dat$modules
for (module in dat$modules) {
  original_and_shadow_modules = c(original_and_shadow_modules, list(module+1000))
}

# add shadow variables
rfdat = dat$data
rfdat = cbind(rfdat, rfdat[,-1])
colnames(rfdat)[-1] = paste0("X", 1:(ncol(rfdat)-1))
# print(colnames(rfdat))
for (col in (((ncol(rfdat)-1)/2)+1):ncol(rfdat)) {
  rfdat[,col] = sample(rfdat[,col])
}

# create model frame
frame_data1 <- model.frame(pheno ~ ., rfdat[,1:15000])
frame_data2 <- model.frame(pheno ~ ., rfdat[,c(1, 15001:ncol(rfdat))])
frame_data <- cbind(frame_data1, frame_data2)
# frame_data <- model.frame(pheno ~ ., dat$data)

# run rf
rf = simpleRFNetwork(pheno ~ .,
                     data = rfdat,
                     num_trees=num_trees,
                     num_threads=num_threads,
                     splitobject="module",
                     splitmethod=splitmethod,
                     varselection="none",
                     mtry="root",
                     varclusters = original_and_shadow_modules,
                     seed = as.integer(iteration + run*num_iterations)
)



# run var importance
varimp = rf$variableImportance(type = importance, num_threads = num_threads)
first_vim[iteration, ] = varimp

# get max var imp of shadow modules
shadow_varimp = varimp[(length(dat$modules)+1):length(varimp)]
max_shadow_varimp = max(shadow_varimp)

# save results
for (i in 1:length(dat$modules) ) {
  if (varimp[i] > max_shadow_varimp) {
    first_binomresults[i] = first_binomresults[i] + 1
    first_bernresults[iteration, i] = 1
  }
}

# clean variables
rm(rf)
rm(rfdat)

# toc()
# }
  
#### STEP 2 ####

  # classify modules
  first_classifications = rep(0, length(dat$modules))
  upcutoff = qbinom(0.975, num_iterations, 0.5)
  downcutoff = qbinom(0.025, num_iterations, 0.5)
  for (i in 1:length(dat$modules) ) {
    if (first_binomresults[i] > upcutoff) {
      first_classifications[i] = 1
    }
    if (first_binomresults[i] < downcutoff) {
      first_classifications[i] = -1
    }
  }
  
  # aggregate classifications
  aggregated_classifications = first_classifications
  
  print("Second run!")
  run = 1
  
  # remove unimportant modules
  updated_modules2 = dat$modules[(aggregated_classifications != -1)]
  
  # initiate
  second_binomresults = rep(0, length(dat$modules))
  second_bernresults = zeros(num_iterations, length(dat$modules))
  second_vim = zeros(num_iterations, 2*length(dat$modules))
  iteration = 0
  
  # prepare
  second_binomresults[(aggregated_classifications == -1)] = NA
  second_bernresults[, (aggregated_classifications == -1)] = NA
  
  if (isempty(updated_modules2)) {
    local_borutares = list(
      networkID = networkID,
      causalmodules = dat$causal_modules,
      aggregated_classifications = aggregated_classifications,
      first_vim = first_vim,
      first_bernresults = first_bernresults,
      first_binomresults = first_binomresults,
      first_classifications = first_classifications,
      updated_modules2 = updated_modules2,
      second_vim = second_vim,
      second_bernresults = second_bernresults,
      second_binomresults = second_binomresults,
      second_classifications = rep(-1, length(dat$modules)),
      updated_modules3 = NULL,
      third_vim = zeros(num_iterations, 2*length(dat$modules)),
      third_bernresults = zeros(num_iterations, length(dat$modules)),
      third_binomresults = rep(0, length(dat$modules)),
      third_classifications = rep(-1, length(dat$modules)),
      updated_modules4 = NULL,
      fourth_vim = zeros(num_iterations, 2*length(dat$modules)),
      fourth_bernresults = zeros(num_iterations, length(dat$modules)),
      fourth_binomresults = rep(0, length(dat$modules)),
      fourth_classifications = rep(-1, length(dat$modules))
    )
    
    # save results
    load(file = saveroot)
    borutares[[networkID]] = local_borutares
    save(borutares, file = saveroot)
    
    return(local_borutares)
  }
  
  # boruta loop
  while (iteration < num_iterations) {
    # advance iteration
    iteration = iteration + 1
    print(paste("Iteration:", iteration))
    
    # tic()
    
    # add shadow modules
    original_and_shadow_modules = updated_modules2
    for (module in updated_modules2) {
      original_and_shadow_modules = c(original_and_shadow_modules, list(module+1000))
    }
    
    # add shadow variables
    rfdat = dat$data
    rfdat = cbind(rfdat, rfdat[,-1])
    colnames(rfdat)[-1] = paste0("X", 1:(ncol(rfdat)-1))
    # print(colnames(rfdat))
    for (col in (((ncol(rfdat)-1)/2)+1):ncol(rfdat)) {
      rfdat[,col] = sample(rfdat[,col])
    }
    
    # run rf
    rf = simpleRFNetwork(pheno ~ .,
                         data = rfdat,
                         num_trees=num_trees,
                         num_threads=num_threads,
                         splitobject="module",
                         splitmethod=splitmethod,
                         varselection="none",
                         mtry="root",
                         varclusters = original_and_shadow_modules,
                         seed = as.integer(iteration + run*num_iterations)
    )
    
    # run var importance
    varimp = rf$variableImportance(type = importance, num_threads = num_threads)
    second_vim[iteration, (aggregated_classifications != -1)] = varimp
    
    # get max var imp of shadow modules
    shadow_varimp = varimp[(length(updated_modules2)+1):length(varimp)]
    max_shadow_varimp = max(shadow_varimp)
    
    # save results
    for (i in 1:length(updated_modules2) ) {
      if (varimp[i] > max_shadow_varimp) {
        second_binomresults[(aggregated_classifications != -1)][i] = second_binomresults[(aggregated_classifications != -1)][i] + 1
        second_bernresults[iteration, (aggregated_classifications != -1)][i] = 1
      }
    }
    
    # clean variables
    rm(rf)
    rm(rfdat)
    
    # toc()
  }
  
  # classify modules
  second_classifications = rep(0, length(dat$modules))
  second_classifications[(aggregated_classifications == -1)] = NA
  upcutoff = qbinom(0.975, num_iterations, 0.5)
  downcutoff = qbinom(0.025, num_iterations, 0.5)
  for (i in 1:length(updated_modules2) ) {
    if (second_binomresults[(aggregated_classifications != -1)][i] > upcutoff) {
      second_classifications[(aggregated_classifications != -1)][i] = 1
    }
    if (second_binomresults[(aggregated_classifications != -1)][i] < downcutoff) {
      second_classifications[(aggregated_classifications != -1)][i] = -1
    }
  }
  
  # aggregate classifications
  interim_classifications = aggregated_classifications
  for (modID in 1:length(updated_modules2)) {
    if (second_classifications[(aggregated_classifications != -1)][modID] == -1) {
      interim_classifications[(aggregated_classifications != -1)][modID] = -1
    }
    if (second_classifications[(aggregated_classifications != -1)][modID] == 1) {
      interim_classifications[(aggregated_classifications != -1)][modID] = 1
    }
  }
  aggregated_classifications = interim_classifications
  
  print("Third run!")
  run = 2
  
  # remove unimportant modules
  updated_modules3 = dat$modules[(aggregated_classifications != -1)]
  
  # initiate
  third_binomresults = rep(0, length(dat$modules))
  third_bernresults = zeros(num_iterations, length(dat$modules))
  third_vim = zeros(num_iterations, 2*length(dat$modules))
  iteration = 0
  
  # prepare
  third_binomresults[(aggregated_classifications == -1)] = NA
  third_bernresults[, (aggregated_classifications == -1)] = NA
  
  if (isempty(updated_modules3)) {
    local_borutares = list(
      networkID = networkID,
      causalmodules = dat$causal_modules,
      aggregated_classifications = aggregated_classifications,
      first_vim = first_vim,
      first_bernresults = first_bernresults,
      first_binomresults = first_binomresults,
      first_classifications = first_classifications,
      updated_modules2 = updated_modules2,
      second_vim = second_vim,
      second_bernresults = second_bernresults,
      second_binomresults = second_binomresults,
      second_classifications = second_classifications,
      updated_modules3 = updated_modules3,
      third_vim = third_vim,
      third_bernresults = third_bernresults,
      third_binomresults = third_binomresults,
      third_classifications = rep(-1, length(dat$modules)),
      updated_modules4 = updated_modules4,
      fourth_vim = zeros(num_iterations, 2*length(dat$modules)),
      fourth_bernresults = zeros(num_iterations, length(dat$modules)),
      fourth_binomresults = rep(0, length(dat$modules)),
      fourth_classifications = rep(-1, length(dat$modules))
    )
    
    # save results
    load(file = saveroot)
    borutares[[networkID]] = local_borutares
    save(borutares, file = saveroot)
    
    return(local_borutares)
  }
  
  # boruta loop
  while (iteration < num_iterations) {
    # advance iteration
    iteration = iteration + 1
    print(paste("Iteration:", iteration))
    
    # tic()
    
    # add shadow modules
    original_and_shadow_modules = updated_modules3
    for (module in updated_modules3) {
      original_and_shadow_modules = c(original_and_shadow_modules, list(module+1000))
    }
    
    # add shadow variables
    rfdat = dat$data
    rfdat = cbind(rfdat, rfdat[,-1])
    colnames(rfdat)[-1] = paste0("X", 1:(ncol(rfdat)-1))
    # print(colnames(rfdat))
    for (col in (((ncol(rfdat)-1)/2)+1):ncol(rfdat)) {
      rfdat[,col] = sample(rfdat[,col])
    }
    
    # run rf
    rf = simpleRFNetwork(pheno ~ .,
                         data = rfdat,
                         num_trees=num_trees,
                         num_threads=num_threads,
                         splitobject="module",
                         splitmethod=splitmethod,
                         varselection="none",
                         mtry="root",
                         varclusters = original_and_shadow_modules,
                         seed = as.integer(iteration + run*num_iterations)
    )
    
    # run var importance
    varimp = rf$variableImportance(type = importance, num_threads = num_threads)
    third_vim[iteration, (aggregated_classifications != -1)] = varimp
    
    # get max var imp of shadow modules
    shadow_varimp = varimp[(length(updated_modules3)+1):length(varimp)]
    max_shadow_varimp = max(shadow_varimp)
    
    # save results
    for (i in 1:length(updated_modules3) ) {
      if (varimp[i] > max_shadow_varimp) {
        third_binomresults[(aggregated_classifications != -1)][i] = third_binomresults[(aggregated_classifications != -1)][i] + 1
        third_bernresults[iteration, (aggregated_classifications != -1)][i] = 1
      }
    }
    
    # clean variables
    rm(rf)
    rm(rfdat)
    
    # toc()
  }
  
  # classify modules
  third_classifications = rep(0, length(dat$modules))
  third_classifications[(aggregated_classifications == -1)] = NA
  upcutoff = qbinom(0.975, num_iterations, 0.5)
  downcutoff = qbinom(0.025, num_iterations, 0.5)
  for (i in 1:length(updated_modules3) ) {
    if (third_binomresults[(aggregated_classifications != -1)][i] > upcutoff) {
      third_classifications[(aggregated_classifications != -1)][i] = 1
    }
    if (third_binomresults[(aggregated_classifications != -1)][i] < downcutoff) {
      third_classifications[(aggregated_classifications != -1)][i] = -1
    }
  }
  
  # aggregate classifications
  interim_classifications = aggregated_classifications
  for (modID in 1:length(updated_modules3)) {
    if (third_classifications[(aggregated_classifications != -1)][modID] == -1) {
      interim_classifications[(aggregated_classifications != -1)][modID] = -1
    }
    if (third_classifications[(aggregated_classifications != -1)][modID] == 1) {
      interim_classifications[(aggregated_classifications != -1)][modID] = 1
    }
  }
  aggregated_classifications = interim_classifications
  
  # remove unimportant modules
  updated_modules4 = dat$modules[(aggregated_classifications != -1)]
  
  print("Fourth run!")
  run = 3
  
  # initiate
  fourth_binomresults = rep(0, length(dat$modules))
  fourth_bernresults = zeros(num_iterations, length(dat$modules))
  fourth_vim = zeros(num_iterations, 2*length(dat$modules))
  iteration = 0
  
  # prepare
  fourth_binomresults[(aggregated_classifications == -1)] = NA
  fourth_bernresults[, (aggregated_classifications == -1)] = NA
  
  if (isempty(updated_modules4)) {
    local_borutares = list(
      networkID = networkID,
      causalmodules = dat$causal_modules,
      aggregated_classifications = aggregated_classifications,
      first_vim = first_vim,
      first_bernresults = first_bernresults,
      first_binomresults = first_binomresults,
      first_classifications = first_classifications,
      updated_modules2 = updated_modules2,
      second_vim = second_vim,
      second_bernresults = second_bernresults,
      second_binomresults = second_binomresults,
      second_classifications = second_classifications,
      updated_modules3 = updated_modules3,
      third_vim = third_vim,
      third_bernresults = third_bernresults,
      third_binomresults = third_binomresults,
      third_classifications = third_classifications,
      updated_modules4 = updated_modules4,
      fourth_vim = fourth_vim,
      fourth_bernresults = fourth_bernresults,
      fourth_binomresults = fourth_binomresults,
      fourth_classifications = rep(-1, length(dat$modules))
    )
    
    # save results
    load(file = saveroot)
    borutares[[networkID]] = local_borutares
    save(borutares, file = saveroot)
    
    return(local_borutares)
  }
  
  # boruta loop
  while (iteration < num_iterations) {
    # advance iteration
    iteration = iteration + 1
    print(paste("Iteration:", iteration))
    
    # tic()
    
    # add shadow modules
    original_and_shadow_modules = updated_modules4
    for (module in updated_modules4) {
      original_and_shadow_modules = c(original_and_shadow_modules, list(module+1000))
    }
    
    # add shadow variables
    rfdat = dat$data
    rfdat = cbind(rfdat, rfdat[,-1])
    colnames(rfdat)[-1] = paste0("X", 1:(ncol(rfdat)-1))
    # print(colnames(rfdat))
    for (col in (((ncol(rfdat)-1)/2)+1):ncol(rfdat)) {
      rfdat[,col] = sample(rfdat[,col])
    }
    
    # run rf
    rf = simpleRFNetwork(pheno ~ .,
                         data = rfdat,
                         num_trees=num_trees,
                         num_threads=num_threads,
                         splitobject="module",
                         splitmethod=splitmethod,
                         varselection="none",
                         mtry="root",
                         varclusters = original_and_shadow_modules,
                         seed = as.integer(iteration + run*num_iterations)
    )
    
    # run var importance
    varimp = rf$variableImportance(type = importance, num_threads = num_threads)
    fourth_vim[iteration, (aggregated_classifications != -1)] = varimp
    
    # get max var imp of shadow modules
    shadow_varimp = varimp[(length(updated_modules4)+1):length(varimp)]
    max_shadow_varimp = max(shadow_varimp)
    
    # save results
    for (i in 1:length(updated_modules4) ) {
      if (varimp[i] > max_shadow_varimp) {
        fourth_binomresults[(aggregated_classifications != -1)][i] = fourth_binomresults[(aggregated_classifications != -1)][i] + 1
        fourth_bernresults[iteration, (aggregated_classifications != -1)][i] = 1
      }
    }
    
    # clean variables
    rm(rf)
    rm(rfdat)
    
    # toc()
  }
  
  # classify modules
  fourth_classifications = rep(0, length(dat$modules))
  fourth_classifications[(aggregated_classifications == -1)] = NA
  upcutoff = qbinom(0.95, num_iterations, 0.5)
  downcutoff = qbinom(0.025, num_iterations, 0.5)
  for (i in 1:length(updated_modules4) ) {
    if (fourth_binomresults[(aggregated_classifications != -1)][i] > upcutoff) {
      fourth_classifications[(aggregated_classifications != -1)][i] = 1
    }
    if (fourth_binomresults[(aggregated_classifications != -1)][i] < downcutoff) {
      fourth_classifications[(aggregated_classifications != -1)][i] = -1
    }
  }
  
  # aggregate classifications
  interim_classifications = aggregated_classifications
  for (modID in 1:length(updated_modules4)) {
    if (fourth_classifications[(aggregated_classifications != -1)][modID] == -1) {
      interim_classifications[(aggregated_classifications != -1)][modID] = -1
    }
    if (fourth_classifications[(aggregated_classifications != -1)][modID] == 1) {
      interim_classifications[(aggregated_classifications != -1)][modID] = 1
    }
  }
  aggregated_classifications = interim_classifications
  
  local_borutares = list(
    networkID = networkID,
    causalmodules = dat$causal_modules,
    aggregated_classifications = aggregated_classifications,
    first_vim = first_vim,
    first_bernresults = first_bernresults,
    first_binomresults = first_binomresults,
    first_classifications = first_classifications,
    updated_modules2 = updated_modules2,
    second_vim = second_vim,
    second_bernresults = second_bernresults,
    second_binomresults = second_binomresults,
    second_classifications = second_classifications,
    updated_modules3 = updated_modules3,
    third_vim = third_vim,
    third_bernresults = third_bernresults,
    third_binomresults = third_binomresults,
    third_classifications = third_classifications,
    updated_modules4 = updated_modules4,
    fourth_vim = fourth_vim,
    fourth_bernresults = fourth_bernresults,
    fourth_binomresults = fourth_binomresults,
    fourth_classifications = fourth_classifications
  )
  
  # save results
  load(file = saveroot)
  borutares[[networkID]] = local_borutares
  save(borutares, file = saveroot)
  
  return(local_borutares)
}
