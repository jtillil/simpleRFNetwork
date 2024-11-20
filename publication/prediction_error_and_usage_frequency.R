setwd(getSrcDirectory(function(){})[1])
source("./source_files.R")

library(simpleRFNetwork)
library(SeqNet)
library(parallel)
library(pracma)
library(tictoc)
library(Rfast)
library(matrixcalc)

library(ridge)
library(glmnet)

library(ggplot2)

library(ranger)

#### setup

# set scenarios
n_networks = c(100)
n_genes = c(1000)
n_samples = c(1000)
n_disease_modules = c(1, 2)
main_disease_gene = c(F)
# main_disease_gene = c(T, F)
# average_beta = c(1)
average_beta = c(0.5, 1, 2)
prop_disease_genes = c(0.5)

scenarios = expand.grid(
  n_networks = n_networks,
  n_genes = n_genes,
  n_samples = n_samples,
  n_disease_modules = n_disease_modules,
  main_disease_gene = main_disease_gene,
  average_beta = average_beta,
  prop_disease_genes = prop_disease_genes
)
scenarios = rbind(scenarios, c(100, 1000, 1000, 0, F, 0, 0.5))
importance = "permutation"

#### prediction error RF

pred_err = function(network) {
  # print(ID_network)
  #
  # network = dat[[ID_network]]
  rfdat = network$data[1:500,]
  preddat = as.matrix(network$data[501:1000,-1])
  predlabel = network$data[501:1000,1]
  modules = network$modules

  rf = simpleRFNetwork(
    pheno ~ .,
    data = rfdat,
    num_trees=500,
    num_threads=1,
    splitobject="module",
    splitmethod=splitmethod,
    varselection="none",
    mtry="root",
    varclusters = modules,
    seed = 1L
  )

  print(length(rf$trees[[1]]$split_clusterIDs))

  res = list()
  res$method = splitmethod
  # err = rf$predictionErrorForestAndTrees()
  # err = rf$predictionErrorForest()
  pred = rf$predict(preddat)
  res$prederr = sum(pred != predlabel) / 500

  split_counts = numeric(length(modules))
  for (treeID in 1:length(rf$trees)) {
    clusterIDs = rf$trees[[treeID]]$split_clusterIDs
    for (cluster in clusterIDs) {
      if (!is.na(cluster)) {
        split_counts[cluster] = split_counts[cluster] + 1
      }
    }
  }
  res$split_group_counts = split_counts

  res$modules = modules
  res$causal_modules = network$causal_modules

  return(res)
}

# for (i in 1:1) {
for (i in 1:nrow(scenarios)) {
  # read scenario
  scenario = scenarios[i,]
  print(scenario)

  for (method in c("LDA", "PCA", "logridge1")) {
    print(method)
  # for (method in c("LDA")) {
    datroot = paste0(
      # "./resclassif",
      "./data/ndclassif",
      "_nn", 100,
      "_ng", 1000,
      "_ns", 1000,
      "_ndm", scenario$n_disease_modules,
      "_mdg", 0,
      "_pdg", 0.5,
      "_ab", scenario$average_beta,
      ".Rdata"
    )
    load(datroot)

    splitmethod = method

    prederr_res = mclapply(dat, pred_err, mc.cores = 60)

    saveroot = paste0(
      "./results/prederrres_",
      splitmethod,
      "_ndm", scenario$n_disease_modules,
      "_ab", scenario$average_beta,
      ".Rdata"
    )
    save(prederr_res, file = saveroot)
  }
}

#### calculate ranger prediction errors

for (i in 1:nrow(scenarios)) {
  # read scenario
  scenario = scenarios[i,]
  print(scenario)

  for (method in c("ranger")) {
    print(method)
    datroot = paste0(
      "./data/ndclassif",
      "_nn", 100,
      "_ng", 1000,
      "_ns", 1000,
      "_ndm", scenario$n_disease_modules,
      "_mdg", 0,
      "_pdg", 0.5,
      "_ab", scenario$average_beta,
      ".Rdata"
    )
    load(datroot)

    prederr_res = list()
    for (i in 1:100) {
      print(i)
      network = dat[[i]]
      rangerrf = ranger(
        dependent.variable.name = "pheno",
        data = network$data[1:500, ],
        num.trees = 500,
        seed = i
      )
      pred = predict(rangerrf, network$data[501:1000, -1])$predictions
      prederr = sum(pred != network$data[501:1000, 1]) / 500
      res = list()
      res$prederr = prederr
      prederr_res[[i]] = res
    }

    saveroot = paste0(
      "./results/prederrres_",
      method,
      "_ndm", scenario$n_disease_modules,
      "_ab", scenario$average_beta,
      ".Rdata"
    )
    save(prederr_res, file = saveroot)
  }
}

#### calculate group lasso prediction errors

library(gglasso)

for (i in 1:nrow(scenarios)) {
  # read scenario
  scenario = scenarios[i,]
  print(scenario)
  
  for (method in c("grplasso")) {
    print(method)
    datroot = paste0(
      "./data/ndclassif",
      "_nn", 100,
      "_ng", 1000,
      "_ns", 1000,
      "_ndm", scenario$n_disease_modules,
      "_mdg", 0,
      "_pdg", 0.5,
      "_ab", scenario$average_beta,
      ".Rdata"
    )
    load(datroot)
    
    prederr_res = list()
    for (i in 1:100) {
      print(i)
      
      tic()
      
      # read X and y
      X = as.matrix(dat[[i]]$data[1:500, -1])
      y = 2*as.numeric(dat[[i]]$data[1:500, 1]) - 3
      
      preddat = as.matrix(dat[[i]]$data[501:1000, -1])
      predlabels = 2*as.numeric(dat[[i]]$data[501:1000, 1]) - 3
      
      # create var and group vectors
      var_gglasso = c()
      group_gglasso = c()
      for (nmodule in 1:length(dat[[i]]$modules)) {
        for (variable in dat[[i]]$modules[[nmodule]]) {
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
      print(paste("GRP lasso for Network Nr", i, "optimize lambda."))
      gr_cv <- gglasso::cv.gglasso(x=Xb, y=y, group=groupb, 
                                   pred.loss="loss", 
                                   loss = "logit",
                                   # intercept = F, 
                                   nfolds=5)
      
      # predict data
      pred = predict(gr_cv$gglasso.fit, preddatb)
      
      # print(gr$non0)
      
      print(paste("GRP lasso for Network Nr", i, "finished!"))
      toc()
      
      # pred = predict(rangerrf, network$data[501:1000, -1])$predictionss
      prederr = sum(pred != predlabels) / 500
      res = list()
      res$prederr = prederr
      prederr_res[[i]] = res
    }
    
    saveroot = paste0(
      "./results/prederrres_",
      method,
      "_ndm", scenario$n_disease_modules,
      "_ab", scenario$average_beta,
      ".Rdata"
    )
    save(prederr_res, file = saveroot)
  }
}

#### plot prediction errors

predictiondat <- data.frame(
  Method = rep(rep(c("grplasso", "ranger", "LDA", "logridge1", "PCA"), each = 100), 6),
  ID = rep(rep(1:100, times = 5), 6),
  Prederr = rep(0, 3000),
  ndm = c(rep(1, 1500), rep(2, 1500)),
  ndm_plot = c(rep("1 disease module per network", 1500), rep("2 disease modules per network", 1500)),
  ab = rep(c(rep(0.5, 500), rep(1, 500), rep(2, 500)), times = 2),
  ab_plot = rep(c(rep("beta = 0.5", 500), rep("beta = 1", 500), rep("beta = 2", 500)), times = 2)
)

for (Method in c("grplasso", "ranger", "LDA", "logridge1", "PCA")) {
  for (ndm in c(1, 2)) {
    for (ab in c(0.5, 1, 2)) {
      if (file.exists(paste0(
        "./serverresults/serverres_24_11_18/prederrres_",
        Method,
        "_ndm", ndm,
        "_ab", ab,
        ".Rdata"
      ))) {
        # load prederr_res
        load(paste0(
          "./serverresults/serverres_24_11_18/prederrres_",
          Method,
          "_ndm", ndm,
          "_ab", ab,
          ".Rdata"
        ))
        print(Method)
        
        for (ID in 1:100) {
          if (!is.na(prederr_res[[ID]][1])) {
            predictiondat[
              predictiondat$Method == Method &
                predictiondat$ab == ab &
                predictiondat$ndm == ndm &
                predictiondat$ID == ID, "Prederr"] = prederr_res[[ID]]$prederr
          }
        }
        
        # print(Method)
        # print(ndm)
        # print(ab)
        # print(median(predictiondat[
        #     predictiondat$Method == Method &
        #     predictiondat$ab == ab &
        #     predictiondat$ndm == ndm, "Prederr"]))
      }
    }
  }
}

predictiondat$Prederr[predictiondat$Method == "grplasso"] = predictiondat$Prederr[predictiondat$Method == "grplasso"]/100

predictiondat$Method[predictiondat$Method == "grplasso"] = "Group Lasso"
predictiondat$Method[predictiondat$Method == "ranger"] = "RF"
predictiondat$Method[predictiondat$Method == "logridge1"] = "Group RF: Ridge"
predictiondat$Method[predictiondat$Method == "LDA"] = "Group RF: LDA"
predictiondat$Method[predictiondat$Method == "PCA"] = "Group RF: PCA"
predictiondat$Method = factor(predictiondat$Method, levels = c("Group Lasso", "RF", "Group RF: LDA", "Group RF: Ridge", "Group RF: PCA"), ordered = TRUE)

library(gridExtra)

ggplot(predictiondat, aes(x = Method, y = Prederr, fill = Method)) +
  # geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  geom_boxplot() +
  labs(
    # title = "Number of detected disease modules per network",
    x = "Method",
    y = "Prediction error") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  # scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  scale_fill_discrete(name = "Splitmethod", labels = c("LDA", "Ridge", "PCA"), guide = "none") +
  # legend(c("LDA", "Ridge", "PCA")) +
  facet_grid(ndm_plot ~ ab_plot, scales = "free")

ggsave("figures/boxplot_Prediction_Error.pdf", width = 7, height = 5)

#### plot usage frequency by module size

Usagedat <- data.frame(
  Method = character(),
  Size = double(),
  Usagepercent = double()
)

for (method in c("LDA", "logridge1", "PCA")) {
  print(method)
  load(paste0(
    "results/prederrres_",
    method,
    "_ndm", 0,
    "_ab", 0,
    ".Rdata"
  ))
  datroot = paste0(
    "./data/ndclassif",
    "_nn", 100,
    "_ng", 1000,
    "_ns", 1000,
    "_ndm", 0,
    "_mdg", 0,
    "_pdg", 0.5,
    "_ab", 0,
    ".Rdata"
  )
  load(datroot)
  
  for (i in 1:length(borutares)) {
    print(i)
    res = prederr_res[[i]]
    network = dat[[i]]
    for (j in 1:length(network$modules)) {
      # Usagedat[nrow(Usagedat) + 1,] = c(method, as.numeric(lengths(network$modules)[j]), res$split_group_counts[j] * length(network$modules) / 500)
      Usagedat[nrow(Usagedat) + 1,] = c(method, as.numeric(lengths(network$modules)[j]), res$split_group_counts[j] / 500)
    }
  }
}
Usagedat$Method[Usagedat$Method == "logridge1"] = "Ridge"
Usagedat$Method = factor(Usagedat$Method, levels = c("LDA", "Ridge", "PCA"))
Usagedat$Size = as.numeric(Usagedat$Size)
Usagedat$Usagepercent = as.numeric(Usagedat$Usagepercent)

p = ggplot(Usagedat, aes(Size, Usagepercent)) +
  geom_point() +
  scale_x_continuous(breaks = 25*(0:4)) +
  theme_bw() +
  xlab("Module size") +
  ylab("Average number of times a module is used per tree") +
  facet_grid(rows = vars(Method), scales = "free")
plot(p)

ggsave("scatter_Usage_Null.pdf", width = 7, height = 5)
