setwd(getSrcDirectory(function(){})[1])
source("./source_files.R")

library(gglasso)
library(MLGL)

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
# scenarios = rbind(scenarios, c(100, 3000, 1000, 0, F, 0, 1))

# go through scenarios
for (i in 1:nrow(scenarios)) {
  print(i)
  
  # read scenario
  scenario = scenarios[i,]
  
  # load data
  datroot = paste0(
    "./data/ndclassif",
    # "_datasc", data_scenario,
    "_nn", scenario$n_networks,
    "_ng", scenario$n_genes,
    "_ns", scenario$n_samples,
    "_ndm", scenario$n_disease_modules,
    "_mdg", scenario$main_disease_gene,
    "_pdg", scenario$prop_disease_genes,
    "_ab", scenario$average_beta,
    ".Rdata"
  )
  print(getwd())
  print(datroot)
  load(datroot)
  print(summary(as.factor(dat[[1]]$data[,1])))
  
  # save root
  saveroot = paste0(
    "./results/grplasso",
    "_nn", scenario$n_networks,
    "_ng", scenario$n_genes,
    "_ns", scenario$n_samples,
    "_ndm", scenario$n_disease_modules,
    "_mdg", scenario$main_disease_gene,
    "_pdg", scenario$prop_disease_genes,
    "_ab", scenario$average_beta,
    ".Rdata"
  )
  grplassores = list()
  # save(grplassores, file = saveroot)
  
  # run grp lasso
  grplassores = lapply(
    1:length(dat),
    function(i) {
      tic()
      
      # read X and y
      X = as.matrix(dat[[i]]$data[, -1])
      y = 2*as.numeric(dat[[i]]$data[, 1]) - 3
      
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
      
      # optimize lambda
      print(paste("GRP lasso for Network Nr", i, "optimize lambda."))
      gr_cv <- gglasso::cv.gglasso(x=Xb, y=y, group=groupb, 
                          pred.loss="L2", 
                          # intercept = F, 
                          nfolds=5)
      # x11(); plot(gr_cv)
      # paste(gr_cv$lambda.min, gr_cv$lambda.1se)
      
      # calc weight
      weight <- as.numeric(sqrt(table(groupb)))
      
      if (i %in% c(14, 15)) {
        # browser()
      }
      
      # perform gglasso
      # gr = gglasso(Xb, y, groupb, pf = weight, lambda = gr_cv$lambda.1se+0.1, intercept = F, loss = "logit"
      print(paste("GRP lasso for Network Nr", i, "started."))
      gr = gglasso::gglasso(x = Xb,
                            y = y,
                            # var = var_gglasso,
                            group = groupb,
                            lambda = gr_cv$lambda.min,
                            # intercept = F,
                            loss="logit")
      # gr = MLGL::overlapgglasso(X = X,
      #       y = y,
      #       var = var_gglasso,
      #       group = group_gglasso,
      #       lambda = gr_cv$lambda.min,
      #       loss="logit",
      #       intercept = F)
      
      # print(gr$non0)
      
      print(paste("GRP lasso for Network Nr", i, "finished!"))
      toc()
      
      # selected = unique(unlist(gr$group))
      selected = unique(groupb[gr$beta != 0])
      # selected = unique(gr$group[[1]])
      
      return(list(
        # selecteds10 = unique(gr$group$s10),
        # selecteds30 = unique(gr$group$s30),
        # selecteds20 = unique(gr$group$s20),
        selected = selected,
        causal = dat[[i]]$causal_modules,
        number_correctly_selected = sum(dat[[i]]$causal_modules %in% selected)
      ))
    }
  )
  
  # save grp lasso
  save(grplassores, file = saveroot)
  
  # clear dat
  rm(dat)
}
