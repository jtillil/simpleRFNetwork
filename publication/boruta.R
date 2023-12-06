boruta <- function(dat, networkID, splitmethod, importance, num_trees, num_threads, num_iterations, seed, saveroot) {
  # seed
  set.seed(seed)
  
  # initiate
  binomresults = rep(0, length(dat$modules))
  bernresults = zeros(num_iterations, length(dat$modules))
  iteration = 0
  
  print(paste("This is", splitmethod, "Boruta"))
  print(saveroot)
  
  # boruta loop
  while (iteration < num_iterations) {
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
    rfdat = dat$data[1:500,]
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
      varclusters = original_and_shadow_modules
    )
    
    # run var importance
    varimp = rf$variableImportance(type = importance, num_threads = num_threads)
    
    # get max var imp of shadow modules
    shadow_varimp = varimp[(length(dat$modules)+1):length(varimp)]
    max_shadow_varimp = max(shadow_varimp)
    
    # classify modules
    for (i in 1:length(dat$modules) ) {
      if (varimp[i] > max_shadow_varimp) {
        binomresults[i] = binomresults[i] + 1
        bernresults[iteration, i] = 1
      }
    }
    
    # clean variables
    rm(rf)
    rm(rfdat)
    
    # toc()
  }
  
  # classify modules
  module_classifications = rep(0, length(dat$modules))
  cutoff = qbinom(0.95, num_iterations, 0.5)
  for (i in 1:length(dat$modules) ) {
    if (binomresults[i] > cutoff) {
      module_classifications[i] = 1
    }
  }
  
  # save results
  load(file = saveroot)
  borutares[[networkID]] = list(
    networkID = networkID,
    bernresults = bernresults,
    binomresults = binomresults,
    classification = module_classifications,
    causalmodules = dat$causal_modules
  )
  save(borutares, file = saveroot)
  
  return(list(
    bernresults = bernresults,
    binomresults = binomresults,
    classification = module_classifications,
    causalmodules = dat$causal_modules
  ))
}
