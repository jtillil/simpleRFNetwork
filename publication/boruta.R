boruta <- function(dat, networkID, splitmethod, importance, num_trees, num_threads, num_iterations, seed, saveroot) {
  # seed
  set.seed(seed)
  
  # initiate
  first_binomresults = rep(0, length(dat$modules))
  first_bernresults = zeros(num_iterations, length(dat$modules))
  first_vim = zeros(num_iterations, 2*length(dat$modules))
  iteration = 0
  
  print(paste("This is", splitmethod, "Boruta"))
  print(saveroot)
  print("First run!")
  
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
      varclusters = original_and_shadow_modules,
      seed = as.integer(iteration)
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
  }
  
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
  
  # remove unimportant modules
  updated_modules = dat$modules[(first_classifications != -1)]
  
  print("Second run!")
  
  # initiate
  second_binomresults = rep(0, length(dat$modules))
  second_bernresults = zeros(num_iterations, length(dat$modules))
  second_vim = zeros(num_iterations, 2*length(dat$modules))
  iteration = 0
  
  # prepare
  second_binomresults[(first_classifications == -1)] = NA
  second_bernresults[, (first_classifications == -1)] = NA
  
  # boruta loop
  while (iteration < num_iterations) {
    # advance iteration
    iteration = iteration + 1
    print(paste("Iteration:", iteration))
    
    # tic()
    
    # add shadow modules
    original_and_shadow_modules = updated_modules
    for (module in updated_modules) {
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
                         varclusters = original_and_shadow_modules,
                         seed = as.integer(iteration)
    )
    
    # run var importance
    varimp = rf$variableImportance(type = importance, num_threads = num_threads)
    second_vim[iteration, (first_classifications != -1)] = varimp
    
    # get max var imp of shadow modules
    shadow_varimp = varimp[(length(updated_modules)+1):length(varimp)]
    max_shadow_varimp = max(shadow_varimp)
    
    # save results
    for (i in 1:length(updated_modules) ) {
      if (varimp[i] > max_shadow_varimp) {
        second_binomresults[(first_classifications != -1)][i] = second_binomresults[(first_classifications != -1)][i] + 1
        second_bernresults[iteration, (first_classifications != -1)][i] = 1
      }
    }
    
    # clean variables
    rm(rf)
    rm(rfdat)
    
    # toc()
  }
  
  # classify modules
  second_classifications = rep(0, length(dat$modules))
  second_classifications[(first_classifications != -1)] = NA
  upcutoff = qbinom(0.95, num_iterations, 0.5)
  
  for (i in 1:length(updated_modules) ) {
    if (second_binomresults[(first_classifications != -1)][i] > upcutoff) {
      second_classifications[(first_classifications != -1)][i] = 1
    }
  }
  
  # save results
  load(file = saveroot)
  borutares[[networkID]] = list(
    networkID = networkID,
    first_vim = first_vim,
    first_bernresults = first_bernresults,
    first_binomresults = first_binomresults,
    first_classification = first_classifications,
    updated_modules = updated_modules,
    second_vim = second_vim,
    second_bernresults = second_bernresults,
    second_binomresults = second_binomresults,
    second_classification = second_classifications,
    causalmodules = dat$causal_modules
  )
  save(borutares, file = saveroot)
  
  return(list(
    networkID = networkID,
    first_vim = first_vim,
    first_bernresults = first_bernresults,
    first_binomresults = first_binomresults,
    first_classification = first_classifications,
    updated_modules = updated_modules,
    second_vim = second_vim,
    second_bernresults = second_bernresults,
    second_binomresults = second_binomresults,
    second_classification = second_classifications,
    causalmodules = dat$causal_modules
  ))
}
