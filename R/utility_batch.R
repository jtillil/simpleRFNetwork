## Find coefficients for linear combination of variables
findBestSplitCoefs_batch = function(split_clusterIDList, best_splitList, data_valList, responseList) {
  
  ## Start timing for combined linear combination time measurement
  tic()

  ## Init lists
  inputs <- NULL
  outputs <- NULL

  ## Fill lists
  for (id in 1:length(split_clusterIDList)) {
    inputs[[id]] <- layer_input(shape = dim(data_valList[[id]])[[2]])
    outputs[[id]] <- inputs[[id]] %>%
      layer_dense(
        units = 1,
        activation = "sigmoid",
        kernel_initializer='glorot_uniform',
        bias_initializer='glorot_uniform')
  }

  ## Create batch model
  model <- keras_model(
    inputs = inputs,
    outputs = outputs
  )

  ## Compile batch model
  model %>%
    compile(
      optimizer = tf$keras$optimizers$Adam(learning_rate = 0.1),
      loss = "binary_crossentropy",
      run_eagerly = TRUE
    )
  
  ## Train batch model
  model %>%
    fit(
      data_valList,
      responseList,
      epochs = 100,
      verbose = 0,
      callbacks = c(callback_early_stopping(
        monitor = "loss",
        min_delta = 0.01,
        mode = "min",
        patience = 10
      ))
    )

  ## Read coefficients and values
  coefList <- NULL
  valueList <- NULL
  for (id in 1:length(split_clusterIDList)) {
    coefList[[id]] <- as.numeric(model$layers[[length(split_clusterIDList) + id]]$weights[[1]]$numpy())
    valueList[[id]] <- -as.numeric(model$layers[[length(split_clusterIDList) + id]]$weights[[2]]$numpy())

    ## Restrict coefficients to norm 1 for a unique solution up to factor -1
    coef_norm <- Norm(as.matrix(coefList[[id]]), "F")
    valueList[[id]] <- valueList[[id]]/coef_norm
    coefList[[id]] <- coefList[[id]]/coef_norm

    ## Restrict coefficients to positive split value for a completely unique solution
    if (valueList[[id]] < 0) {
      valueList[[id]] <- -valueList[[id]]
      coefList[[id]] <- -coefList[[id]]
    }
  }
  
  ## Stop timing for combined linear combination time measurement
  linearcomb_time <- toc(quiet = TRUE)
  
  skipvec <- rep(FALSE, length(split_clusterIDList))
  ## Process all cluster results
  for (id in 1:length(split_clusterIDList)) {
    ## Count classes in childs
    idx <- as.matrix(data_valList[[id]])%*%coefList[[id]] <= valueList[[id]]
    class_counts_left <- tabulate(responseList[[id]][idx])
    class_counts_right <- tabulate(responseList[[id]][!idx])
  
    ## Skip cluster if one child empty
    if (sum(class_counts_left) == 0 | sum(class_counts_right) == 0) {
      skipvec[id] <- TRUE
    }
    
    if (!skipvec[id]) {
      ## Decrease of impurity
      decrease <- sum(class_counts_left^2)/sum(class_counts_left) + 
        sum(class_counts_right^2)/sum(class_counts_right)
      
      ## Use this cluster for the split if decrease better than from earlier clusters
      best_splitList[[id]]$clusterID <- split_clusterIDList[id]
      best_splitList[[id]]$coefficients <- coefList[[id]]
      best_splitList[[id]]$value <- valueList[[id]]
      best_splitList[[id]]$decrease <- decrease
      best_splitList[[id]]$linearcomb_time <- as.numeric(linearcomb_time$toc - linearcomb_time$tic)
    } 
  }
  
  if (sum(skipvec) < length(split_clusterIDList)) {
    ## Sort by decrease and extract best split
    best_id <- 0
    best_decrease <- 0
    for (id in (1:length(split_clusterIDList))[!skipvec]) {
      if (best_splitList[[id]]$decrease > best_decrease) {
        best_id <- id
        best_decrease <- best_splitList[[id]]$decrease
      }
    }

    ## Return best split
    return(best_splitList[[best_id]])

  } else {
    ## Create empty split
    empty_split <- NULL
    empty_split$clusterID <- -1
    empty_split$linearcomb_time <- as.numeric(linearcomb_time$toc - linearcomb_time$tic)
    
    ## Return empty split
    return(empty_split)
  }
}