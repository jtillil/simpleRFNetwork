## Initialize batch tree growth
## @splitNode_batch
grow_batch <- function(treeIDs) {

  ## Initialize bootstraps
  trees <<- mclapply(trees, function(x) {
    x$init_grow(replace)
    x
  }, mc.cores = num_threads)

  ## Initialize splitNode_batch
  splitNode_batch(treeIDs, 1)
}

## Handle batch node splitting
## @splitNode_batch
## @splitNodeInternal_batch
## @tree[[id]]$estimate
splitNode_batch <- function(treeIDs, nodeID) {

  ## Possible split clusters, maximum mtry
  if (length(varclusters) <= mtry) {
    possible_split_clusterIDs <- 1:length(varclusters)
  } else {
    possible_split_clusterIDs <- sample(
      x = 1:length(varclusters),
      size = mtry,
      replace = FALSE)
  }
  
  ## Split node
  split <- splitNodeInternal_batch(nodeID, possible_split_clusterIDs)

  ## Calculate node depth, size and gini impurity
  depth <- 1
  current_nodeID <- nodeID
  while(current_nodeID != 1) {
    current_nodeID <- go_to_parent(child_nodeIDs, current_nodeID)
    depth <- depth + 1
  }
  depths[nodeID] <<- as.integer(depth)
  sizes[nodeID] <<- length(sampleIDs[[nodeID]])
  impurities[nodeID] <<- 1-(
    (sum(data$subset(sampleIDs[[nodeID]], 1) == "1")/sizes[nodeID])^2+
    (sum(data$subset(sampleIDs[[nodeID]], 1) == "0")/sizes[nodeID])^2
  )
  
  if (!is.null(split$clusterID)) {
    ## Read timing metrics
    linearcomb_times[[nodeID]] <<- split$linearcomb_times
    node_times[nodeID] <<- split$node_time

    # ## Assign split
    split_clusterIDs[[nodeID]] <<- split$clusterID
    # split_selectedVarIDs[[nodeID]] <<- split$selectedVarIDs
    split_coefficients[[nodeID]] <<- split$coefficients
    split_values[[nodeID]] <<- split$value
    
    ## Create child nodes
    left_child <- length(sampleIDs) + 1
    right_child <- length(sampleIDs) + 2
    child_nodeIDs[[nodeID]] <<- c(left_child, right_child)
    
    ## For each sample in node, assign to left or right child
    if (varselection == "none") {
      idx <- as.matrix(data$subset(sampleIDs[[nodeID]], varclusters[[split$clusterID]] + 1)) %*% split$coefficients <= split$value
    } else {
      idx <- as.matrix(data$subset(sampleIDs[[nodeID]], split_selectedVarIDs[[nodeID]] + 1)) %*% split$coefficients <= split$value
    }
    
    sampleIDs[[left_child]] <<- sampleIDs[[nodeID]][idx]
    sampleIDs[[right_child]] <<- sampleIDs[[nodeID]][!idx]
    
    ## Recursively call split node batch on child nodes
    splitNode_batch(treeIDs, left_child)
    splitNode_batch(treeIDs, right_child)
  } else {
    ## Read timing metrics
    linearcomb_times[[nodeID]] <<- NA
    node_times[nodeID] <<- NA

    ## Compute estimate for terminal node
    split_clusterIDs[[nodeID]] <<- NA
    # split_selectedVarIDs[[nodeID]] <<- NA
    split_coefficients[[nodeID]] <<- NA
    split_values[[nodeID]] <<- estimate(nodeID)
  }

}

## Initiate finding of best split in batch for one node
## @findBestSplit
splitNodeInternal_batch = function(nodeID, possible_split_clusterIDs) {
  ## Check node size, stop if minimum reached
  if (length(sampleIDs[[nodeID]]) <= min_node_size) {
    return(NULL)
  }

  ## Get response
  response <- data$subset(sampleIDs[[nodeID]], 1)
  
  ## Stop if node is pure
  if (length(unique(response)) == 1) {
    return(NULL)
  }

  ## IF LDA: stop if node has only one observation left for one class
  if (splitmethod == "LDA" & (length(response[response == 1]) <= 1 | length(response[response == 0]) <= 1)) {
    return(NULL)
  }
  
  ## Find best split
  return(findBestSplit(nodeID, possible_split_clusterIDs, response))
}

## Find best split in batch
findBestSplit_batch = function(nodeID, possible_split_clusterIDs, response) {
  
  ## Initialize
  best_split <- NULL
  best_split$clusterID <- -1
  # best_split$selectedVarIDs <- -1
  best_split$coefficients <- -1
  best_split$value <- -1
  best_split$decrease <- -1
  best_split$linearcomb_time <- -1
  
  ## Initialize array for individual linear combination time measurement
  linearcomb_times <- c()
  
  ## Start timing for node splitting time measurement
  tic()
  
  ## For all possible variable clusters
  for (split_clusterID in possible_split_clusterIDs) {
    
    ## Read data values from samples in current node
    data_values <- data$subset(sampleIDs[[nodeID]], varclusters[[split_clusterID]] + 1)

    ## IF LDA: calculate covariance matrix and skip if singular
    if (splitmethod == "LDA") {
      ## Calculate mean of both covariance matrices due to homoscedasticity
      mat <- 0.5*(cova(as.matrix(data_values[response == 0,]), center=TRUE, large=FALSE) +
                  cova(as.matrix(data_values[response == 1,]), center=TRUE, large=FALSE))
      ## Condition matrix by adding 1e-10 to diagonal elements that are 0
      sapply(1:ncol(mat), function(j) {
        if (mat[j,j] == 0) {
          mat[j,j] <<- 1e-10
        }
      })
      ## Check if singular
      if (!is.positive.definite(mat)) {
        next
      }
    } else {
      mat <- NULL
    }
    
    ## IF CART: read IQR scaled data values from samples in current node
    if (splitmethod == "CART" | splitmethod == "CART_fast") {
      IQR_data_values <- Data$new(data = IQR_data$subset(sampleIDs[[nodeID]], varclusters[[split_clusterID]]))
    } else {
      IQR_data_values <- NULL
    }
    
    ## Find best split
    best_split = findBestSplitCoefs_batch(split_clusterID, best_split, data_values, IQR_data_values, response, mat)
    
    ## Save time measurement for single linear combination
    linearcomb_times <- c(linearcomb_times, best_split$linearcomb_time)
  }
  
  ## Stop timing for node splitting time measurement
  node_time <- toc(quiet = TRUE)
  
  if (best_split$clusterID < 0) {
    ## Stop if no good split found
    result <- NULL
    result$clusterID <- NULL
    result$linearcomb_times <- linearcomb_times
    result$node_time <- as.numeric(node_time$toc - node_time$tic)
    return(result)
  } else {
    ## Return best split
    result <- NULL
    result$clusterID <- as.integer(best_split$clusterID)
    # result$selectedVarIDs <- NULL
    result$coefficients <- best_split$coefficients
    result$value <- best_split$value
    result$linearcomb_times <- linearcomb_times
    result$node_time <- as.numeric(node_time$toc - node_time$tic)
    return(result)
  }
}

## Find coefficients for linear combination of variables
findBestSplitCoefs_batch = function(treeIDs, split_clusterID, best_split, data_values, response, mat=NULL) {
  
    ## Start timing for individual linear combination time measurement
    tic()

    ## Init lists
    inputs <- NULL
    denses <- NULL
    outputs <- NULL
    # models <- NULL

    ## Fill lists
    for (id in 1:length(treeIDs)) {
        inputs[[id]] <- layer_input(shape = dim(data_values)[[2]])
        denses[[id]] <- layer_dense(
            units = 1,
            activation = "sigmoid",
            kernel_initializer='glorot_uniform',
            bias_initializer='glorot_uniform')
        outputs[[id]] <- inputs[[id]] %>%
            denses[[id]]
        # models[[id]] <- keras_model(
        #     inputs = inputs[[id]],
        #     outputs = outputs[[id]]
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
            list(as.matrix(dat), as.matrix(dat)),
            list(as.numeric(resp)-1, as.numeric(resp)-1),
            epochs = 100,
            verbose = 1,
            callbacks = c(callback_early_stopping(
                monitor = "loss",
                min_delta = 0.01,
                mode = "min",
                patience = 10
            ))
        )

    ## Read coefficients and value
    coefList <- NULL
    valueList <- NULL
    for (id in 1:length(treeIDs)) {
        coefList[[id]] <- as.numeric(model$layers[[length(treeIDs) + id]]$weights[[1]]$numpy())
        valueList[[id]] <- -as.numeric(model$layers[[length(treeIDs) + id]]$weights[[2]]$numpy())

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
    
    ## Stop timing for individual linear combination time measurement
    linearcomb_time <- toc(quiet = TRUE)
    
    ## Count classes in childs
    idx <- as.matrix(data_values)%*%coefficients <= value
    class_counts_left <- tabulate(response[idx])
    class_counts_right <- tabulate(response[!idx])
    
    ## Skip cluster if one child empty
    if (sum(class_counts_left) == 0 | sum(class_counts_right) == 0) {
        best_split$linearcomb_time <- as.numeric(linearcomb_time$toc - linearcomb_time$tic)
        return(best_split)
    }
    
    if (splitrule == "Gini") {
        ## Decrease of impurity
        decrease <- sum(class_counts_left^2)/sum(class_counts_left) + 
        sum(class_counts_right^2)/sum(class_counts_right)
    } else {
        stop("Unknown splitrule.")
    }
    
    ## Use this cluster for the split if decrease better than from earlier clusters
    if (decrease > best_split$decrease) {
        best_split$clusterID <- split_clusterID
        best_split$coefficients <- coefficients
        best_split$value <- value
        best_split$decrease <- decrease
        best_split$linearcomb_time <- as.numeric(linearcomb_time$toc - linearcomb_time$tic)
    } else {
        best_split$linearcomb_time <- as.numeric(linearcomb_time$toc - linearcomb_time$tic)
    }
    
    return(best_split)
}