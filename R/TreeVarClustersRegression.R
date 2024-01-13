
##' @title Regression tree class
##' @description Subclass for regression tree.
##' Contains all fields and methods used special for regression trees.
TreeVarClustersRegression <- setRefClass("TreeVarClustersRegression", 
  contains = "TreeVarClusters",
  fields = list(),
  methods = list(
 
    ## Initiate finding of best split for one node
    ## @findBestSplit
    splitNodeInternal = function(nodeID, possible_split_clusterIDs) {
      ## Check node size, stop if minimum reached
      if (length(sampleIDs[[nodeID]]) <= min_node_size) {
        return(NULL)
      }
      
      ## Get response
      response <- data$subset(sampleIDs[[nodeID]], 1)
      
      ## Stop if node is pure
      if (var(response) == 0) {
        return(NULL)
      }
      
      ## IF LDA: stop if node has only one observation left for one class
      # if (splitmethod == "LDA" & (length(response[response == 1]) <= 1 | length(response[response == 0]) <= 1)) {
      #   return(NULL)
      # }
      
      ## Find best split
      return(findBestSplit(nodeID, possible_split_clusterIDs, response))
    },
    
    ## Find best split
    ## @findBestSplitCoefs
    findBestSplit = function(nodeID, possible_split_clusterIDs, response) {
      
      ## Initialize
      best_split <- NULL
      best_split$clusterID <- -1
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
        
        ## Find best split
        best_split = findBestSplitCoefs(split_clusterID, best_split, data_values, response, mat)
        
        ## Save time measurement for single linear combination
        linearcomb_times <- c(linearcomb_times, best_split$linearcomb_time)
      }
      
      if (best_split$varID < 0) {
        ## Stop if no good split found
        return(NULL)
      } else {
        ## Return best split
        result <- NULL
        result$varID <- as.integer(best_split$varID)
        result$value <- best_split$value
        return(result)
      }      
    },
    
    ## Find coefficients for linear combination of variables
    findBestSplitCoefs = function(split_clusterID, best_split, data_values, IQR_data_values, response, mat=NULL) {
      
      ## Coerce all but the most frequent factor level to a single one
      ## Irrelevant, if exactly two factors
      ## Skip cluster if less than two levels in response
      # if (nlevels(droplevels(response)) < 2) {
      #   return(best_split)
      # } else {
      #   bin_response <- response
      #   most_freq_idx <- response==names(which.max.random(summary(response)))
      #   bin_response[!most_freq_idx] <- response[!most_freq_idx][1]
      # }
      
      ## Start timing for individual linear combination time measurement
      tic()
      
      ## Find a linear combination
      if (splitmethod == "ride0") {
        
        res <- ride0(data_values, response)
        
      } else if (splitmethod == "ridgeauto") {
        
        res <- ridgeauto(data_values, response)
        
      } else if (splitmethod == "ridge1e10") {
        
        res <- ridge1e10(data_values, response)
        
      }
      
      ## Restrict coefficients to norm 1 for a unique solution up to factor -1
      coef_norm <- norm(as.matrix(res[-1]), "F")
      value <- res[1]/coef_norm
      coefficients <- res[-1]/coef_norm
      
      ## Restrict coefficients to positive split value for a completely unique solution
      if (value < 0) {
        value <- -value
        coefficients <- -coefficients
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
      
      if (splitrule == "Variance") {
        ## Decrease of impurity TODO
        decrease <- 
      } else {
        stop("Unknown splitrule for regression var cluster trees.")
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
    },
    
    estimate = function(nodeID) {      
      ## Return mean response
      return(mean(data$subset(sampleIDs[[nodeID]], 1)))
    }, 
    
    getNodePrediction = function(nodeID) {
      return(split_values[nodeID])
    }, 
    
    # ATTENTION: renamed from predictionError
    OOBpredictionError = function(pred = NULL) {
      if (is.null(pred)) {
        pred <- predictOOB()
      }
      sum((pred - data$subset(oob_sampleIDs, 1))^2, na.rm = TRUE) / length(oob_sampleIDs)
    })
    
)