
##' @title Tree class
##' @description Virtual class for Random forest tree where split variables are clusters. 
##' Contains all fields and methods used in all "VarClusters" tree subclasses.
TreeVarClusters <- setRefClass("TreeVarClusters",
  fields = list(
    unordered_factors = "character",
    ## Tree-specific parameters
    mtry = "integer", 
    min_node_size = "integer",
    splitrule = "character",
    data = "Data",
    IQR_data ="Data",
    sampleIDs = "list",
    oob_sampleIDs = "integer",
    child_nodeIDs = "list",
    ## Split-specific parameters
    split_clusterIDs = "integer",
    # split_selectedVarIDs = "integer",
    split_values = "numeric",
    split_coefficients = "list",
    ## Module parameters
    varclusters = "list",
    splitobject = "character",
    splitmethod = "character",
    varselection = "character",
    ## Performance Parameters
    linearcomb_times = "list",
    node_times = "numeric",
    depths = "integer",
    sizes = "integer"
    ),
  methods = list(
    
    ## Initiate bootstrap
    ## @splitNode recursively
    grow = function(replace) {
      ## Bootstrap 
      num_samples <- data$nrow
      if (replace) {
        num_bootstrap_samples <- num_samples
      } else {
        num_bootstrap_samples <- num_samples * 0.6321
      }
      bootstrap_sample <- sample(num_samples, num_bootstrap_samples, replace = replace)
      oob_sampleIDs <<- (1:num_samples)[-bootstrap_sample]
      
      ## Assign bootstrap samples to root node
      sampleIDs <<- list(bootstrap_sample)
      
      ## Call recursive splitting function on root node
      splitNode(1)
    },
    
    ## Handle node splitting
    ## @splitNode
    ## @splitNodeInternal
    ## @estimate
    splitNode = function(nodeID) {
      
      ## Possible split clusters, maximum mtry
      if (length(varclusters) <= mtry) {
        possible_split_clusterIDs <- 1:length(varclusters)
      } else {
        possible_split_clusterIDs <- sample(x = 1:length(varclusters),
                                            size = mtry,
                                            replace = FALSE)
      }
      
      ## Split node
      split <- splitNodeInternal(nodeID, possible_split_clusterIDs)

      print(split)
      print(split$node_time)
      
      ## Read performance metrics
      linearcomb_times[[nodeID]] <<- split$linearcomb_times
      node_times[nodeID] <<- split$node_time
      sizes[nodeID] <<- length(sampleIDs[[nodeID]])
      
      ## Calculate node depth
      depth <- 1
      current_nodeID <- nodeID
      while(current_nodeID != 1) {
        current_nodeID <- go_to_parent(child_nodeIDs, current_nodeID)
        depth <- depth + 1
      }
      depths[nodeID] <<- as.integer(depth)
      
      if (!is.null(split$clusterID)) {
        ## Assign split
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
        
        ## Recursively call split node on child nodes
        splitNode(left_child)
        splitNode(right_child)
      } else {
        ## Compute estimate for terminal node
        split_clusterIDs[[nodeID]] <<- NA
        # split_selectedVarIDs[[nodeID]] <<- NA
        split_coefficients[[nodeID]] <<- NA
        split_values[[nodeID]] <<- estimate(nodeID)
      }
    },
    
    ## virtual
    splitNodeInternal = function(nodeID, possible_split_varIDs) { 
      ## references subclass
    },
    
    ## virtual
    estimate = function(nodeID) {
      ## references subclass
    }, 
    
    ## virtual
    getNodePrediction = function(nodeID) {
      ## references subclass
    },
    
    ## predict on data with the tree
    ## @getNodePrediction
    predict = function(predict_data) {
      ## Initialize
      num_samples_predict <- predict_data$nrow
      predictions <- list()
      
      ## For each sample start in root and drop down tree
      for (i in 1:num_samples_predict) {
        nodeID <- 1
        while(TRUE) {
          ## Break if terminal node
          if (nodeID > length(child_nodeIDs) || is.null(child_nodeIDs[[nodeID]])) {
            break
          }
          
          ## Move to child
          if (varselection == "none") {
            value <- as.matrix(data$subset(i, varclusters[[split_clusterIDs[nodeID]]] + 1)) %*% split_coefficients[[nodeID]]
          } else {
            value <- as.matrix(data$subset(i, split_selectedVarIDs[[nodeID]] + 1)) %*% split_coefficients[[nodeID]]
          }
          if (value <= split_values[nodeID]) {
            nodeID <- child_nodeIDs[[nodeID]][1]
          } else {
            nodeID <- child_nodeIDs[[nodeID]][2]
          }
        }
        
        ## Add to prediction
        predictions[[i]] <- getNodePrediction(nodeID)
      }
      return(simplify2array(predictions))
    }, 
    
    ## predict OOB data with thee tree
    ## @getNodePrediction
    predictOOB = function() {
      ## Initialize
      num_samples_predict <- length(oob_sampleIDs)
      predictions <- list()
      
      ## For each sample start in root and drop down tree
      for (i in 1:num_samples_predict) {
        nodeID <- 1
        while(TRUE) {
          ## Break if terminal node
          if (nodeID > length(child_nodeIDs) || is.null(child_nodeIDs[[nodeID]])) {
            break
          }
          
          ## Move to child
          if (length(split_levels_left[[nodeID]]) == 0) {
            ## Ordered splitting
            if (varselection == "none") {
              value <- as.matrix(data$subset(i, varclusters[[split_clusterIDs[nodeID]]] + 1)) %*% split_coefficients[[nodeID]]
            } else {
              value <- as.matrix(data$subset(i, split_selectedVarIDs[[nodeID]] + 1)) %*% split_coefficients[[nodeID]]
            }
            if (value <= split_values[nodeID]) {
              nodeID <- child_nodeIDs[[nodeID]][1]
            } else {
              nodeID <- child_nodeIDs[[nodeID]][2]
            }
          } else {
            ## Unordered splitting
            value <- data$subset(oob_sampleIDs[i], split_varIDs[nodeID])
            if (value %in% split_levels_left[[nodeID]]) {
              nodeID <- child_nodeIDs[[nodeID]][1]
            } else {
              nodeID <- child_nodeIDs[[nodeID]][2]
            }
          }
        }
        
        ## Add to prediction
        predictions[[i]] <- getNodePrediction(nodeID)
      }
      return(simplify2array(predictions))
    }, 
    
    ## virtual
    ## @predictOOB
    predictionError = function(pred = NULL) {
      ## references subclass
    },
    
    ## permute and predict OOB data with the tree
    ## @getNodePrediction
    permuteAndPredictOOB = function(permuted_varID) {
      ## Initialize
      num_samples_predict <- length(oob_sampleIDs)
      predictions <- list()
      permutations <- sample(num_samples_predict)

      ## For each sample start in root and drop down tree
      for (i in 1:num_samples_predict) {
        nodeID <- 1
        while(TRUE) {
          ## Break if terminal node
          if (nodeID > length(child_nodeIDs) || is.null(child_nodeIDs[[nodeID]])) {
            break
          }

          ## Move to child
          if (split_varIDs[nodeID] == permuted_varID) {
            value <- as.matrix(data$subset(oob_sampleIDs[permutations[i]], varclusters[[split_clusterIDs[nodeID]]] + 1)) %*% split_coefficients[[nodeID]]
          } else {
            value <- as.matrix(data$subset(oob_sampleIDs[i], varclusters[[split_clusterIDs[nodeID]]] + 1)) %*% split_coefficients[[nodeID]]
          }
          if (value <= split_values[nodeID]) {
            nodeID <- child_nodeIDs[[nodeID]][1]
          } else {
            nodeID <- child_nodeIDs[[nodeID]][2]
          }
        }

        ## Add to prediction
        predictions[[i]] <- getNodePrediction(nodeID)
      }
      return(simplify2array(predictions))
    },
    
    ## calculate variable importance via prediction errors
    ## @predictionError
    ## @permuteAndPredictOOB
    variableImportance = function(type = "permutation") {
      if (type == "permutation") {
        
        ## Prediction error without any permutation
        oob_error <- predictionError()

        ## For each variable, prediction error after permutation
        res <- sapply(2:data$ncol, function(varID) {
          pred <- permuteAndPredictOOB(varID)
          oob_error_perm <- predictionError(pred)
          oob_error_perm - oob_error
        })
        names(res) <- data$names[-1]
        res
        
      } else if (type == "surrogate_splits") {
        
        ## 
        
      } else {
        stop("Unknown variable importance type")
      }
    }
  )
)
