
##' @title Classification tree class
##' @description Subclass for classification tree where split variables are clusters.
##' Contains all fields and methods used special for classification trees.
##' @import e1071
##' @import Rfast
TreeVarClustersClassification <- setRefClass("TreeVarClustersClassification",
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
      
      ## Stop if node is pure
      unique_response <- unique(data$subset(sampleIDs[[nodeID]], 1))
      if (length(unique_response) == 1) {
        return(NULL)
      }
      
      ## Find best split, stop if no decrease of impurity
      return(findBestSplit(nodeID, possible_split_clusterIDs))
    }, 
    
    ## Try to order factors and finds best split
    ## @findBestSplitValuePartition
    ## @findBestSplitValueOrdered
    findBestSplit = function(nodeID, possible_split_clusterIDs) {
      
      ## Initialize
      best_split <- NULL
      best_split$clusterID <- -1
      best_split$selectedVarIDs <- -1
      best_split$coefficients <- -1
      best_split$value <- -1
      best_split$decrease <- -1
      
      ## Get response
      response <- data$subset(sampleIDs[[nodeID]], 1)
      
      ## For all possible variable clusters (ignoring mtry)
      for (i in 1:length(possible_split_clusterIDs)) {
        
        ## Set current cluster ID
        split_clusterID <- possible_split_clusterIDs[i]
        
        ## Read data values from samples in current node
        data_values <- data$subset(sampleIDs[[nodeID]], varclusters[[split_clusterID]] + 1)
        
        ## Read IQR scaled data values from samples in current node if needed
        if (splitmethod == "CART") {
          IQR_data_values <- IQR_data$subset(sampleIDs[[nodeID]], varclusters[[split_clusterID]] + 1)
        }
        
        ## Select variables
        if (varselection=="half_lowest_p" | varselection=="signif_p") {
          ## Obtain p values from logistic regression
          p_vals <- lapply(varclusters[[split_clusterID]] + 1,
                           function(x) {
                             p <- summary(glm(response ~ data_values[,x],
                                              family=binomial(link="logit")))$coefficients[2,4] 
                           })
          if (varselection=="half_lowest_p") {
            ## Get ranks of variables, sorted by p value
            ranks <- order(p_vals)
            ## Use only the variables with below average p value
            best_split$selectedVarIDs <- varclusters[[split_clusterID]][ranks[1:round(length(ranks)/2)]]
          } else if (varselection=="signif_p") {
            ## Use only the variables with significant p value
            best_split$selectedVarIDs <- varclusters[[split_clusterID]][p_vals < 0.15]
          }
          data_values <- data_values[,best_split$selectedVarIDs]
        }
        
        ## Find best split
        best_split = findBestSplitCoefs(split_clusterID, data_values, best_split, response)
        
        ## Assign split_levels_left for compatibility with cluster-less version
        if (best_split$clusterID == split_clusterID) {
          split_levels_left[[nodeID]] <<- list()
        }
      }
      
      if (best_split$clusterID < 0) {
        ## Stop if no good split found
        return(NULL)
      } else {
        ## Return best split
        result <- NULL
        result$clusterID <- as.integer(best_split$clusterID)
        result$selectedVarIDs <- best_split$selectedVarIDs
        result$coefficients <- best_split$coefficients
        result$value <- best_split$value
        return(result)
      }
    },
    
    ## Find Gini-optimal coefficients for linear combination of variables
    findBestSplitCoefs = function(split_clusterID, data_values, best_split, response) {
      
      ## Coerce all but the most frequent factor level to a single one
      ## Irrelevant, if exactly two factors
      ## Skip cluster if less than two levels in response
      if (nlevels(droplevels(response)) < 2) {
        return(best_split)
      } else {
        bin_response <- response
        most_freq_idx <- response==names(which.max.random(summary(response)))
        bin_response[!most_freq_idx] <- response[!most_freq_idx][1]
      }
      
      ## Find a linear combination
      if (splitmethod == "SVM_linear") {
        
        ## Calculate SVM plane
        svmfit <- svm(y=bin_response, x=data_values, kernel="linear", scale=FALSE)
        
        ## Read coefficients and value
        coefficients <- drop(t(svmfit$coefs)%*%as.matrix(data_values)[svmfit$index,])
        value <- svmfit$rho
        
      } else if (splitmethod == "SVM_nonparametric") {
        
        stop("SVM_nonparametric not implemented yet.")
        
      } else if (splitmethod == "SVM_Gini") {
        
        stop("SVM_Gini not implemented yet.")
        
      } else if (splitmethod == "LDA") {
        
        ## Calculate class means
        mean0 <- colmeans(as.matrix(scale(data_values[response == 0,])))
        mean1 <- colmeans(as.matrix(scale(data_values[response == 1,])))
        
        ## Calculate coefficients and value
        ## Calculate mean of both covariance matrices due to homoscedasticity
        coefficients <- spdinv(0.5*(cova(as.matrix(scale(data_values[response == 1,]))) +
                                    cova(as.matrix(scale(data_values[response == 0,]))))) %*% (mean1 - mean0)
        value <- coefficients %*% (0.5*(mean1 - mean0))
        
      } else if (splitmethod == "QDA") {
        
        stop("QDA not implemented yet.")
        
      } else if (splitmethod == "Gini_optimal") {
        
        ## Calculate Gini-optimal plane
        par <- optim(
          par = c(1, rep(1/ncol(data_values), ncol(data_values))),
          fn = function(par) {
            return(gini_impurity(data_values,
                                 response,
                                 par))
          },
          method="Nelder"
        )$par
        
        ## Restrict parameters to norm 1 for a unique solution
        par <- par/norm(par, type="2")
        
        ## Read coefficients and value
        coefficients <- par[2:length(par)]
        value <- par[1]
        
      } else if (splitmethod == "CART") {
        
        ## Use data_values centered around 0 and divided by their interquartile range
        ## IQR_data_values
        
        ## CARE FOR WHEN IQR VALS GET SUBSETTED
        
        ## Compute first split as best univariate split
        split <- findBestUnivariateSplit(nodeID,
                                         varclusters[[split_clusterID]] + 1,
                                         IQR_data_values,
                                         response)
        coefficients <- numeric(IQR_data_values$ncol)
        coefficients[split$varID] <- 1/split$value
        value <- 1
        Gini_impurity_nplus1 <- gini_impurity(IQR_data_values,
                                              response,
                                              c(value, coefficients))
        Gini_impurity_n <- Gini_impurity_nplus1 + 1000
        
        ## Set threshold
        epsilon <- 0.0001
        
        ## Perform updates until improvement below threshold
        while ((Gini_impurity_n - Gini_impurity_nplus1) > epsilon) {
          ## Cycle through all variables and search for an improved split by varying their coefficient
          sapply(1:ncol(data_values),
                 function(varID) {
                   ## Compute current split values for all observations
                   v <- 
                   ## For gamma equals -0.25, 0, 0.25
                   for (gamma in c(-0.25, 0, 0.25)) {
                     ## For delta equals ...
                     for (delta in c(1,2)) {
                       ## Find best
                       print("2")
                     }
                   }
                 }
          )
          
          ## Calculate new Gini impurity
          Gini_impurity_n <- Gini_impurity_nplus1
          Gini_impurity_nplus1 <- gini_impurity(IQR_data_values,
                                                response,
                                                c(value, coefficients))
        }
        
      }
      
      ## Count classes in childs
      idx <- as.matrix(data_values)%*%coefficients <= value
      class_counts_left <- tabulate(response[idx])
      class_counts_right <- tabulate(response[!idx])
      
      ## Skip cluster if one child empty
      if (sum(class_counts_left) == 0 | sum(class_counts_right) == 0) {
        return(best_split)
      }
      
      if (splitrule == "Gini") {
        ## Decrease of impurity
        decrease <- sum(class_counts_left^2)/sum(class_counts_left) + 
          sum(class_counts_right^2)/sum(class_counts_right)
      } else {
        stop("Unknown splitrule.")
      }
      
      ## Use this svm split if decrease better than from earlier clusters
      if (decrease > best_split$decrease) {
        best_split$clusterID <- split_clusterID
        best_split$coefficients <- coefficients
        best_split$value <- value
        best_split$decrease <- decrease
      } else {
        best_split$selectedVarIDs <- -1
      }
      
      return(best_split)
    },
    
    ## Find most common label of samples in node
    estimate = function(nodeID) {
      ## Return class with maximal count, random at ties
      class_counts <- table(data$subset(sampleIDs[[nodeID]], 1))
      which.max.random(class_counts)
    },
    
    ## Return split value for a node
    ## If terminal node, this is the class prediction for that node
    getNodePrediction = function(nodeID) {
      return(split_values[nodeID])
    },
    
    ## Calculate proportion of wrong OOB predictions
    ## @predictOOB
    predictionError = function(pred = NULL) {
      if (is.null(pred)) {
        pred <- predictOOB()
      }
      sum(pred != as.numeric(data$subset(oob_sampleIDs, 1)), na.rm = TRUE) / length(oob_sampleIDs)
    },
    
    findBestUnivariateSplit = function(nodeID, possible_split_varIDs, data_values, response) {
      ## Initialize
      best_split <- NULL
      best_split$decrease <- -1
      best_split$varID <- -1
      best_split$value <- -1
      
      ## For all possible variables
      for (i in 1:length(possible_split_varIDs)) {
        split_varID <- possible_split_varIDs[i]
        data_values <- data_values$subset(sampleIDs[[nodeID]], split_varID)
        
        ## If not ordered, use partition splitting
        if (!is.numeric(data_values) & !is.ordered(data_values)) {
          best_split = findBestSplitValuePartition(split_varID, data_values, best_split, response)
          
          ## Set split levels left
          if (best_split$varID == split_varID) {
            split_levels_left[[nodeID]] <<- best_split$values_left
          }
        } else {
          best_split = findBestSplitValueOrdered(split_varID, data_values, best_split, response)
          
          ## Set split levels left (empty if ordered splitting)
          if (unordered_factors == "order_split") {
            if (best_split$varID == split_varID) {
              split_levels_left[[nodeID]] <<- unique(data_values[data_values <= best_split$value])
              
              if (is.factor(data_values)) {
                ## Use same splits as in partition
                ints <- as.integer(factor(split_levels_left[[nodeID]], levels = levels(data$subset(sampleIDs[[nodeID]], split_varID))))
                if (sum(2^(ints-1)) >= 2^(max(as.numeric(data$subset(sampleIDs[[nodeID]], split_varID))) - 1)) {
                  split_levels_left[[nodeID]] <<- unique(data_values[data_values > best_split$value])
                }
              }
            }
          } else {
            if (best_split$varID == split_varID) {
              split_levels_left[[nodeID]] <<- list()
            }
          }
        }
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
    
    findBestSplitValueOrdered = function(split_varID, data_values, best_split, response) {
      ## For all possible splits
      possible_split_values <- unique(data_values)
      for (j in 1:length(possible_split_values)) {
        split_value <- possible_split_values[j]
        
        ## Count classes in childs
        idx <- data_values <= split_value
        class_counts_left <- tabulate(response[idx])
        class_counts_right <- tabulate(response[!idx])
        
        ## Skip if one child empty
        if (sum(class_counts_left) == 0 | sum(class_counts_right) == 0) {
          next
        }
        
        if (splitrule == "Gini") {
          ## Decrease of impurity
          decrease <- sum(class_counts_left^2)/sum(class_counts_left) + 
            sum(class_counts_right^2)/sum(class_counts_right)
        } else {
          stop("Unknown splitrule.")
        }
        
        ## Use this split if better than before
        if (decrease > best_split$decrease) {
          best_split$value <- split_value
          best_split$varID <- split_varID
          best_split$decrease <- decrease
        }
      }
      return(best_split)
    },
    
    findBestSplitValuePartition = function(split_varID, data_values, best_split, response) {
      ## For all possible splits
      possible_split_values <- sort(unique(data_values))
      
      ## For all 2^(n-1)-1 2-partitions
      num_partitions <- 2^(length(possible_split_values) - 1) - 1
      for (j in 1:num_partitions) {
        ## Convert number to logic vector
        left_idx <- as.bitvect(j, length = length(possible_split_values))
        values_left <- possible_split_values[left_idx]
        
        ## Count classes in childs
        idx <- data_values %in% values_left
        class_counts_left <- tabulate(response[idx])
        class_counts_right <- tabulate(response[!idx])
        
        ## Skip if one child empty
        if (sum(class_counts_left) == 0 | sum(class_counts_right) == 0) {
          next
        }
        
        if (splitrule == "Gini") {
          ## Decrease of impurity
          decrease <- sum(class_counts_left^2)/sum(class_counts_left) + 
            sum(class_counts_right^2)/sum(class_counts_right)
        } else {
          stop("Unknown splitrule.")
        }
        
        ## Use this split if better than before
        if (decrease > best_split$decrease) {
          best_split$values_left <- values_left
          best_split$varID <- split_varID
          best_split$decrease <- decrease
        }
      }
      return(best_split)
    }
  )
)
