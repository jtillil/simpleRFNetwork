
##' @title Classification tree class
##' @description Subclass for classification tree where split variables are clusters.
##' Contains all fields and methods used special for classification trees.
##' @import e1071
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
        
        ## Read data values from samples in nodes of current 
        data_values <- data$subset(sampleIDs[[nodeID]], varclusters[[split_clusterID]] + 1)
        
        ## Select variables
        if (varselection!="none") {
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
        best_split = findBestSplit(split_clusterID, data_values, best_split, response)
        
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
    
    ## Find best Gini split for clusters via SVM
    findBestSplit = function(split_clusterID, data_values, best_split, response) {
      
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
      
      ## Optimize linear combination
      if (splitmethod == "SVM_linear") {
        ## Calculate SVM plane
        svmfit <- svm(y=bin_response, x=data_values, kernel="linear", scale=FALSE)
        ## Read coefficients and value
        coefficients <- drop(t(svmfit$coefs)%*%as.matrix(data_values)[svmfit$index,])
        value <- svmfit$rho
      } else if (splitmethod == "Gini_optimal") {
        ## Calculate gini-optimal plane
        par <- optim(
          par=runif(ncol(data_values) + 1),
          fn=function(par) {
            select_idx <- as.matrix(data_values) %*% par[2:length(par)] > par[1]
            N1 <- sum(select_idx)
            N2 <- sum(!select_idx)
            gini <- -(
                    (sum(response[select_idx] == "1")/N1)^2+
                    (sum(response[select_idx] == "0")/N1)^2+
                    (sum(response[!select_idx] == "1")/N2)^2+
                    (sum(response[!select_idx] == "0")/N2)^2)
            return(gini)
          },
          method="L-BFGS-B"
        )$par
        ## Read coefficients and value
        coefficients <- par[2:length(par)]
        value <- par[1]
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
    })
)