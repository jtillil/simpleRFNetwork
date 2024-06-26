
##' @title Classification tree class
##' @description Subclass for classification tree where split variables are clusters.
##' Contains all fields and methods used special for classification trees.
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

      ## Get response
      response <- data$subset(sampleIDs[[nodeID]], 1)
      
      ## Stop if node is pure
      if (length(unique(response)) <= 1) {
        return(NULL)
      }

      ## IF LDA: stop if node has only one observation left for one class
      if (splitmethod == "LDA" & (length(response[response == 1]) <= 1 | length(response[response == 0]) <= 1)) {
        return(NULL)
      }
      
      ## Find best split
      return(findBestSplit(nodeID, possible_split_clusterIDs, response))
    }, 
    
    ## Find best split
    ## @findBestSplitCoefs
    findBestSplit = function(nodeID, possible_split_clusterIDs, response) {
      
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
      
      if (splitmethod != "Gini_stoch_optimal") {
        ## For all possible variable clusters
        for (split_clusterID in possible_split_clusterIDs) {
          
          ## Read data values from samples in current node
          data_values <- data$subset(sampleIDs[[nodeID]], varclusters[[split_clusterID]] + 1)
          
          # print(varclusters[[split_clusterID]] + 1)
          # print(response)
          # print(response == 0)
          # print(data_values)
          # print(data_values[response == 0,])
          # print(as.matrix(data_values[response == 0,]))
          # print(cova(as.matrix(data_values[response == 0,]), center=TRUE, large=FALSE))
          
          ## Calculate mean of both covariance matrices due to homoscedasticity
          mat <- 0.5*(cova(as.matrix(data_values[response == 0,]), center=TRUE, large=FALSE) +
                        cova(as.matrix(data_values[response == 1,]), center=TRUE, large=FALSE))
          
          ## Skip if bad matrix
          if (!is.square.matrix(mat) | any(is.na(mat))) {
            next
          }
          
          ## Condition matrix by adding 1e-10 to diagonal elements that are 0
          sapply(1:ncol(mat), function(j) {
            if (mat[j,j] == 0) {
              mat[j,j] <<- 1e-10
            }
          })
          
          ## Skip if singular
          if (!is.positive.definite(mat)) {
            next
          }
          
          ## IF LDA_weighted calculate different covariance matrix and skip if singular
          if (splitmethod == "LDA_weighted") {
            ## Calculate class sizes
            N0 <- sum(response == 0)
            N1 <- sum(response == 1)
            N <- N0 + N1

            ## Calculate mean of both covariance matrices due to homoscedasticity
            if (N0 != 1 & N1 != 1) {
              mat <- (N0/N) * cova(as.matrix(data_values[response == 0,]), center=TRUE, large=FALSE) +
                     (N1/N) * cova(as.matrix(data_values[response == 1,]), center=TRUE, large=FALSE)
            } else if (N0 == 1) {
              mat <- cova(as.matrix(data_values[response == 1,]), center=TRUE, large=FALSE)
            } else if (N1 == 1) {
              mat <- cova(as.matrix(data_values[response == 0,]), center=TRUE, large=FALSE)
            } else {
              next
            }

            ## Skip if matrix contains NAs
            if (any(is.na(mat))) {
              next
            }

            ## Condition matrix by adding 1e-10 to diagonal elements that are 0
            # sapply(1:ncol(mat), function(j) {
            #   if (mat[j,j] == 0) {
            #     mat[j,j] <<- 1e-10
            #   }
            # })

            ## Skip if matrix is singular
            if (!is.positive.definite(mat)) {
              next
            }
          }
          
          ## IF CART: read IQR scaled data values from samples in current node
          if (splitmethod == "CART" | splitmethod == "CART_fast") {
            IQR_data_values <- Data$new(data = IQR_data$subset(sampleIDs[[nodeID]], varclusters[[split_clusterID]]))
          } else {
            IQR_data_values <- NULL
          }
          
          ## Select variables
          # if (varselection=="half_lowest_p" | varselection=="signif_p") {
          #   ## Obtain p values from logistic regression
          #   p_vals <- lapply(varclusters[[split_clusterID]] + 1,
          #                    function(x) {
          #                      p <- summary(glm(response ~ data_values[,x],
          #                                       family=binomial(link="logit")))$coefficients[2,4] 
          #                    })
          #   if (varselection=="half_lowest_p") {
          #     ## Get ranks of variables, sorted by p value
          #     ranks <- order(p_vals)
          #     ## Use only the variables with below average p value
          #     best_split$selectedVarIDs <- varclusters[[split_clusterID]][ranks[1:round(length(ranks)/2)]]
          #   } else if (varselection=="signif_p") {
          #     ## Use only the variables with significant p value
          #     best_split$selectedVarIDs <- varclusters[[split_clusterID]][p_vals < 0.15]
          #   }
          #   data_values <- data_values[,best_split$selectedVarIDs]
          # }
          
          ## Find best split
          best_split = findBestSplitCoefs(split_clusterID, best_split, data_values, IQR_data_values, response, mat)
        
          ## Save time measurement for single linear combination
          linearcomb_times <- c(linearcomb_times, best_split$linearcomb_time)
        }
      } else {
        ## Init lists
        data_valList <- NULL
        responseList <- NULL

        ## For all possible variable clusters
        for (id in 1:length(possible_split_clusterIDs)) {
          ## Read data values from samples in current node
          data_valList[[id]] <- as.matrix(data$subset(sampleIDs[[nodeID]], varclusters[[possible_split_clusterIDs[id]]] + 1))
          ## Put in response
          responseList[[id]] <- as.numeric(response) - 1
        }

        ## Find best split in batch model training
        best_split <- findBestSplitCoefs_batch(
          split_clusterIDList = possible_split_clusterIDs,
          best_splitList = vector("list", length(possible_split_clusterIDs)),
          data_valList = data_valList,
          responseList = responseList
        )

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
    },
    
    ## Find coefficients for linear combination of variables
    findBestSplitCoefs = function(split_clusterID, best_split, data_values, IQR_data_values, response, mat) {
      
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
      if (splitmethod == "univariate") {
        
        res <- univariate_split(data_values, response)
        
      } else if (splitmethod == "univariate_fast") {
        
        res <- univariate_split_fast(data_values, response)
        
      } else if (splitmethod == "SVM") {
        
        res <- SVM(data_values, response)
        
      } else if (splitmethod == "SVM_nonparametric") {
        
        stop("SVM_nonparametric not implemented yet.")
        
      } else if (splitmethod == "SVM_Gini") {
        
        stop("SVM_Gini not implemented yet.")
        
      } else if (splitmethod == "LDA" | splitmethod == "LDA_weighted") {
        
        res <- LDA(data_values, response, mat)
        
      } else if (splitmethod == "QDA") {
        
        stop("QDA not implemented yet.")
        
      } else if (splitmethod == "Nelder") {
        
        res <- Nelder(data_values, response)
        
      } else if (splitmethod == "SANN") {
        
        res <- SANN(data_values, response)
        
      } else if (splitmethod == "NN") {
        
        res <- stoch_optim(data_values, response)
        
      } else if (splitmethod == "CART") {
        
        res <- CART(IQR_data_values, data_values, response)
        
      } else if (splitmethod == "CART_fast") {
        
        res <- CART_fast(IQR_data_values, data_values, response)
        
      } else if (splitmethod == "logridge0") {
        
        res <- logridge0(data_values, response)
        
      } else if (splitmethod == "logridge1") {
        
        res <- logridge1(data_values, response)
        
      } else if (splitmethod == "logridgeauto") {
        
        res <- logridgeauto(data_values, response)
        
      } else if (splitmethod == "logridge1e10") {
        
        res <- logridge1e10(data_values, response)
        
      } else if (splitmethod == "PCA") {
        
        res <- PCA(data_values, response)
        
      }
      
      if (!is.numeric(res)) {
        browser()
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
    OOBPredictionErrorTree = function(pred = NULL) {
      if (is.null(pred)) {
        pred <- predictOOB()
      }
      sum(pred != as.numeric(data$subset(oob_sampleIDs, 1)), na.rm = TRUE) / length(oob_sampleIDs)
    }
  )
)
