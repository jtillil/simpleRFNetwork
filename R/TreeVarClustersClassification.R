
##' @title Classification tree class
##' @description Subclass for classification tree where split variables are clusters.
##' Contains all fields and methods used special for classification trees.
##' @import e1071
##' @import Rfast
##' @import keras
##' @import tensorflow
##' @import tictoc
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
      # best_split$selectedVarIDs <- -1
      best_split$coefficients <- -1
      best_split$value <- -1
      best_split$decrease <- -1
      best_split$linearcomb_time <- -1
      
      ## Get response
      response <- data$subset(sampleIDs[[nodeID]], 1)
      
      ## Initialize array for individual linear combination time measurement
      linearcomb_times <- c()
      
      ## Start timing for node splitting time measurement
      tic()
      
      ## For all possible variable clusters
      for (i in 1:length(possible_split_clusterIDs)) {
        
        ## Set current cluster ID
        split_clusterID <- possible_split_clusterIDs[i]
        
        ## Read data values from samples in current node
        data_values <- data$subset(sampleIDs[[nodeID]], varclusters[[split_clusterID]] + 1)
        
        ## Read IQR scaled data values from samples in current node if needed
        if (splitmethod == "CART" | splitmethod == "CART_fast") {
          IQR_data_values <- IQR_data$subset(sampleIDs[[nodeID]], varclusters[[split_clusterID]])
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
        best_split = findBestSplitCoefs(split_clusterID, best_split)
        
        ## Save time measurement for single linear combination
        linearcomb_times <- c(linearcomb_times, best_split$linearcomb_time)
        
        ## Assign split_levels_left for compatibility with cluster-less version
        if (best_split$clusterID == split_clusterID) {
          split_levels_left[[nodeID]] <<- list()
        }
      }
      
      ## Stop timing for node splitting time measurement
      node_time <- toc(quiet = TRUE)
      
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
        result$linearcomb_times <- linearcomb_times
        result$node_time <- as.numeric(node_time$toc - node_time$tic)
        return(result)
      }
    },
    
    ## Find Gini-optimal coefficients for linear combination of variables
    findBestSplitCoefs = function(split_clusterID, best_split) {
      
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
        
        ## Initiate
        Gini_impurity_start <- 9999
        best_val <- 0
        best_varID <- 0
        
        ## Iterate over all variables
        sapply(1:ncol(data_values), function(varID) {
          ## Read value candidates
          val_candidates <- data_values[,varID]
          sapply(val_candidates, function(val) {
            ## Compute new Gini impurity
            coefficients_start <- numeric(ncol(data_values))
            coefficients_start[varID] <- 1
            Gini_impurity_val <- gini_impurity(data_values,
                                               response,
                                               c(val, coefficients_start))
            
            ## Compare to current best Gini impurity and set value, varID if smaller
            if (Gini_impurity_val < Gini_impurity_start) {
              Gini_impurity_start <<- Gini_impurity_val
              best_val <<- val
              best_varID <<- varID
            }
          })
        })
        
        ## Coerce best uni-variate split into coefficients and value
        coefficients <- numeric(ncol(data_values))
        coefficients[best_varID] <- 1
        value <- best_val
        
      } else if (splitmethod == "univariate_fast") {
        
        ## Set fraction of subset variables
        nu <- 0.1
        
        ## Initiate
        Gini_impurity_start <- 9999
        best_val <- 0
        best_varID <- 0
        
        ## Iterate over all variables
        sapply(1:ncol(data_values), function(varID) {
          ## Sample value candidates
          val_candidates <- sample(data_values[,varID], round(nu*nrow(data_values)))
          sapply(val_candidates, function(val) {
            ## Compute new Gini impurity
            coefficients_start <- numeric(ncol(data_values))
            coefficients_start[varID] <- 1
            Gini_impurity_val <- gini_impurity(data_values,
                                               response,
                                               c(val, coefficients_start))
            
            ## Compare to current best Gini impurity and set value, varID if smaller
            if (Gini_impurity_val < Gini_impurity_start) {
              Gini_impurity_start <<- Gini_impurity_val
              best_val <<- val
              best_varID <<- varID
            }
          })
        })
        
        ## Coerce best uni-variate split into coefficients and value
        coefficients <- numeric(ncol(data_values))
        coefficients[best_varID] <- 1
        value <- best_val
        
      } else if (splitmethod == "SVM_linear") {
        
        ## Calculate SVM plane
        svmfit <- svm(y=bin_response,
                      x=data_values,
                      kernel="linear",
                      scale=FALSE)
        
        ## Read coefficients and value
        coefficients <- drop(t(svmfit$coefs)%*%as.matrix(data_values)[svmfit$index,])
        value <- svmfit$rho
        
      } else if (splitmethod == "SVM_nonparametric") {
        
        stop("SVM_nonparametric not implemented yet.")
        
      } else if (splitmethod == "SVM_Gini") {
        
        stop("SVM_Gini not implemented yet.")
        
      } else if (splitmethod == "LDA") {
        
        ## Calculate class means
        mean0 <- colmeans(as.matrix(data_values[response == 0,]))
        mean1 <- colmeans(as.matrix(data_values[response == 1,]))
        
        ## Calculate coefficients and value
        ## Calculate mean of both covariance matrices due to homoscedasticity
        coefficients <- spdinv(0.5*(cova(as.matrix(data_values[response == 1,])) +
                                    cova(as.matrix(data_values[response == 0,])))) %*% (mean1 - mean0)
        value <- sum(coefficients * (0.5*(mean1 + mean0)))
        
      } else if (splitmethod == "QDA") {
        
        stop("QDA not implemented yet.")
        
      } else if (splitmethod == "Gini_optimal") {
        
        ## Find first split as best uni-variate split
        ## Set fraction of subset variables
        nu <- 0.1
        
        ## Initiate
        Gini_impurity_start <- 9999
        best_val <- 0
        best_varID <- 0
        
        ## Iterate over all variables
        sapply(1:ncol(data_values), function(varID) {
          ## Sample value candidates
          val_candidates <- sample(data_values[,varID], round(nu*nrow(data_values)))
          sapply(val_candidates, function(val) {
            ## Compute new Gini impurity
            coefficients_start <- numeric(ncol(data_values))
            coefficients_start[varID] <- 1
            Gini_impurity_val <- gini_impurity(data_values,
                                               response,
                                               c(val, coefficients_start))
            
            ## Compare to current best Gini impurity and set value, varID if smaller
            if (Gini_impurity_val < Gini_impurity_start) {
              Gini_impurity_start <<- Gini_impurity_val
              best_val <<- val
              best_varID <<- varID
            }
          })
        })
        
        ## Coerce best uni-variate split into parameters
        coefficients <- numeric(ncol(data_values))
        coefficients[best_varID] <- 1
        param <- c(best_val, coefficients)
        
        ## Calculate Gini-optimal plane with univariate initial split
        par <- optim(
          par = param,
          fn = function(par) {
            return(gini_impurity(data_values,
                                 response,
                                 par))
          },
          method="Nelder"
        )$par
        
        ## Read coefficients and value
        coefficients <- par[2:length(par)]
        value <- par[1]
        
      } else if (splitmethod == "Gini_stoch_optimal") {
        
        ## Create, compile and fit Keras model
        model <- keras_model_sequential() %>%
          layer_dense(units = 1,
                      activation = "sigmoid",
                      input_shape = ncol(data_values))
        model %>%
          compile(optimizer = "adam",
                  # loss = "binary_crossentropy",
                  loss = gini_loss,
                  run_eagerly = TRUE)
        model %>%
          fit(as.matrix(data_values),
              as.numeric(response),
              epochs = 100,
              verbose = 0)
        
        ## Read coefficients and value
        coefficients <- as.numeric(model$layers[[1]]$weights[[1]]$numpy())
        value <- as.numeric(model$layers[[1]]$weights[[2]]$numpy())
        
      } else if (splitmethod == "CART") {
        
        ## Use data_values centered around 0 and divided by their interquartile range
        ## IQR_data_values
        
        ## CARE FOR WHEN IQR VALS GET SUBSETTED
        
        ## CARE THAT WITHOUT CONSIDERATION OF THE U SEARCH ADMISSION TO CHILDS IS EITHER LEFT OR RIGHT, NOT ALWAYS LEFT
        
        ## Find first split as best uni-variate split
        ## Initiate
        Gini_impurity_start <- 9999
        best_val <- 0
        best_varID <- 0
        
        ## Iterate over all variables
        sapply(1:IQR_data_values$ncol, function(varID) {
          ## Read value candidates
          val_candidates <- IQR_data_values$column(varID)
          sapply(val_candidates, function(val) {
            ## Compute new Gini impurity
            coefficients_start <- numeric(IQR_data_values$ncol)
            coefficients_start[varID] <- 1
            Gini_impurity_val <- gini_impurity(IQR_data_values$data,
                                               response,
                                               c(val, coefficients_start))
            
            ## Compare to current best Gini impurity and set value, varID if smaller
            if (Gini_impurity_val < Gini_impurity_start) {
              Gini_impurity_start <<- Gini_impurity_val
              best_val <<- val
              best_varID <<- varID
            }
          })
        })
        
        ## Coerce best uni-variate split into coefficients and value
        coefficients <- numeric(IQR_data_values$ncol)
        coefficients[best_varID] <- 1/best_val
        value <- 1
        
        ## Compute starting Gini impurity
        Gini_impurity_nplus1 <- gini_impurity(IQR_data_values$data,
                                              response,
                                              c(value, coefficients))
        Gini_impurity_n <- 9999
        
        ## Set convergence threshold
        epsilon <- 0.01
        
        ## Perform updates until improvement below threshold
        while ((Gini_impurity_n - Gini_impurity_nplus1) > epsilon) {
          ## Set Gini impurity of last cycle
          Gini_impurity_n <- Gini_impurity_nplus1
          
          ## Cycle through all variables and search for an improved split by varying their coefficient
          sapply(1:IQR_data_values$ncol,
                 function(varID) {
                   ## Compute current split values for all observations
                   v <- as.matrix(IQR_data_values$data) %*% coefficients
                   
                   ## For gamma equals -0.25, 0, 0.25
                   sapply(c(-0.25, 0, 0.25), function(gamma) {
                     ## Compute candidates
                     u <- (v - value) / (IQR_data_values$column(varID) + gamma)
                     
                     ## For every candidate
                     sapply(u, function(u_n) {
                       ## Convert to candidate coefficients and value
                       coefficients_u <- coefficients
                       coefficients_u[varID] <- coefficients_u[varID] - u_n
                       value_u <- value + u_n * gamma
                       
                       if (gini_impurity(IQR_data_values$data,
                                         response,
                                         c(value_u, coefficients_u)
                       ) < Gini_impurity_nplus1) {
                         
                         ## Update coefficients, value, split values
                         coefficients <<- coefficients_u
                         value <<- value_u
                         v <<- as.matrix(IQR_data_values$data) %*% coefficients
                         
                         ## Calculate new Gini impurity
                         Gini_impurity_nplus1 <<- gini_impurity(IQR_data_values$data,
                                                                response,
                                                                c(value, coefficients))
                       }
                     })
                   })
                 }
          )
        }
        
        ## Rescale coefficients and value
        axis_points <- value/coefficients*IQR_vals
        translated_axis_points <- NULL
        for (varID in 1:length(coefficients)) {
          translated_axis_points[varID] <- axis_points[varID] + Mean_vals[varID] + sum(axis_points[varID]/axis_points[-varID] * Mean_vals[-varID])
        }
        coefficients <- 1/translated_axis_points
        value <- 1
        
      } else if (splitmethod == "CART_fast") {
        
        ## Only change compared to CART splitmethod is sampling of observations during update steps to reduce search time
        
        ## Use data_values centered around 0 and divided by their interquartile range
        ## IQR_data_values
        
        ## CARE FOR WHEN IQR VALS GET SUBSETTED
        
        ## CARE THAT WITHOUT CONSIDERATION OF THE U SEARCH ADMISSION TO CHILDS IS EITHER LEFT OR RIGHT, NOT ALWAYS LEFT
        
        ## Find first split as best uni-variate split
        ## Set fraction of subset variables
        nu <- 0.1
        
        ## Initiate
        Gini_impurity_start <- 9999
        best_val <- 0
        best_varID <- 0
        
        sapply(1:IQR_data_values$ncol, function(varID) {
          ## Sample value candidates
          val_candidates <- sample(IQR_data_values$column(varID), round(nu*IQR_data_values$nrow))
          sapply(val_candidates, function(val) {
            ## Compute new Gini impurity
            coefficients_start <- numeric(IQR_data_values$ncol)
            coefficients_start[varID] <- 1
            Gini_impurity_val <- gini_impurity(IQR_data_values$data,
                                               response,
                                               c(val, coefficients_start))
            
            ## Compare to current best Gini impurity and set value, varID if smaller
            if (Gini_impurity_val < Gini_impurity_start) {
              Gini_impurity_start <<- Gini_impurity_val
              best_val <<- val
              best_varID <<- varID
            }
          })
        })
        
        ## Coerce best uni-variate split into coefficients and value
        coefficients <- numeric(IQR_data_values$ncol)
        coefficients[best_varID] <- 1/best_val
        value <- 1
        
        ## Compute starting Gini impurity
        Gini_impurity_nplus1 <- gini_impurity(IQR_data_values$data,
                                              response,
                                              c(value, coefficients))
        Gini_impurity_n <- 9999
        
        ## Set convergence threshold
        epsilon <- 0.01
        
        ## Perform updates until improvement below threshold
        while ((Gini_impurity_n - Gini_impurity_nplus1) > epsilon) {
          ## Set Gini impurity of last cycle
          Gini_impurity_n <- Gini_impurity_nplus1
          
          ## Cycle through all variables and search for an improved split by varying their coefficient
          sapply(1:IQR_data_values$ncol,
                 function(varID) {
                   ## Sample subset
                   subset_data_values <- IQR_data_values$data[sample(1:nrow(data_values), round(nu*nrow(data_values))),]
                   
                   ## Compute current split values for subset observations
                   v <- as.matrix(subset_data_values) %*% coefficients
                   
                   ## For gamma equals -0.25, 0, 0.25
                   sapply(c(-0.25, 0, 0.25), function(gamma) {
                     ## Compute candidates
                     u <- (v - value) / (subset_data_values[,varID] + gamma)
                     
                     ## For every candidate
                     sapply(u, function(u_n) {
                       ## Convert to candidate coefficients and value
                       coefficients_u <- coefficients
                       coefficients_u[varID] <- coefficients_u[varID] - u_n
                       value_u <- value + u_n * gamma
                       
                       if (gini_impurity(IQR_data_values$data,
                                         response,
                                         c(value_u, coefficients_u)
                       ) < Gini_impurity_nplus1) {
                         
                         ## Update coefficients, value, split values
                         coefficients <<- coefficients_u
                         value <<- value_u
                         v <<- as.matrix(IQR_data_values$data) %*% coefficients
                         
                         ## Calculate new Gini impurity
                         Gini_impurity_nplus1 <<- gini_impurity(IQR_data_values$data,
                                                                response,
                                                                c(value, coefficients))
                       }
                     })
                   })
                 }
          )
        }
        
        ## Rescale coefficients and value
        axis_points <- value/coefficients*IQR_vals
        translated_axis_points <- NULL
        for (varID in 1:length(coefficients)) {
          translated_axis_points[varID] <- axis_points[varID] + Mean_vals[varID] + sum(axis_points[varID]/axis_points[-varID] * Mean_vals[-varID])
        }
        coefficients <- 1/translated_axis_points
        value <- 1
        
      }
      
      ## Restrict coefficients to norm 1 for a unique solution
      coef_norm <- Norm(as.matrix(coefficients), "F")
      value <- value/coef_norm
      coefficients <- coefficients/coef_norm
      
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
      
      ## Use this svm split if decrease better than from earlier clusters
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
    predictionError = function(pred = NULL) {
      if (is.null(pred)) {
        pred <- predictOOB()
      }
      sum(pred != as.numeric(data$subset(oob_sampleIDs, 1)), na.rm = TRUE) / length(oob_sampleIDs)
    }
  )
)
