##' Implements Random Forests (Breiman 2001) with emphasis on simplicity. 
##' Uses reference classes and only plain \code{R}. 
##' Not optimized for computation speed. 
##' Allows rapid prototyping of RF-type algorithms.
##' 
##' @title simpleRFNetwork
##' @param formula Object of class \code{formula} or \code{character} describing the model to fit.
##' @param data Training data of class \code{data.frame}.
##' @param num_trees Number of trees.
##' @param mtry Number of variables to possibly split at in each node. One of "root", "half", "full".
##' @param min_node_size Minimal node size. Default 1 for classification, 5 for regression, 3 for survival and 10 for probability estimation.
##' @param replace Sample with replacement. Default TRUE.
##' @param probability Grow a probability forest. Default FALSE.
##' @param splitrule Splitrule to use in trees. Default "Gini" for classification and probability forests, "Variance" for regression forests and "Logrank" for survival forests.
##' @param splitobject Object for node splitting. Use "single_variable" for standard procedure and "module" for use of modules as vector-valued predictor variables in splitting.
##' @param splitmethod Method for node splitting if splitobject is "module". One of "univarate", "univariate_fast", "CART", "CART_fast", "LDA", "SVM", "Nelder", "SANN".
##' @param varselection Variable selection in case multiple variables are used for node splitting. Default "none" uses all variables of a module.
##' @param varclusters List of numeric vectors that contain the IDs of genes in the respective modules.
##' @param unordered_factors How to handle unordered factor variables. One of "ignore", "order_once", "order_split" and "partition" with default "ignore".
##' @param num_threads Number of threads used for mclapply. Default is 1.
##' @examples
##' \donttest{
##' library(simpleRFNetwork) 
##' library(SeqNet)
##' 
##' # Generate Network Data
##' 
##' testdat <- genGeneNetworkData(
##'   num_networks = 1,
##'   num_genes = 500,
##'   num_modules = NULL,
##'   max_genes_per_module = 100,
##'   sd_genes_per_module = 25,
##'   num_observations = 1500,
##'   num_causal_modules = 1,
##'   prop_causal_genes = 0.5,
##'   total_effect_size = 10,
##'   effect_intercept = -1,
##'   causal_genes_randomly_distributed = FALSE,
##'   num_threads = 12,
##'   seed = 1
##' )
##' 
##' # Create Random Forest
##' 
##' rf <- simpleRFNetwork(pheno ~ .,
##'                       data=testdat[[1]]$data[1001:1500,],
##'                       num_trees=50,
##'                       num_threads=1,
##'                       splitobject="module",
##'                       splitmethod="LDA",
##' # alternative splitmethods: univariate_fast, CART_fast, LDA, SVM, Nelder, SANN
##'                       varselection="none",
##'                       mtry="root",
##'                       varclusters=testdat[[1]]$modules)
##' 
##' # Forest Diagnostics
##' 
##' rf$forest_time
##' 
##' rf$trees[[1]]$varclusters
##' rf$trees[[1]]$split_clusterIDs
##' 
##' rf$trees[[1]]$split_coefficients
##' rf$trees[[1]]$split_values
##' 
##' rf$trees[[1]]$linearcomb_times
##' rf$trees[[1]]$node_times
##' rf$trees[[1]]$sizes
##' rf$trees[[1]]$depths
##'
##' # Permutation Importance and OOB Prediction Error
##' 
##' rf$variableImportance(num_threads = 1)
##' rf$predictionErrorForestAndTrees()
##' 
##' # Generalization Error
##' 
##' sum(rf$predict(testdat[[1]]$data[1:1000,-1]) != testdat[[1]]$data[1:1000,1]) / 1000
##' }
##' 
##' @author Marvin N. Wright, Johannes Tillil
##' @references
##' Breiman, L. (2001). Random forests. Mach Learn, 45(1), 5-32. \cr
##' @export
simpleRFNetwork <- function(
  ## Standard parameters
  formula,
  data,
  num_trees = 50, 
  mtry = NULL, 
  min_node_size = NULL, 
  replace = TRUE,
  probability = FALSE,
  splitrule = NULL,
  unordered_factors = "ignore",
  num_threads = 1,
  ## Module Parameters
  varclusters = NULL,
  splitobject = NULL,
  splitmethod = NULL,
  varselection = NULL,
  seed = 1L
  ) {
  
  model.data <- model.frame(formula, data)
  
  ## Treetype
  if (class(model.data[, 1]) == "factor") {
    if (probability) {
      treetype <- "Probability" 
    } else {
      treetype <- "Classification"
    }
  } else if (class(model.data[, 1]) == "numeric") {
    treetype <- "Regression"
  } else if (class(model.data[, 1]) == "Surv") {
    treetype <- "Survival"
  } else {
    stop("Unkown response type.")
  }
  
  ## Check parameters
  if (mtry == "full") {
    mtry <- length(varclusters)
  } else if (mtry == "half") {
    mtry <- ceiling(length(varclusters) * 0.5)
  } else if (mtry == "root") {
    mtry <- sqrt(length(varclusters))
  } else {
    mtry <- sqrt(length(varclusters))
  }
  if (mtry > length(varclusters)) {
    stop("Mtry cannot be larger than number of variable groups.")
  }
  if (is.null(min_node_size)) {
    if (treetype == "Classification") {
      min_node_size <- 1
    } else if (treetype == "Probability") {
      min_node_size <- 10
    } else if (treetype == "Regression") {
      min_node_size <- 5
    } else if (treetype == "Survival") {
      min_node_size <- 3
    }
  }
  
  ## Splitrule
  if (is.null(splitrule)) {
    if (treetype == "Classification") {
      splitrule <- "Gini"
    } else if (treetype == "Probability") {
      splitrule <- "Gini"
    } else if (treetype == "Regression") {
      splitrule <- "Variance"
    } else if (treetype == "Survival") {
      splitrule <- "Logrank"
    }
  }
  
  ## Splitobject
  if (is.null(splitobject)) {
    splitobject <- "single_variable"
  }
  if (!(splitobject %in% c("single_variable",
                           "module"))) {
    stop("Unknown value for splitobject.")
  }
  
  ## Splitmethod
  if (is.null(splitmethod)) {
    splitmethod <- "Gini_optimal"
  }
  if (!(splitmethod %in% c(## CLASSIFICATION
                           "univariate",
                           "univariate_fast",
                           "SVM",
                           "SVM_nonparametric",    # Todo
                           "SVM_Gini",             # Todo
                           "LDA",
                           "LDA_weighted",
                           "QDA",                  # Todo
                           "Nelder",
                           "SANN",
                           "NN",
                           "CART",
                           "CART_fast",
                           ## REGRESSION
                           "logridge0",
                           "logridgeauto",
                           "logridge1e10",
                           "ridge0",
                           "ridge1",
                           "ridge1e10"))) {
    stop("Unknown value for splitmethod.")
  }
  
  ## Variable Selection
  if (treetype == "Classification") {
    if (is.null(varselection)) {
      varselection <- "none"
    }
    if (!(varselection %in% c("none",
                              "half_lowest_p",
                              "signif_p"))) {
      stop("Unknown value for varselection.")
    }
  }
  
  ## Variable Clusters
  ## 
  ## varclusters: 
  ##    list of clusters
  ##
  ## objects inside lists:
  ##    vector of variable IDs inside cluster
  
  ## Unordered factors
  if (!(unordered_factors %in% c("ignore", "order_once", "order_split", "partition"))) {
    stop("Unknown value for unordered_factors.")
  }
  covariate_levels <- list()
  
  if (unordered_factors == "order_once") {
    ## Reorder factor columns depending on response type
    model.data <- reorder.factor.columns(model.data)
    
    ## Save levels
    covariate_levels <- lapply(model.data[, -1], levels)
  }
  else if (unordered_factors == "ignore") {
    ## Just set to ordered if "ignore"
    character.idx <- sapply(model.data[, -1], is.character)
    ordered.idx <- sapply(model.data[, -1], is.ordered)
    factor.idx <- sapply(model.data[, -1], is.factor)
    recode.idx <- character.idx | (factor.idx & !ordered.idx)
    model.data[, -1][, recode.idx] <- lapply(model.data[, -1][, recode.idx], as.ordered)
    
    ## Save levels
    covariate_levels <- lapply(model.data[, -1], levels)
  }
  
  ## Create forest object
  if (treetype == "Classification") {
    ## Create data object
    dat <- Data$new(data = model.data)
    ## Create interquartile range normalized data object
    if (splitmethod == "CART" | splitmethod == "CART_fast") {
      IQR_dat <- model.data[,-1]
      means <- colMeans(as.matrix(IQR_dat))
      sapply(1:ncol(IQR_dat),
             function(i) {
               if (IQR(IQR_dat[,i])==0) {
                 IQR_dat[,i] <<- 0
               } else {
                 IQR_dat[,i] <<- (IQR_dat[,i] - means[i])/IQR(IQR_dat[,i])
               }
             })
      IQR_dat <- Data$new(data = IQR_dat)
    } else {
      IQR_dat <- Data$new(data = data.frame())
    }
    ## Create forest
    forest <- ForestClassification$new(## Standard parameters
                                       num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
                                       min_node_size = as.integer(min_node_size), 
                                       replace = replace, splitrule = splitrule,
                                       data = dat,
                                       IQR_data = IQR_dat,
                                       formula = formula, unordered_factors = unordered_factors, 
                                       covariate_levels = covariate_levels,
                                       response_levels = levels(model.data[, 1]),
                                       ## Module Parameters
                                       varclusters = varclusters,
                                       splitobject = splitobject,
                                       splitmethod = splitmethod,
                                       varselection = varselection,
                                       seed = seed)
  } else if (treetype == "Probability") {
    forest <- ForestProbability$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
                                   min_node_size = as.integer(min_node_size), 
                                   replace = replace, splitrule = splitrule,
                                   data = Data$new(data = model.data), 
                                   formula = formula, unordered_factors = unordered_factors,
                                   covariate_levels = covariate_levels,
                                   response_levels = levels(model.data[, 1]),
                                   ## Module Parameters
                                   varclusters = varclusters,
                                   splitobject = splitobject,
                                   splitmethod = splitmethod,
                                   varselection = varselection,
                                   seed = seed)
  } else if (treetype == "Regression") {
    forest <- ForestRegression$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
                                   min_node_size = as.integer(min_node_size), 
                                   replace = replace, splitrule = splitrule,
                                   data = Data$new(data = model.data), 
                                   formula = formula, unordered_factors = unordered_factors, 
                                   covariate_levels = covariate_levels,
                                   ## Module Parameters
                                   varclusters = varclusters,
                                   splitobject = splitobject,
                                   splitmethod = splitmethod,
                                   varselection = varselection,
                                   seed = seed)
  } else if (treetype == "Survival") {
    idx.death <- model.data[, 1][, 2] == 1
    timepoints <- sort(unique(model.data[idx.death, 1][, 1]))
    forest <- ForestSurvival$new(num_trees = as.integer(num_trees), mtry = as.integer(mtry), 
                                 min_node_size = as.integer(min_node_size), 
                                 replace = replace, splitrule = splitrule,
                                 data = Data$new(data = model.data), 
                                 formula = formula, unordered_factors = unordered_factors, 
                                 covariate_levels = covariate_levels,
                                 timepoints = timepoints)
  } else {
    stop("Unkown tree type.")
  }

  ## Grow forest
  forest$grow(num_threads = num_threads)

  ## Return forest
  return(forest) 
}
