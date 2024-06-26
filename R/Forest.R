
##' @title Forest class
##' @description Virtual class for Random forest. 
##' Contains all fields and methods used in all Forest subclasses.
Forest <- setRefClass("Forest", 
  fields = list(
    ## Standard parameters
    num_trees = "integer", 
    mtry = "integer", 
    min_node_size = "integer", 
    splitrule = "character",
    unordered_factors = "character",
    data = "Data",
    IQR_data = "Data",
    predict_data = "Data",
    formula = "formula",
    trees = "list",
    treetype = "character",
    replace = "logical", 
    covariate_levels = "list",
    forest_time = "numeric",
    ## Module parameters
    varclusters = "list",
    splitobject = "character",
    splitmethod = "character",
    varselection = "character",
    seed = "integer"),
  methods = list(
    
    grow = function(num_threads) { 

      ## Start timing for forest growth
      tic()
      
      ## Init trees
      temp <- lapply(trees, function(x) {
        ## Standard parameters
        x$mtry <- mtry
        x$min_node_size <- min_node_size
        x$splitrule <- splitrule
        x$unordered_factors <- unordered_factors
        x$data <- data
        x$IQR_data <- IQR_data
        ## Module parameters
        x$varclusters <- varclusters
        x$splitobject <- splitobject
        x$splitmethod <- splitmethod
        x$varselection <- varselection
      })

      # if (splitmethod == "Gini_stoch_optimal") {

      #   ## Divide trees into batches
      #   treeIDs <- 1:num_trees
      #   treeBatchList <- NULL
      #   i <- 0
      #   while(length(treeIDs) > 0) {
      #     i <- i+1
      #     treeBatchList[[i]] <- treeIDs[1:num_threads]
      #     treeIDs <- treeIDs[-(1:num_threads)]#
      #   }
      #   treeBatchList[[i]] <- treeBatchList[[i]][!is.na(treeBatchList[[i]])]

      #   ## Grow trees in batches
      #   ## @grow_batch
      #   for (treeBatch in treeBatchList) {
      #     grow_batch(treeBatch)
      #   }

      # } else {
        ## Set up parallel reproducibility
        RNGkind("L'Ecuyer-CMRG")
        set.seed(seed)
        mc.reset.stream()
        
        ## Grow trees
        if (Sys.info()["sysname"]=="Windows") {
          ## On Windows
          cl <- makeCluster(num_threads)
          trees <<- parLapply(cl, X=trees, fun=function(x) {
            x$grow(replace)
            x
          })
        } else {
          ## On Unix
          trees <<- mclapply(trees, function(x) {
            x$grow(replace)
            x
          }, mc.cores = num_threads)
        }
      # }
      
      # trees <<- lapply(trees, function(x) {
      #   x$grow(replace)
      #   x
      # })

      ## Stop timing for forest growth
      forest_time <- toc(quiet = TRUE)
      forest_time <<- as.numeric(forest_time$toc - forest_time$tic)

    },
    
    predict = function(newdata = data, num_threads = 1) {
      ## Save prediction data in model
      # predict_data <<- Data$new(data = newdata)
      
      ## Predict in trees
      if (Sys.info()["sysname"]=="Windows") {
        ## On Windows
        cl <- makeCluster(num_threads)
        predictions <- parLapply(cl, X=trees, fun=function(x) {
          x$predict(newdata)
        })
      } else {
        ## On Unix
        predictions <- simplify2array(mclapply(trees, function(x) {
          x$predict(newdata)
        }, mc.cores = num_threads))
      }
      
      ## Aggregate predictions
      return(aggregatePredictions(predictions))
    }, 
    
    aggregatePredictions = function(predictions) {
      ## Empty virtual function
    }, 
    
    predictionErrorForest = function() {
      ## Empty virtual function
    },

    predictionErrorTrees = function() {
      ## Empty virtual function
    },

    predictionErrorForestAndTrees = function() {
      ## Empty virtual function
    },
    
    variableImportance = function(type = "permutation", num_threads = 1) {
      ## Calculate tree VIM
      
      if (Sys.info()["sysname"]=="Windows") {
        ## On Windows
        cl <- makeCluster(num_threads)
        vim_trees <- parLapply(cl, X=trees, fun=function(x) {
          x$variableImportance(type)
        })
      } else {
        ## On Unix
        vim_trees <- mclapply(trees, function(x) {
          x$variableImportance(type)
        }, mc.cores = num_threads)
      }

      # vim_trees
      
      ## Aggregate over trees
      return(rowMeans(simplify2array(vim_trees)))
    },

    readSplitClusters = function(num_threads = 1) {
      ## Read split_clusterIDs
      return(mclapply(
        trees,
        function(tree) tree$split_clusterIDs,
        mc.cores = 1
      ))
    },

    readSplitValues = function(num_threads = 1) {
      ## Read split_values
      return(mclapply(
        trees,
        function(tree) tree$split_values,
        mc.cores = 1
      ))
    },

    readSplitCoefficients = function(num_threads = 1) {
      ## Read split_coefficients
      return(mclapply(
        trees,
        function(tree) tree$split_coefficients,
        mc.cores = 1
      ))
    },

    readDepths = function(num_threads = 1) {
      ## Read depths
      return(mclapply(
        trees,
        function(tree) tree$depths,
        mc.cores = 1
      ))
    },

    readSizes = function(num_threads = 1) {
      ## Read sizes
      return(mclapply(
        trees,
        function(tree) tree$sizes,
        mc.cores = 1
      ))
    },

    readLinearcombTimes = function(num_threads = 1) {
      ## Read linearcomb times
      return(mclapply(
        trees,
        function(tree) tree$linearcomb_times,
        mc.cores = 1
      ))
    },

    readImpurities = function(num_threads = 1) {
      ## Read linearcomb times
      return(mclapply(
        trees,
        function(tree) tree$impurities,
        mc.cores = 1
      ))
    },

    readOOBSamples = function(num_threads = 1) {
      ## Read OOB samples
      return(mclapply(
        trees,
        function(tree) tree$oob_sampleIDs,
        mc.cores = 1
      ))
    },

    readChildNodes = function(num_threads = 1) {
      ## Read child nodes
      return(mclapply(
        trees,
        function(tree) tree$child_nodeIDs,
        mc.cores = 1
      ))
    },
    
    show = function() {
      cat("simpleRF Forest\n")
      cat("Type:                            ", treetype, "\n")
      cat("Splitrule:                       ", splitrule, "\n")
      cat("Splitmethod:                     ", splitmethod, "\n")
      cat("Variable Selection:              ", varselection, "\n")
      cat("Number of trees:                 ", num_trees, "\n")
      cat("Sample size:                     ", data$nrow, "\n")
      cat("Number of independent variables: ", data$ncol-1, "\n")
      cat("Mtry:                            ", mtry, "\n")
      cat("Target node size:                ", min_node_size, "\n")
      cat("Replace                          ", replace, "\n")
      # cat("Unordered factor handling        ", unordered_factors, "\n")
      # cat("OOB prediction error:            ", predictionError(), "\n")
      cat("Forest execution time:           ", forest_time, "\n")
    }, 
    
    print = function() {
      show()
    })
)
