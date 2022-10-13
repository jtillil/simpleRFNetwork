
##' @title Forest class
##' @description Virtual class for Random forest. 
##' Contains all fields and methods used in all Forest subclasses.
##' @importFrom parallel mclapply
##' @importFrom parallel makeCluster
##' @importFrom parallel parLapply
##' @import methods
##' @import tictoc
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
    varselection = "character"),
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
      
      # trees <<- lapply(trees, function(x) {
      #   x$grow(replace)
      #   x
      # })

      ## Stop timing for forest growth
      forest_time <- toc(quiet = TRUE)
      forest_time <<- as.numeric(forest_time$toc - forest_time$tic)

    }, 
    
    predict = function(newdata) {
      model.data <- model.frame(formula, newdata)

      ## Recode factors if forest grown 'order_once' mode
      if (unordered_factors == "order_once" & length(covariate_levels) > 0) {
        model.data[, -1] <- mapply(function(x, y) {
          if(is.null(y)) {
            x
          } else {
            new.levels <- setdiff(levels(x), y)
            factor(x, levels = c(y, new.levels), ordered = TRUE)
          }
        }, model.data[, -1], covariate_levels, SIMPLIFY = FALSE)
      }

      ## Save prediction data in model
      predict_data <<- Data$new(data = model.data)
      
      ## Predict in trees
      predictions <- simplify2array(lapply(trees, function(x) {
        x$predict(predict_data)
      }))
      
      ## Aggregate predictions
      return(aggregatePredictions(predictions))
    }, 
    
    aggregatePredictions = function(predictions) {
      ## Empty virtual function
    }, 
    
    predictionError = function() {
      ## Empty virtual function
    },
    
    variableImportance = function(type = "permutation", num_threads = 1) {
      ## Calculate tree VIM
      vim_trees <- mclapply(trees, function(x) {
        x$variableImportance(type)
      }, mc.cores = num_threads)
      
      ## Aggregate over trees
      rowMeans(simplify2array(vim_trees))
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
      cat("OOB prediction error:            ", predictionError(), "\n")
      cat("Forest execution time:           ", forest_time, "\n")
    }, 
    
    print = function() {
      show()
    })
)
