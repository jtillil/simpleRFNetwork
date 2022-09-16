##' Generates random networks of genes, their associations, expression data and
##' binary labels where the effect measures of individual genes can be set by
##' the user.
##' 
##' @title simpleRFNetwork
##' @param num_networks Integer, number of networks to generate.
##' @param num_genes Integer, number of genes per network.
##' @param num_modules Integer, number of modules per network.
##' @param num_observations Integer, number of expression data observations per network.
##' @param num_causal_modules Integer, how many modules should be causal for the phenotype.
##' @param num_causal_genes Integer, how many genes in each causal module should be causal for the phenotype. Can also be "all".
##' @param effect_measure Float, standardized effect measure of causal genes in the network.
##' @examples
##' \donttest{
##' library(simpleRFNetwork)
##' library(SeqNet) 
##' 
##' 
##' }
##' 
##' @author Johannes Tillil
##' @references
##' Breiman, L. (2001). Random forests. Mach Learn, 45(1), 5-32. \cr
##' @import stats
##' @import SeqNet
##' @export
genGeneNetworkData <- function(
  num_networks,
  num_genes,
  num_modules,
  num_observations,
  num_causal_modules,
  num_causal_genes,
  effect_measure
  ) {
  
  ## Set seed
  set.seed(1)
  
  ## Generate networks, modules and associations between genes
  networkdat <- lapply(1:num_networks,
                       function(i) {
                         rn <- random_network(num_genes, num_modules)
                         rn <- gen_partial_correlations(rn)
                         exprdat <- gen_rnaseq(num_observations, rn)
                         return(list(
                           exprdat = exprdat,
                           modules = rn$modules
                         ))
                       })
  
  ## Sample causal modules and genes and set gene effects accordingly
  if (num_causal_modules > 0) {
    sapply(1:num_networks,
           function(i) {
             effects <- numeric(num_genes)
             causal_modules <- sample(x = 1:num_modules,
                                      size = num_causal_modules,
                                      replace = FALSE)
             causal_genes <- c()
             if (num_causal_genes == "all") {
               sapply(causal_modules,
                      function(j) {
                        causal_genes <<- c(causal_genes,
                                           networkdat[[i]]$modules[[j]]$nodes)
                      })
             } else {
               sapply(causal_modules,
                      function(j) {
                        causal_genes <<- c(causal_genes,
                                           sample(x = networkdat[[i]]$modules[[j]]$nodes,
                                                  size = num_causal_genes,
                                                  replace = FALSE))
                      })
             }
             causal_genes <- unique(causal_genes)
             effects[causal_genes] <- effect_measure
             networkdat[[i]]$effects <<- effects
           })
  } else {
    sapply(1:num_networks,
           function(i) {
             networkdat[[i]]$effects <<- numeric(num_genes)
           })
  }
  
  ## Sample phenotype from gene effects and combine phenotype and expression data into one data frame for training
  return(lapply(1:num_networks,
                function(i) {
                  probs <- 1/(1 + exp(-scale(as.matrix(networkdat[[i]]$exprdat$x)) %*% networkdat[[i]]$effects))
                  res <- data.frame(pheno = as.factor(sapply(
                    1:num_observations,
                    function(i) {
                      sample(x = c(0,1),
                             size = 1,
                             prob = c(1-probs[i], probs[i]))
                    }
                  )))
                  res <- cbind(res, data.frame(networkdat[[i]]$exprdat$x))
                  return(res)
                }))
  
}