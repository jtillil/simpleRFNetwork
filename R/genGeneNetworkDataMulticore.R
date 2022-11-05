##' Generates random networks of genes, their associations, expression data and
##' binary labels where the effect measures of individual genes can be set by
##' the user.
##' 
##' @title simpleRFNetwork
##' @param num_networks Integer, number of networks to generate.
##' @param num_genes Integer, number of genes per network.
##' @param num_modules Integer, number of modules per network. Can also be NULL for random number of modules.
##' @param num_observations Integer, number of expression data observations per network.
##' @param num_causal_modules Integer, how many modules should be causal for the phenotype.
##' @param num_causal_genes Integer, how many genes in each causal module should be causal for the phenotype. Can also be "all".
##' @param effect_size Float, standardized effect size of causal genes in the network.
##' @param effect_intercept Float, standardized intercept effect of all genes in the network.
##' @param causal_genes_randomly_distributed Boolean, if TRUE causal genes will be sampled randomly from all available genes.
##' @param num_threads Integer, number of cores to parallelize on.
##' @param seed Integer, initial seed for the L'Ecuyer-CMRG random number streams.
##' @examples
##' \donttest{
##' library(simpleRFNetwork)
##' 
##' 
##' }
##' 
##' @author Johannes Tillil
##' @import stats
##' @import SeqNet
##' @import parallel
##' @export
genGeneNetworkDataMulticore <- function(
  num_networks,
  num_genes,
  num_modules = NULL,
  num_observations,
  num_causal_modules,
  num_causal_genes,
  effect_size = 1,
  effect_intercept = -1,
  causal_genes_randomly_distributed = FALSE,
  num_threads = 1,
  seed = 1
) {

  ## Set up parallel reproducibility
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  mc.reset.stream()
  
  ## Start parallel computing
  return(mclapply(
    1:num_networks,
    function(i) {
      ## Generate networks, modules and associations between genes
      if (is.null(num_modules)) {
        rn <- random_network(num_genes)
      } else {
        rn <- random_network(num_genes, num_modules)
      }
      rn <- gen_partial_correlations(rn)
      exprdat <- scale(log(gen_rnaseq(num_observations, rn)$x + 1))
      modules <- rn$modules
      num_modules <- length(rn$modules)

      ## Sample causal modules and genes and set gene effects accordingly
      if (causal_genes_randomly_distributed) {
        ## Sample causal genes
        causal_genes <- sample(
          1:num_genes,
          num_causal_genes,
          replace = FALSE
        )
        ## Count causal modules
        causal_modules <- NULL
        sapply(
          1:length(modules),
          function(j) {
            if (sum(causal_genes %in% modules[[j]]$nodes) > 0) {
              causal_modules <<- c(causal_modules, j)
            }
          }
        )
        ## Save sampled values
        causal_genes <- unique(causal_genes)
        effects <- numeric(num_genes)
        effects[causal_genes] <- effect_size
        effects <- effects
        causal_modules <- causal_modules
        causal_genes <- causal_genes
      } else if (num_causal_modules > 0) {
        ## Sample causal modules
        causal_modules <- sample(
          x = 1:num_modules,
          size = min(num_causal_modules, num_modules),
          replace = FALSE)
        ## Sample causal genes
        causal_genes <- NULL
        if (num_causal_genes == "all") {
          sapply(1:min(num_causal_modules, num_modules),
                function(j) {
                  causal_genes <<- c(
                    causal_genes,
                    modules[[causal_modules[j]]]$nodes)
                })
        } else {
          sapply(
            1:min(num_causal_modules, num_modules),
            function(j) {
              causal_genes <<- c(
                causal_genes,
                sample(
                  x = modules[[causal_modules[j]]]$nodes,
                  size = num_causal_genes,
                  replace = FALSE))
            })
        }
        ## Save sampled values
        causal_genes <- unique(causal_genes)
        effects <- numeric(num_genes)
        effects[causal_genes] <- effect_size
        effects <- effects
        causal_modules <- causal_modules
        causal_genes <- causal_genes
      } else {
        causal_modules <- NULL
        causal_genes <- NULL
        effects <- numeric(num_genes)
      }
      
      ## Sample phenotype from gene effects and combine phenotype and expression data into one data frame for training
      probs <- 1/(1 + exp(-as.matrix(exprdat) %*% effects - effect_intercept))
      res <- data.frame(pheno = as.factor(sapply(
        1:num_observations,
        function(i) {
          sample(
            x = c(0,1),
            size = 1,
            prob = c(1-probs[i], probs[i]))
        }
      )))
      res <- cbind(res, data.frame(exprdat))
      return(list(
        data = res, 
        modules = lapply(modules, function(x){x$nodes}),
        causal_modules = causal_modules,
        causal_genes = causal_genes,
        effects = effects
      ))
    },
    mc.cores = num_threads
  ))
}
