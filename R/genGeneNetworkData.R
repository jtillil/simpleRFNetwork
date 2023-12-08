##' Generates random networks of genes, their associations, expression data and
##' binary labels where the effect measures of individual genes can be set by
##' the user.
##' 
##' @title genGeneNetworkData
##' @param num_networks Integer, number of networks to generate.
##' @param num_genes Integer, number of genes per network.
##' @param num_modules Integer, number of modules per network. Can also be NULL for random number of modules.
##' @param max_genes_per_module Integer, maximum module size.
##' @param sd_genes_per_module Float, standard deviation of module sizes.
##' @param num_observations Integer, number of expression data observations per network.
##' @param num_causal_modules Integer, how many modules should be causal for the phenotype.
##' @param prop_causal_genes Float, proportion of genes in each causal module that should be causal for the phenotype. If "causal_genes_randomly_distributed == TRUE" refers to proportion of all genes that should be causal.
##' @param total_effect_size Float, standardized effect size of causal genes in the network.
##' @param effect_intercept Float, standardized intercept effect of all genes in the network.
##' @param causal_genes_randomly_distributed Boolean, if TRUE causal genes will be sampled randomly from all available genes.
##' @param num_threads Integer, number of cores to parallelize on.
##' @param seed Integer, initial seed for the L'Ecuyer-CMRG random number streams.
##' @param effect_type Character, type of gene-effect on phenotype. Default is "linear". One of "linear", "quadratic".
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
##' }
##' 
##' @author Johannes Tillil
##' @export
# genGeneNetworkData <- function(
#   num_networks,
#   num_genes,
#   num_modules = NULL,
#   max_genes_per_module,
#   sd_genes_per_module,
#   num_observations,
#   num_causal_modules,
#   prop_causal_genes,
#   total_effect_size = 10,
#   effect_intercept = -1,
#   effect_type = "linear",
#   # effect_error_sd = 0,
#   causal_genes_randomly_distributed = FALSE,
#   num_threads = 1,
#   seed # DEPRECATED, NOT USED
# ) {
genGeneNetworkData <- function(
  n_networks = 1,
  n_genes = 1000,
  n_samples = 1000,
  max_genes_per_module = 100,
  sd_genes_per_module = 25,
  n_disease_modules = 2,
  main_disease_gene = F,
  prop_disease_genes = 1,
  average_beta = 1,
  effect_intercept = -1,
  num_threads = 1
) {
  ## Set up parallel reproducibility
  RNGkind("L'Ecuyer-CMRG")
  set.seed(1)
  mc.reset.stream()
  
  ## Start parallel computing
  return(mclapply(
    1:n_networks,
    function(i) {
      set.seed(i)
      
      ## Generate networks, modules and associations between genes
      # if (is.null(n_modules)) {
        rn <- random_network(
          n_genes, 
          max_module_size = max_genes_per_module,
          sd_module_size = sd_genes_per_module)
      # } else {
      #   rn <- random_network(
      #     num_genes, 
      #     num_modules, 
      #     max_module_size = max_genes_per_module,
      #     sd_module_size = sd_genes_per_module)
      # }
      rn <- gen_partial_correlations(rn)
      exprdat <- scale(log(gen_rnaseq(n_samples, rn)$x + 1))
      modules <- lapply(rn$modules, function(x){x$nodes})
      num_modules <- length(modules)

      ## Sample causal modules and genes and set gene effects accordingly
      # if (causal_genes_randomly_distributed) {
      #   ## Sample causal genes
      #   causal_genes <- sample(
      #     1:num_genes,
      #     ceiling(prop_causal_genes*num_genes),
      #     replace = FALSE
      #   )
      #   ## Count causal modules
      #   causal_modules <- NULL
      #   sapply(
      #     1:length(modules),
      #     function(j) {
      #       if (sum(causal_genes %in% modules[[j]]) > 0) {
      #         causal_modules <<- c(causal_modules, j)
      #       }
      #     }
      #   )
      #   ## Save sampled values
      #   causal_genes <- unique(causal_genes)
      #   effects <- numeric(num_genes)
      #   effects[causal_genes] <- 1
      #   causal_modules <- causal_modules
      #   causal_genes <- causal_genes
      # } else if (num_causal_modules > 0 & prop_causal_genes > 0) {
      
      ## UPDATED sample disease modules
      # disease module candidates
      module.length <- lengths(modules)
      mod.len.q1 <- floor(quantile(module.length, probs = 0.25))
      mod.candidate <- which(module.length <= mod.len.q1)
      
      # first module
      mod.signal.1st <- sample(mod.candidate, 1)
      causal_modules = mod.signal.1st
      gene.sig <- modules[[mod.signal.1st]]
      
      # second module
      if(n_disease_modules > 1) {
        mod.intersect <- sapply(mod.candidate,
                                function(mod){length(intersect(gene.sig, modules[[mod]]))>0})
        mod.mutual <- mod.candidate[!mod.intersect]
        mod.signal.2nd <- sample(mod.mutual, 1)
        causal_modules = c(causal_modules, mod.signal.2nd)
      }
      
      ## LEGACY Sample causal modules
      # if (sum(lengths(modules) < 30) >= min(num_causal_modules, num_modules)) {
      #   causal_modules <- sample(
      #     x = (1:num_modules)[lengths(modules) < 30],
      #     size = min(num_causal_modules, num_modules),
      #     replace = FALSE)
      # } else {
      #   causal_modules <- order(lengths(modules))[1:min(num_causal_modules, num_modules)]
      # }
      
      ## Sample causal genes
      causal_genes <- NULL
      sapply(
        1:min(num_causal_modules, num_modules),
        function(j) {
          ## Read required amount of genes
          num_required_genes <- ceiling(prop_disease_genes*length(modules[[causal_modules[j]]]))
          ## Sample first causal gene
          sampled_genes <- sample(
            x = 1:length(modules[[causal_modules[j]]]),
            size = 1,
            replace = FALSE
          )
          ## Read adjacency matrix for module
          adj_mat <- get_adjacency_matrix(rn)[modules[[causal_modules[j]]], modules[[causal_modules[j]]]]
          ## Search for associated genes in the module
          while (length(sampled_genes) < num_required_genes) {
            ## For all candidates not yet added to causal genes
            for (candidate_id in sample((1:nrow(adj_mat))[-sampled_genes])) {
              ## If still required AND connected to sampled_genes
              if (
                length(sampled_genes) < num_required_genes &
                sum(adj_mat[candidate_id, sampled_genes]) > 0
              ) {
                sampled_genes <- c(sampled_genes, candidate_id)
              }
            }
          }
          causal_genes <<- c(causal_genes, modules[[causal_modules[j]]][sampled_genes])
        }
      )
      
      ## Save sampled values
      causal_genes <- unique(causal_genes)
      effects <- numeric(n_genes)
      effects[causal_genes] <- average_beta
      
      # } else {
      #   causal_modules <- NULL
      #   causal_genes <- NULL
      #   effects <- numeric(num_genes)
      # }
      
      ## LEGACY Sample phenotype from gene effects and combine phenotype and expression data into one data frame for training
      # if (sum(effects) != 0) {
      #   effects <- effects * (total_effect_size / sum(effects))
      # }
      # # if (effect_error_sd > 0) {
      # #   effects <- effects + rnorm(length(effects), 0, effect_error_sd)
      # # }
      # if (effect_type == "linear") {
      #   probs <- 1/(1 + exp(-as.matrix(exprdat) %*% effects - effect_intercept))
      # } else if (effect_type == "quadratic") {
      #   probs <- 1/(1 + exp(-(sign(as.matrix(exprdat) %*% effects) * (as.matrix(exprdat) %*% effects)^2) - effect_intercept))
      # } else if (effect_type == "cubic") {
      #   probs <- 1/(1 + exp(-((as.matrix(exprdat) %*% effects)^3) - effect_intercept))
      # } else if (effect_type == "root") {
      #   probs <- 1/(1 + exp(-(sign(as.matrix(exprdat) %*% effects) * sqrt(abs(as.matrix(exprdat) %*% effects))) - effect_intercept))
      # }
      # res <- data.frame(pheno = as.factor(sapply(
      #   1:num_observations,
      #   function(i) {
      #     sample(
      #       x = c(0,1),
      #       size = 1,
      #       prob = c(1-probs[i], probs[i]))
      #   }
      # )))
      
      probs <- 1/(1 + exp(-as.matrix(exprdat) %*% effects - effect_intercept))
      res <- data.frame(pheno = as.factor(sapply(
        1:n_samples,
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
        modules = modules,
        causal_modules = causal_modules,
        causal_genes = causal_genes,
        effects = effects
      ))
    },
    mc.cores = num_threads
  ))
}
