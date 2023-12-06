##' Generates random networks of genes, their associations, expression data and
##' binary labels where the effect measures of individual genes can be set by
##' the user.
##' 
##' @title genGeneNetworkDataClassification
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
##' testdat <- genGeneNetworkDataClassification(
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
genGeneNetworkDataClassification <- function(
    n_networks = 1,
    n_genes = 1000,
    n_samples = 1000,
    max_genes_per_module = 100,
    sd_genes_per_module = 25,
    n_disease_modules = 2,
    main_disease_gene = F,
    average_beta = 1,
    num_threads = 1
) {
  ## Start parallel computing
  networks = mclapply(
    1:n_networks,
    function(i) {
      # seed
      set.seed(i)
      
      # generate network
      network <- random_network(
        n_genes,
        max_module_size = max_genes_per_module,
        sd_module_size = sd_genes_per_module
      )
      network <- gen_partial_correlations(network)
      
      # generate RNA-Seq data and take log-transformation
      x.total <- gen_rnaseq(n_samples, network)
      x.total <- log(x.total$x + 1)
      x.total <- scale(x.total, center = TRUE, scale = TRUE)
      
      # disease module candidates
      module.length <- sapply(network$modules, function(module) length(module$nodes))
      mod.len.q1 <- floor(quantile(module.length, probs = 0.25))
      mod.candidate <- which(module.length <= mod.len.q1)
      
      # first module
      mod.signal.1st <- sample(mod.candidate, 1)
      mod.sig = mod.signal.1st
      gene.sig <- network$modules[[mod.signal.1st]]$nodes
      
      # second module
      if(n_disease_modules > 1) {
        mod.intersect <- sapply(mod.candidate,
                                function(mod){length(intersect(gene.sig, network$modules[[mod]]$nodes))>0})
        mod.mutual <- mod.candidate[!mod.intersect]
        mod.signal.2nd <- sample(mod.mutual, 1)
        mod.sig = c(mod.sig, mod.signal.2nd)
        gene.sig <- c(gene.sig, network$modules[[mod.signal.2nd]]$nodes)
      }
      
      #### sample phenotype ####
      
      # read disease gene data
      x.disease = x.total[, gene.sig]
      
      # binary phenotype
      pheno = as.factor(rbinom(n_samples, 1, sigmoid(average_beta*rowMeans(x.disease) - 1)))
      
      # return
      return(list(
        data = cbind(pheno, as.data.frame(x.total)),
        modules = lapply(network$modules, function(x) x$nodes),
        causal_modules = mod.sig,
        causal_genes = gene.sig
      ))
    },
    mc.cores = num_threads
  )
}
