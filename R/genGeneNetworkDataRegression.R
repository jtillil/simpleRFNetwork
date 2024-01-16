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
genGeneNetworkDataRegression <- function(
    n_networks = 1,
    n_genes = 1000,
    n_samples = 1000,
    max_genes_per_module = 100,
    sd_genes_per_module = 25,
    disease_modules = T,
    # n_disease_modules = 2,  ALWAYS 3 DISEASE MODULES
    # main_disease_gene = F,  NO MAIN DISEASE GENES
    # average_beta = 1,       PARAMETERS PRE-SPECIFIED
    prop_disease_genes = 1,
    num_threads = 1
) {
  ## Set up parallel reproducibility
  RNGkind("L'Ecuyer-CMRG")
  set.seed(1)
  mc.reset.stream()
  
  ## Start parallel computing
  networks = mclapply(
    1:n_networks,
    function(i) {
      # seed
      set.seed(i)
      
      #### Phenotype ####
      
      # error
      e <- rnorm(n_samples, mean = 0, sd = 0.1)
      
      if (disease_modules) {
        # while counter
        while_idx = 1
        
        # ensure enough modules
        while (TRUE) {
          # generate network
          network <- random_network(
            n_genes, 
            max_module_size = max_genes_per_module,
            sd_module_size = sd_genes_per_module)
          network <- gen_partial_correlations(network)
          
          # generate RNA-Seq data and take log-transformation
          x.total <- gen_rnaseq(n_samples, network)
          # TODO: log2??? changed from log2 to log
          x.total <- log(x.total$x + 1)
          x.total <- scale(x.total, center = TRUE, scale = TRUE)
          
          modules <- lapply(network$modules, function(x){x$nodes})
          num_modules <- length(modules)
          
          # disease module candidates
          module.length <- lengths(modules)
          mod.len.q1 <- floor(quantile(module.length, probs = 0.25))
          mod.candidate <- which(module.length <= mod.len.q1)
          if (length(mod.candidate) < 3) {
            while_idx = while_idx + 1
            print(paste("while iteration:", while_idx))
            next
          }
          
          # first module
          mod.signal.1st <- sample(mod.candidate, 1)
          causal_modules = mod.signal.1st
          gene.sig <- network$modules[[mod.signal.1st]]$nodes
          
          # second module
          mod.intersect <- sapply(mod.candidate,
                                  function(mod){length(intersect(gene.sig, network$modules[[mod]]$nodes))>0})
          mod.mutual <- mod.candidate[!mod.intersect]
          if (length(mod.mutual) < 2) {
            while_idx = while_idx + 1
            print(paste("while iteration:", while_idx))
            next
          }
          mod.signal.2nd <- sample(mod.mutual, 1)
          causal_modules = c(causal_modules, mod.signal.2nd)
          gene.sig <- c(gene.sig, network$modules[[mod.signal.2nd]]$nodes)
          
          # third module
          mod.intersect <- sapply(mod.candidate,
                                  function(mod){length(intersect(gene.sig, network$modules[[mod]]$nodes))>0})
          mod.mutual <- mod.candidate[!mod.intersect]
          if (length(mod.mutual) < 1) {
            while_idx = while_idx + 1
            print(paste("while iteration:", while_idx))
            next
          }
          mod.signal.3rd <- sample(mod.mutual, 1)
          causal_modules = c(causal_modules, mod.signal.3rd)
          gene.sig <- c(gene.sig, network$modules[[mod.signal.3rd]]$nodes)
          
          break
        }
        
        # sample causal genes
        causal_genes = lapply(
          1:3,
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
            adj_mat <- get_adjacency_matrix(network)[modules[[causal_modules[j]]], modules[[causal_modules[j]]]]
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
            return(modules[[causal_modules[j]]][sampled_genes])
          }
        )
        
        # read disease gene data
        x.disease <- lapply(causal_genes,
                            function(genes) rowMeans(x.total[, genes]))
        x.disease <- as.data.frame(x.disease)
        
        # regression phenotype
        pheno <- 0.25 * exp(4*x.disease[, 1]) + 4 / (1+exp(-20*(x.disease[, 2]-0.5))) + 3*x.disease[, 3] + e
        
        # return
        return(list(
          data = cbind(pheno, as.data.frame(x.total)),
          modules = modules,
          causal_modules = causal_modules,
          causal_genes = causal_genes
        ))
      } else {
        # generate network
        network <- random_network(n_genes)
        network <- gen_partial_correlations(network)
        
        # generate RNA-Seq data and take log-transformation
        x.total <- gen_rnaseq(n_samples, network)
        x.total <- log2(x.total$x + 1)
        x.total <- scale(x.total, center = TRUE, scale = TRUE)
        
        # random data not from the gene network
        x.disease <- matrix(runif(n_samples*3, 0, 1), ncol = 3)
        
        # regression phenotype
        pheno <- 0.25 * exp(4*x.disease[, 1]) + 4 / (1+exp(-20*(x.disease[, 2]-0.5))) + 3*x.disease[, 3] + e
        
        # return
        return(list(
          data = cbind(pheno, as.data.frame(x.total)),
          modules = lapply(network$modules, function(x) x$nodes),
          causal_modules = c(),
          causal_genes = c()
        ))
      }
    },
    mc.cores = num_threads
  )
}
