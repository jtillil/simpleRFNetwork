boruta <- function(dat, splitmethod, num_trees, num_threads) {
  
  while (not_all_modules_classified) {
    rf = simpleRFNetwork(pheno ~ .,
                         data = dat$data[501:1000,],
                         num_trees=num_trees,
                         num_threads=num_threads,
                         splitobject="module",
                         splitmethod=splitmethod,
                         # alternative splitmethods: univariate_fast, CART_fast, LDA, SVM, Nelder, SANN
                         varselection="none",
                         mtry="root",
                         varclusters = dat$modules
    )
    
  }
  
  
  return(rf)
}