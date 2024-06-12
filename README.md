## simpleRFNetwork
Johannes Tillil

### Description
R package for the publication "Tillil, Hu: Group-based random forest for disease module detection". Based upon the package simpleRF by Marvin N. Wright. Contains an implementation of RF that handles node splitting with a linear combination of multiple variables instead of a single one. Find example code for usage in the descriptions of the functions exported by this package.

### Publication
To reproduce the results of "Tillil, Hu: Group-based random forest for disease module detection", run "MAIN.R" in the "/publication" folder. Due to the nature of the Boruta variable selection method to grow a larger number of forests and the unoptimized implementation in R, executing the entire code unparallelized will take months. Reproduction is recommended on a dedicated calculation server.

To just reproduce the figures, run "generate_figures.R" in the "/publication" folder. All necessary output files are provided in "/publication/results".

### Functions
#### simpleRFNetwork
Generates a group-based random forest object that can be used for prediction and data inference. The Boruta variable selection procedure can be used to detect influential groups of predictors for a certain phenotype.

#### genGeneNetworkData
Generates random networks with gene expression data and draws a binary label for each observation using a logistic regression model.

Be careful, this package is not extensively tested!
