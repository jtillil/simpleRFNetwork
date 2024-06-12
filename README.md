## simpleRFNetwork
Johannes Tillil

### Description
R package for the master thesis of Johannes Tillil. Based upon the package simpleRF by Marvin N. Wright. Contains an implementation of RF that handles node splitting with a linear combination of multiple variables instead of a single one. Find example code for usage in the descriptions of the functions exported by this package.

### Publication
To reproduce the results of "Tillil, Hu: Group-based random forest for disease module detection", run "MAIN.R" in the "/publication" folder.
To just reproduce the figures, run "generate_figures.R" in the "/publication" folder. All necessary output files are provided in "/publication/results".

### Functions
#### simpleRFNetwork
Generates a RF object that can be used for prediction and data inference using the grouped variable random forest algorithm proposed in the master thesis of Johannes Tillil.

#### genGeneNetworkData
Generates random networks with gene expression data and draws a binary label for each observation using a logistic regression model.

Be careful, this package is not extensively tested!
