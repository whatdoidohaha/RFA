# Reference-free analysis of genetic architecture

Here we provide 

Reference-free analysis of genetic architecture



## Examples

We provide scripts and tutorials for four types of analyses.

**RFA:** Here, reference-free effects on genetic score and nonspecific epistasis are jointly inferred. Any form of genotype space, including when the number of states varies among sites, is supported. We perform unregularized and regularized regression (using cross-validation to determine the optimal LASSO penalty) and build minimal reference-free models with maximal explanatory power. The joint inference of reference-free effects and nonspecific epistasis requires nonlinear regression, which makes this analysis slow. When the number of genotypes exceeds 100,000, cross-validating a second- or third-order model may take days. Furthermore, we find that the optimzer can get stuck in local optima when the fraction of genotypes in the dynamic range of measurement is very small. The next example illustrates an approximate optimization that is faster and more robust yet largely retains the accuracy.

**RFA-fast:** Instead of jointly inferring nonspecific epistasis with the reference-free effects, this analysis asks the user to specify a sigmoid function for modeling nonspecific epistasis.
