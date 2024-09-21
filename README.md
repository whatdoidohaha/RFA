# Reference-free analysis of genetic architecture

Reference-free analysis is a powerful formalism for analyzing the genetic architecture of proteins and nucleic acids (https://www.nature.com/articles/s41467-024-51895-5). Here we provide scripts and tutorials for performing reference-free analysis on experimental data. The directory */Tutorials* provides annotated workflows and example datasets. */Scripts* contains the function implementations, which we recommend accessing through the annotated workflows. *Paper.tar.gz* contains all the data and scripts in our ppaer. 

## Tutorials
We illustrate four workflows.

**RFA:** Here reference-free effects and nonspecific epistasis are jointly inferred. Any form of genotype space, including when the number of states varies among sites, is supported. We perform unregularized and regularized regression (using cross-validation to determine the optimal LASSO penalty) and build minimal reference-free models with maximal explanatory power. The joint inference of reference-free effects and nonspecific epistasis requires nonlinear regression, which makes this analysis slow. When the number of genotypes exceeds 100,000, cross-validating a second- or third-order model may take days. Furthermore, the optimizer can get stuck in local optima when the fraction of genotypes in the dynamic range of measurement is very small. The next workflow illustrates an approximate method that is faster and more robust yet largely retains the accuracy.

**RFA-fast:** Instead of jointly inferring reference-free effects and nonspecific epistasis, this workflow asks the user to specify the nonspecific epistasis parameters. It then performs generalized linear regression using the specified function as link function. The result is the best estimate of reference-free effects under the specified shape of nonspecific epistasis. The highly optimized R package `glmnet` makes this workflow fast and robust. Once the reference-free effects are inferred, the nonspecific epistasis parameters can be optimized while fixing the effects; this two-step procedure can be iterated to further improve the model fit.

**RFA-binary:** Here reference-free effects and nonspecific epistasis are jointly inferred using scripts tailored for binary genotype spaces. When there are only two states per site, the definition of reference-free effects implies that all effects in a site-combination equal in magnitude and vary only in sign. Therefore, only one effect needs to be explicitly modeled and inferred for each site-combination. This compact encoding has a few shortcomings (see Fig. 3B and associated text in our paper), but for binary genotype spaces the shortcomings are minor and outweighed by computational efficiency of the compact encoding.

**RFA-fast-binary:** The fast, approximate inference in **RFA-fast** is tailored for binary genotype spaces.

We use three published datasets for our tutorials. The avGFP dataset comprehensively samples a binary genotype space over 13 sites in a fluorescent protein (total 8,192 genotypes; see the notes for each dataset for citation and more detail). The phenotype is the average fluorescence at two wavelengths. The CR9114-B dataset is a near-complete sample of a binary genotype space over 16 sites in an antibody (99.7% sampling of 65,536 possible genotypes). The phenotype is the affinity towards the influenza strain B hemagglutinin; only 0.1% of genotypes have a phenotype above the lower bound of 6. The ParD3-ParE3 dataset is a near-complete sample of a genotype space over four sites in the bacterial antitoxin ParD3 (98.2% sampling of 9,360 possible genotypes). The number of states varies across sites from 6 to 12, and the phenotype is the absolute fitness conferred by binding to the toxin ParE3.

## Notes on model interpretation

As in any statistical inference, care should be taken in interpreting the inferred reference-free model. We describe three causes of misinterpretation and suggest safeguard practices.

*Overfitting.* In the experimental datasets we have analyzed, third-order models exhibit almost the same out-of-sample predictive accuracy as second-order models, indicating that third-order effects contribute negligibly to phenotype. Surprisingly, however, some third-order models return large third-order effects that account for a substantial amount of variance at the level of genetic score. This occurs because LASSO regularization is not fully sufficient to prevent the overfitting of the very large number of parameters in third-order models. This can be seen in the significantly better in-sample fit than out-of-sample fit in cross-validation.

Several pratices can safeguard against overfitting. First, overfitting should be diagnosed by comparing the in-sample and out-of-sample fit in cross-validation. It is particularly helpful to plot the observed phenotype against prediction for both in-sample and out-of-sample predictions, which we do in our workflows. Second, the phenotypic contribution of an effect or a set of effects should be measured by separately fitting a model with and without them and comparing the model fit, rather than partitioning the variance in the full model. Third, higher-order effects can be modeled for specific site-combinations where interactions are expected rather than for all site-combinations of the given order. Our workflow provides an easy interface to choose a desired subset of site-combinations. Should the need arise, an algorithm that automatically distinguishes important site-combinations may also be devised.

*Inaccurate account of nonspecific epistasis.* The sigmoid link function captures the primary cause of nonspecific epistasis - phenotype bounding - but the sigmoid curvature between the bounds may not fit a particular dataset well. This can be diagnosed by a nonlinearity in the plot of observed versus predicted phenotype, and would cause specific epistatic interactions to be overestimated in both extent and order. We have only implemented the sigmoid link function but have built an easy interface to implement additional link functions.

*Limited dynamic range of measurement.* Consider an amino acid state that is incompatible with function and causes every genotype containing that state to be at the lower bound. We would know that the state has a strongly negative effect, but we would not be able to assign an exact value to its effect because any value that is sufficiently negative would imply the same phenotypes; the value inferred by the model will be highly sensitive to the regularization method. The dynamic range of measurement should therefore be inspected before attempting a quantitative interpretation of the inferred effects. If the range is severely limited, the only way to improve the inference is to repeat the experiment with better resolution.
