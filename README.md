# Reference-free analysis of genetic architecture

Here we provide 

Reference-free analysis of genetic architecture


## Tutorials
We provide tutorials for four types of analyses.

**RFA:** Here, reference-free effects on genetic score and nonspecific epistasis are jointly inferred. Any form of genotype space, including when the number of states varies among sites, is supported. We perform unregularized and regularized regression (using cross-validation to determine the optimal LASSO penalty) and build minimal reference-free models with maximal explanatory power. The joint inference of reference-free effects and nonspecific epistasis requires nonlinear regression, which makes this analysis slow. When the number of genotypes exceeds 100,000, cross-validating a second- or third-order model may take days. Furthermore, we find that the optimzer can get stuck in local optima when the fraction of genotypes in the dynamic range of measurement is very small. The next example illustrates an approximate optimization that is faster and more robust yet largely retains the accuracy.

**RFA-fast:** Instead of jointly inferring nonspecific epistasis with reference-free effects, this analysis asks the user to specify a sigmoid function for modeling nonspecific epistasis and then performs generalized linear regression using the sigmoid as link function. The result is the best estimate of reference-free effects under the fixed shape of nonspecific epistasis. The highly optimized R package `glmnet` is used to make this procedure fast and robust. Once the reference-free effects are inferred, the sigmoid can be optimized while fixing the effects; this two-step procedure can be iterated to further improve the model fit.

**RFA-binary:** Here, joint inference of reference-free effects and nonspecific epistasis is performed using scripts tailored for binary genotype spaces. When there are only two states per site, the zero-mean property of reference-free effects causes all effects in a site-combination to equal in magnitude and vary only in sign. Therefore, only one effect needs to be explicitly modeled and inferred for each site-combination. (Simplex encoding and graph Fourier transform exploit the zero-mean property for multiple states to always achieve this compact encoding.) This compact encoding makes the model difficult to intuit; it also causes the phenotype of any genotype to depend on every model parameter and therefore the inference to be more sensitive to missing genotypes, as detailed in Fig. 1E of our paper. However, the impact on sensitivity is minor when the number of states is small and is outweighed by the speed and memory benefit of compact encoding.

**RFA-binary-fast:** The fast, approximate inference in **RFA-fast** is performed specifically for binary genotype spaces.

We provide three published datasets for tutorials. The **avGFP** dataset comprehensively samples a binary genotype space over 13 sites in a fluorescent protein (see the notes in each dataset directory for citation and more detail). The phenotype is the average fluorescence at two wavelengths. The **CR9114-B** dataset is a near-complete sample of a binary genotype space over 16 sites in an antibody. The phenotype is the affinity towards the influenza strain B hemagglutinin. Only 0.1% of genotypes lie above the lower phenotype bound of 6. The **ParD3-ParE3** dataset is a near-complete sample of a genotype space over four sites in the bacterial antitoxin ParD3. The number of states varies across sites from 6 to 12. The phenotype is the absolute fitness conferred by binding to the toxin ParE3.

## Notes on model interpretation




