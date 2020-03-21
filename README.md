# Bayesian latent factor analysis for spatially correlated data: Sampler.

The fuction implements the Gibbs sampling alghorithm described in Hogan and Tchernis (2004).
In particular it builds five Metropolis Hastings steps, one for each model's spatial specification highlighted in the paper. The main purpose of the function is to model spatial correlation among geographic areas, including  such information when estimating the model's latent factor. The output of the function are the simulated posterior parameters distributions of <img src="https://render.githubusercontent.com/render/math?math=\boldsymbol{\lambda}">, the vector of factor loadings,  <img src="https://render.githubusercontent.com/render/math?math=\delta_i">, the latent variable in each areas, **d** the expected value of the latent factor and **Y**, posterior replications of the data at hands. 

The empirical illustration is based on the public data-set from "A misura di comune". In particular I focus on the Varese province which entails 138 different areas. I analyse the year 2014 and 2015.  Within this implementation, the latent factor should represents the latent well-being of each areas.

## References:

J. W. Hogan and R. Tchernis.  Bayesian factor analysis for spatially correlated data, with applica-tion to summarizing area-level material deprivation from census data.Journal of the AmericanStatistical Association, 99(466):314–324, 2004