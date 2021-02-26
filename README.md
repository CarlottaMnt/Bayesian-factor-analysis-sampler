# Bayesian latent factor analysis for spatially correlated multivariate data: Sampler.

The fuction implements the Gibbs sampling alghorithm described in Hogan and Tchernis (2004).
In particular it builds five Metropolis Hastings steps, one for each models' spatial specification described in the paper. The main purpose of the function is to model spatial correlation among geographical areas, thus including such information in the models'estimate. The output of the function are the simulated posterior parameters distributions of <img src="https://render.githubusercontent.com/render/math?math=\boldsymbol{\lambda}">, the vector of factor loadings,  <img src="https://render.githubusercontent.com/render/math?math=\delta_i">, the latent variable in each areas, **d** the expected value of the latent factor and **Y**, posterior replications of the data at hands. 

The first empirical illustration is based on the public data-set from "A misura di comune". In particular I focus on the Varese province which entails 138 different areas. I analyse the year 2014 and 2015.  Within this implementation, the latent factor should represents the latent well-being of each areas.

The second empirical illustration relies instead on the Public "Province-BES" dataset. Thus, we extend the analysis to the entire Italian territory assessing the well-being of provinces over the years 2014-2017. 

The main function's input are:

1. **n iter**: number of simulated values.    
1. **burn in** : number of burn in iterations.     
1. **y**: the data frame <img src="https://render.githubusercontent.com/render/math?math={N\times D}"> of scaled elementary indicators.      
1. **inits**: list of initial values for each parameter.       
1. **g, G**, <img src="https://render.githubusercontent.com/render/math?math={\alpha}">, <img src="https://render.githubusercontent.com/render/math?math={\beta}">, <img src="https://render.githubusercontent.com/render/math?math={V_\mu}">: set of prior parameters.       
1. **R, W, w, O**: the adjacency(weight), the distance, the neighbours matrices respectively. Where w is the list indicating neighbouring areas       
1. <img src="https://render.githubusercontent.com/render/math?math={\mu_a}">, <img src="https://render.githubusercontent.com/render/math?math={V_a}">: spatial parameter prior   distribution mean and variance.        
1. **tuning**: the value for the tuning parameter in the Metropolis Hastings step.        
1. **model**: the model to be implemented. Could takes values 0 (Spatial Independence),1 (Marginal specification), 2 (CAR 3A), 3 (CAR 3B). 


## References:

J. W. Hogan and R. Tchernis.  Bayesian factor analysis for spatially correlated data, with applica-tion to summarizing area-level material deprivation from census data.Journal of the American Statistical Association, 99(466):314–324, 2004
