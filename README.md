# Bayesian latent factor analysis for spatially correlated data: Sampler.

The fuction implements the models specified in Hogan and Tchernis (2004).
In particular it builds five Metropolis Hastings steps, one for each model, within a Gibbs sampling alghorithm. 

The empirical illustration is based on the public data-set from "A misura di comune". In particular I focus on the Varese province which entails 138 different areas. I analyse the year 2014 and 2015. The output of the function are simulated posterior parameters distributions, in particular I return the the estimates for <img src="https://render.githubusercontent.com/render/math?math=\boldsymbol{\lambda}">, <img src="https://render.githubusercontent.com/render/math?math=\delta_i">,<img src="https://render.githubusercontent.com/render/math?math= \boldsymbol{Y}^{\textrm{rep}}"> and **d**. 

The main purpose of the function is to model the Varese areas' spatial correlation, allowing the inclusion of such information in the latent factor estimation. Within the implementation, the latent factor should represents the latent well-being of each areas. 


## References:

J. W. Hogan and R. Tchernis.  Bayesian factor analysis for spatially correlated data, with applica-tion to summarizing area-level material deprivation from census data.Journal of the AmericanStatistical Association, 99(466):314–324, 2004