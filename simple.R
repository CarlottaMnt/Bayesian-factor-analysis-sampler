library(matlab)
library(Matrix)
library(mvtnorm)
library(invgamma)
library(tidyverse)
library(matrixcalc)


rep.row <- function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

#inits <- list(lambda=matrix(0,J,L),
#              mu=matrix(0,J,1),
#              Sigma=diag(J),
#              PSI=diag(N),
#              delta=matrix(1,nrow=N,ncol=1),
#              a = 0
#              )
#
simple <- function(n_iter = 1000,
                   burn_in = 500,
                   y,
                   intis,
                   g=0,G=10000,
                   alpha=1/1000,beta=1/1000,
                   V_mu,
                   J=NULL,
                   mu_a=NULL,V_a= 10,
                   R=NULL,
                   W=NULL,
                   O=NULL,
                   gamma=NULL,
                   w=NULL,
                   tuning=NULL,
                   model="0",
	    tuning_lambda = NULL){
  cat("Start Gibbs...\n") 
  param <- inits
  V.hat <- rep(NA,J)
  lambda.hat <- rep(NA,J)
  V.mu <- rep(NA,J)
  sd <- tuning
  mu.hat <- rep(NA,J)
  beta.hat <- rep(NA,J)
  prova <- rep(NA,N)
  y_1 <- rep(NA,N*J)
  o_c <- rep(NA,N)
  W_c=matrix(0,N,N)
  W_d_04=matrix(0,N,N)
  samples <- list(lambda=matrix(NA,n_iter-burn_in,J*L),
                  delta = matrix(NA,n_iter-burn_in,N*L),
                  a = rep(NA,n_iter-burn_in),
                  y_rep = matrix(NA,n_iter-burn_in,N*J),
                  d = matrix(NA,n_iter-burn_in,N*L),
                  D = matrix(NA,n_iter-burn_in,(N*L)*(N*L)),
                  PSI = matrix(NA,n_iter-burn_in,(N*L)*(N*L)),
                  Sigma=matrix(NA,n_iter-burn_in,J),
                  mu=matrix(NA,n_iter-burn_in,J)
                  )
  
  start_time = Sys.time()
  count=0
  z_count = 0
  for (t in 2:(n_iter+1))
  {
    cat("\r", "Iteration ", t-1,"time",Sys.time()- 
          start_time, "acc. ", count/(t-1), "delta", quantile(param$delta, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)), "IQR", IQR(param$delta))
    if ((t-1) %% 500 == 0) cat("\n")
    start_time = Sys.time()
    for (j in 1:J)
    {
    #Step 1
      V.hat[j]= 1/(G^(-1) + t(param$delta)%*%param$delta/
                   (param$Sigma[j,j]))
      lambda.hat[j] = V.hat[j]*(g/G + (t(param$delta)%*%
                                         (y[,j]-ones(N)[,1]*param$mu[j])))/
                                  param$Sigma[j,j]
      
      param$lambda[j,L] = rnorm(1,lambda.hat[j],sqrt(V.hat[j]))
      if (t <= burn_in)
      {
        if (j==1)
        {
         param$lambda[j,L] = rnorm(1,lambda.hat[j],sqrt(V.hat[j]))
         z_count =1
         while (param$lambda[j,L]<=0)
         {
           param$lambda[j,L] = rnorm(1,lambda.hat[j],sqrt(V.hat[j]))
           z_count = z_count+1
           if (z_count > 100)
           param$lambda[j,L] = tuning_lambda
        }
       }
      }
     }
     for (j in 1:J)
     {
    #Step 2
      V.mu[j]= 1/(0.001+(N/param$Sigma[j,j]))
      mu.hat[j]=V.mu[j]*t(ones(N)[,j])%*%(y[,j]-param$delta*param$lambda[j,L])/
        param$Sigma[j,j]
      param$mu[j]= rnorm(1,mu.hat[j],sqrt(V.mu[j]))    
      }     
    for (j in 1:J)
    {
    if ( t > burn_in)
    samples[["lambda"]][t-burn_in-1,j]=param$lambda[j,L]
    samples[["mu"]][t-burn_in-1,j]=t(param$mu[j])  
    }
    #Step 3
      y_1=as.vector(t(y))
      D = solve(param$PSI + diag(sum((param$lambda[,1])^2 /
                                      diag(param$Sigma)),N))		        
      for (s in (0:(N-1)))
      {
        prova[s+1]= sum((param$lambda[,1]/diag(param$Sigma)) * 
                          (y_1[(J*s+1):(J *(s+1))]-param$mu))
      }
      
      d = D %*% prova
      param$delta= t(rmvnorm(1,d,D))
      if ( t > burn_in)
      {
      samples[["delta"]][t-burn_in-1,]=t((param$delta))
      samples[["d"]][t-burn_in-1,]=d
      samples[["D"]][t-burn_in-1,]=as.vector(D)
      }
      
      #Step 4
      for (j in 1:J)
      {
        alpha.hat=(alpha + N )/2
        beta.hat[j]=((t(y[,j] - ones(N)[,j]*param$mu[j]- (param$delta)
                       * param$lambda[j,L]) %*%
                       (y[,j] - ones(N)[,j]*param$mu[j]-
                                    param$delta * param$lambda[j,L]))+ beta)/2
        param$Sigma[j,j]=rinvgamma(1,alpha.hat,beta.hat[j])
      }
    #Step 5 Metropolis Hastings
     #Model 2(Marginal)
      if (model != "0")
      {
        if (model == "1")
        {
          u = rnorm(1,0,1)
          proposed_a = param$a + u*sd
          if (proposed_a <=0)
          {
            A = log(0)
          } else 
          {
          target_1 = function(x,mu_x,W)
           {
            PSI=exp(-x* W)
            PSI[lower.tri(PSI)] <- t(PSI)[lower.tri(PSI)]
            return(dmvnorm(t(param$delta),rep(0,N),PSI,log = TRUE)+
                     dnorm(x,mu_a,V_a, log = TRUE) + dnorm(mu_x,x,tuning,log = TRUE))
           }
         A = target_1(proposed_a,param$a,W)- target_1(param$a,proposed_a,W)  
         }
        }
        #Model 3A
        if (model == "2")
        { 
          Xi <- eigen(R)
          xi <- sort(Xi$values)
          xi_inv <- 1/(xi)
          u = rnorm(1,0,1)
          proposed_a = param$a + u*sd
          if (xi_inv[1]> proposed_a  | proposed_a > xi_inv[N])
          {
            A = log(0)#proposed_a = rnorm(1,param$a,sqrt(tuning))
          } else 
          {
          target_2 <- function(x,mu_x,R,O)
          {
           B =O-x*R
           B_inv = solve(B)
           M = ones(N,1)/sqrt(diag(B_inv))
           M_1 = diag(M[,1])
           PSI = M_1 %*% B_inv %*% M_1
           PSI[lower.tri(PSI)] <- t(PSI)[lower.tri(PSI)]
           return(dmvnorm(t(param$delta),rep(0,N),PSI,log = TRUE)+
                     dnorm(x,mu_a,V_a, log = TRUE) + dnorm(mu_x,x,tuning,log = TRUE))
          }   
          A = (target_2(proposed_a,param$a,R,O))- (target_2(param$a,proposed_a,R,O))  
          }            
        }
        #Model 3B
        if (model=="3")
        { 
          Xi <- eigen(R)
          xi <- sort(Xi$values)
          xi_inv <- 1/(xi)
          u = rnorm(1,0,1)
          proposed_a = param$a + u*sd
          if (xi_inv[1] > proposed_a  | proposed_a > xi_inv[N])
          {
            A = log(0) #proposed_a = rnorm(1,param$a,sqrt(tuning))
          } else
          {
          target_3 <- function(x,mu_x,R,W,O)
          {
            B=O-x*W
            B_inv=solve(B)
            M = ones(N,1)/sqrt(diag(B_inv))
            M_1 = diag(M[,1])
            PSI=M_1 %*% B_inv %*% M_1
            PSI[lower.tri(PSI)] <- t(PSI)[lower.tri(PSI)]
            return(dmvnorm(t(param$delta) , rep(0,N) , PSI , log = TRUE ) +
                     dnorm(x,mu_a,V_a, log = TRUE) + dnorm(mu_x,x,sqrt(tuning),log = TRUE))
          }   
          A =  (target_3(proposed_a, param$a,R,W,O))-
            (target_3(param$a,proposed_a,R,W,O))   
          }
        }
        # cat("Proposal", A)
        if(log(runif(1)) < A)
        {
          param$a = proposed_a
          count = count + 1
        }
        if (t <= burn_in)
        {
        sd = sd*sqrt(1+min(1,(t-1)^(-2/3))*(min(1,exp(A))-0.44))
        }
      }
      if (t >  burn_in)
       samples[["a"]][t-burn_in-1]=param$a
       if (model == "1")
        {
          param$PSI=exp(-param$a* W)
          param$PSI[lower.tri(param$PSI)] <- t(param$PSI)[lower.tri(param$PSI)]
        }
        if (model == "2")
        {
          B=O-param$a*R
          diag(B) <- 1
          B_inv=solve(B)
          M=diag(B_inv)
          M_1=sqrt(diag(1/M))
          param$PSI=M_1 %*% B_inv %*% M_1  
          param$PSI[lower.tri(param$PSI)] <- t(param$PSI)[lower.tri(param$PSI)]
        }
        if (model =="3")
        {
          B=O-param$a*W
          B_inv=solve(B)
          M=diag(B_inv)
          M_1=sqrt(diag(1/M))
          param$PSI=M_1 %*% B_inv %*% M_1
          param$PSI[lower.tri(param$PSI)] <- t(param$PSI)[lower.tri(param$PSI)]
        }   
      if (t > burn_in)
      {
         samples[["PSI"]][t-burn_in-1,]=as.vector(param$PSI)  
      } 	  
      
      if ( t > burn_in)
      samples[["Sigma"]][t-burn_in-1,]=as.vector(diag(param$Sigma))
      for (s in (0:(N-1)))
      {
       samples[["y_rep"]][t-burn_in-1,(J*s+1):(J *(s+1))]= param$mu +
       param$lambda[,1] * param$delta[(s+1)]
      }
  }
  cat("\nDONE\n","Metropolis Hastings Acceptance rate", count/n_iter )
  return(samples)
 }



