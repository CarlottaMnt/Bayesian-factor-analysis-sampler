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
                   mu_a=NULL,V_a=1000,
                   R=NULL,
                   W=NULL,
                   O=NULL,
                   gamma=NULL,
                   w=NULL,
                   tuning=NULL,
                   model="0"){
  cat("Start Gibbs...\n") 
  param <- inits
  V.hat <- rep(NA,J)
  lambda.hat <- rep(NA,J)
  V.mu <- rep(NA,J)
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
                  #D = matrix(NA,n_iter,(N*L)*(N*L)),
                  #Diff=matrix(NA,n_iter,(N*L)*(N*L)),
                  PSI = matrix(NA,n_iter-burn_in,(N*L)*(N*L)),
                  Sigma=matrix(NA,n_iter-burn_in,J),
                  mu=matrix(NA,n_iter-burn_in,J)
                  )
  
  start_time = Sys.time()
  count=0
  for (t in 2:(n_iter+1))
  {
    cat("\r", "Iteration ", t-1,". Last iteration took",Sys.time()- 
          start_time, "\n")
    start_time = Sys.time()
    for (j in 1:J)
    {
      #Step 1
      V.hat[j]= (1/G + ((t(param$delta)%*%param$delta))/
                   param$Sigma[j,j])^(-1)
      lambda.hat[j] = V.hat[j]*(g/G + (t(param$delta)%*%
                                         (y[,j]-ones(N)[,1]*param$mu[j]))/
                                  param$Sigma[j,j])
      param$lambda[j,L] = rnorm(1,lambda.hat[j],sqrt(V.hat[j]))
      cat("j",j,"mean",lambda.hat[j],"variance",V.hat[j])
      if (j==1)
      {
        while (param$lambda[j,L]<=0)
        {
          param$lambda[j,L] = rnorm(1,lambda.hat[j],sqrt(V.hat[j]))
        }
      cat("Exit")
      } 
      cat("j",j,"lambda",param$lambda[j,L])
      #Step 2
      V.mu[j]= param$Sigma[j,j]*V_mu/(param$Sigma[j,j]+V_mu*N)
      mu.hat[j]=V.mu[j]*t((ones(N)[,j]))%*%(y[,j]-param$delta*param$lambda[j,L])/
        param$Sigma[j,j]
      param$mu[j]= rnorm(1,mu.hat[j],sqrt(V.mu[j]))
      cat("j",j,"mu",param$mu[j])
      if ( t > burn_in)
      samples[["lambda"]][t-burn_in-1,j]=param$lambda[j,L]
      samples[["mu"]][t-burn_in-1,]=t(param$mu)
    }
    
      #Step 3
      cat("Start step 3")
      #Lambda = kronecker(diag(N),param$lambda)
      #SSigma=  solve(kronecker(diag(N),param$Sigma))
      y_1=as.vector(t(y))
      D = solve(param$PSI + diag(sum((param$lambda[,1])^2 /
                                      diag(param$Sigma)),N))
      #D = solve(param$PSI+t(Lambda)%*%(SSigma)%*%(Lambda))
      for (s in (0:(N-1)))
      {
        prova[s+1]= sum((param$lambda[,1]/diag(param$Sigma)) * 
                          (y_1[(J*s+1):(J *(s+1))]-param$mu))
      }
      #cat("prova",prova)
      d = D %*% prova
      #d = D %*% t(Lambda)%*% (SSigma)%*%(as.vector(y-rep.row(param$mu,N)))
      cat("End step 3")
      param$delta= t(rmvnorm(1,d,D))
      cat("delta",param$delta)
      if ( t > burn_in)
      {
      samples[["delta"]][t-burn_in-1,]=t((param$delta))
      samples[["d"]][t-burn_in-1,]=d
      }
      #samples[["D"]][t-1,]=as.vector(D)
      #Step 4
      for (j in 1:J)
      {
        alpha.hat=(alpha + N )/2
        beta.hat[j]=(t(y[,j] - ones(N)[,j]* param$mu[j]- (param$delta)
                       * param$lambda[j,L])%*%
                       (y[,j] - ones(N)[,j]*param$mu[j]-
                                    param$delta * param$lambda[j,L])+ beta)/2
        cat("j",j,"beta",beta.hat[j])
        param$Sigma[j,j]=rinvgamma(1,alpha.hat,beta.hat[j])
        cat("Sigma", param$Sigma[j,j])
      }
      if ( t > burn_in)
      samples[["Sigma"]][t-burn_in-1,]=as.vector(diag(param$Sigma))
      #Step 5 Metropolis Hastings
      #Model 2(Marginal)
      if (model != "0")
      {
        if (model == "1")
        {
          proposed_a = rnorm(1,param$a,sqrt(tuning))
          cat("a_new",proposed_a,"a_old",param$a)
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
                     dnorm(x,mu_a,V_a, log = TRUE) + dnorm(mu_x,x,sqrt(tuning),log = TRUE))
           }
         A = target_1(proposed_a,param$a,W)- target_1(param$a,proposed_a,W)  #+ pnorm(proposed_a, 0,sqrt(tuning),log.p = TRUE))-
             #+ pnorm(param$a,0,sqrt(tuning),log.p = TRUE))
         }
        }
        #Model 3A
        if (model == "2")
        { 
          Xi <- eigen(R)
          xi <- sort(Xi$values)
          xi_inv <- 1/(xi)
          proposed_a = rnorm(1,param$a,tuning)
          if (xi_inv[1] > proposed_a  | proposed_a > xi_inv[N])
          {
            A = log(0)#proposed_a = rnorm(1,param$a,sqrt(tuning))
          } else 
          {
          target_2 <- function(x,mu_x,R,O)
          {
            B = O-x*R
            B_inv = solve(B)
            M = diag(B_inv)
            M_1 = sqrt(diag(1/M))
            PSI = M_1 %*% B_inv %*% M_1
            #PSI[lower.tri(PSI)] <- t(PSI)[lower.tri(PSI)]
           return(dmvnorm(t(param$delta),rep(0,N),PSI,log = TRUE)+
                     dnorm(x,mu_a,V_a, log = TRUE) + dnorm(mu_x,x,tuning,log = TRUE))
          }   
          A = (target_2(proposed_a,param$a,R,O))- #+ log(pnorm(xi_inv[N],param$a,sqrt(tuning)) - pnorm(xi_inv[1],param$a,sqrt(tuning)))) -
            (target_2(param$a,proposed_a,R,O))# + log(pnorm(xi_inv[N],proposed_a,sqrt(tuning)) - pnorm(xi_inv[1],proposed_a,sqrt(tuning))))
          }            
        }
        #Model 3B
        if (model=="3")
        { 
          Xi <- eigen(R)
          xi <- sort(Xi$values)
          xi_inv <- 1/(xi)
          proposed_a = rnorm(1,param$a,tuning)
          if (xi_inv[1] > proposed_a  | proposed_a > xi_inv[N])
          {
            A = log(0) #proposed_a = rnorm(1,param$a,sqrt(tuning))
          } else
          {
          target_3 <- function(x,mu_x,R,W,O)
          {
            B=O-x*W
            diag(B) <- diag(O)
            B_inv=solve(B)
            M=diag(B_inv)
            M_1=sqrt(diag(1/M))
            PSI=M_1 %*% B_inv %*% M_1
            return(dmvnorm( t(param$delta) , rep(0,N) , PSI , log = TRUE ) +
                     dnorm(x,mu_a,V_a, log = TRUE) + dnorm(mu_x,x,tuning,log = TRUE))
          }   
          A =  target_3(proposed_a, param$a,R,W,O)-# # log(pnorm(xi_inv[N],param$a,sqrt(tuning))- pnorm(xi_inv[1],param$a,sqrt(tuning))))-
            target_3(param$a,proposed_a,R,W,O)#+log(pnorm(xi_inv[N],proposed_a,sqrt(tuning))- pnorm(xi_inv[1],proposed_a,sqrt(tuning))))   
          }
        }
        #Model 3C
        if (model =="4")
        {
          proposed_a = rnorm(1,param$a,tuning)
          
          if (proposed_a <=0)
          {
            A = log(0) #proposed_a = rnorm(1,param$a,sqrt(tuning))
          } else
          {
          target_4 <- function(x,mu_x,W)
          {
            for (i in 1:N){
              o_c[i]=sum(exp(-x*W[i,-i]))
            }
            O_c=diag(o_c)
            for (i in 1:N)
            {
              W_c[i,-i]= -exp(-x*W[i,-i])
            }
            
            Diff= O_c-W_c
            #diag(Diff) <- 1
            PSI=solve(Diff)
            return(dmvnorm(t(param$delta),rep(0,N),PSI,log = TRUE) +
                     dnorm(x,mu_a,V_a, log = TRUE) + dnorm(mu_x,x,sqrt(tuning),log = TRUE))
          }
          A =  (target_4(proposed_a, param$a,W))-
            (target_4(param$a,proposed_a,W))
          if (A == "NaN")
          { 
           A= log(0)
          }
          }
        }
        #Model 3D
        if (model %in% c("5","6"))
        { 
          Xi <- eigen(R)
          xi <- sort(Xi$values)
          xi_inv <- 1/(xi)
          proposed_a = rnorm(1,param$a,sqrt(tuning))
          if (xi_inv[1] > proposed_a  | proposed_a > xi_inv[N])
          {
            A = log(0) #proposed_a = rnorm(1,param$a,sqrt(tuning))
          } else 
          {
          target_5 <- function(x,mu_x,R,W,w,gamma)
           { 
            for (i in 1:N)
            {
              W_c[i,w[[i]]]= exp(-W[[i]] * gamma)
            }
            B=diag(N)- x * W_c
            #diag(B) <- 1
            B_inv=solve(B)
            A=diag(B_inv)
            A=sqrt(diag(1/A))
            PSI = A %*% B_inv %*% A
            return(dmvnorm(t(param$delta),rep(0,N),PSI,log = TRUE)+
                     dnorm(x,mu_a,V_a, log = TRUE) + dnorm(mu_x,x,sqrt(tuning),log = TRUE))
          }
          
          A =  (target_5(proposed_a,param$a,W,R=R,gamma=1,w))-#+ pnorm(param$a, proposed_a,sqrt(tuning),log.p = TRUE))-
            (target_5(param$a,proposed_a,W,R=R,gamma=1,w))#+pnorm(proposed$a,param$a,sqrt(tuning),log.p = TRUE))
         }
        }
        cat("Proposal", A)
        if(log(runif(1)) < A)
        {
          param$a = proposed_a
          count = count + 1
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
        }
        if (model =="3")
        {
          B=O-param$a*W
          #diag(B) <- diag(O)
          B_inv=solve(B)
          M=diag(B_inv)
          M_1=sqrt(diag(1/M))
          param$PSI=M_1 %*% B_inv %*% M_1
        }
        if (model =="4")
        {
          for (i in 1:N){
            o_c[i]=sum(exp(-param$a*W[i,]))
          }
          O_c = diag(o_c)
          for (i in 1:N[1])
          {
            W_c[i,]= -exp(-param$a*W[i,])
          }
          
          Diff= O_c-W_c
          #diag(Diff) <- 1
          #samples[["Diff"]][t-1,]=Diff
          param$PSI=solve(Diff)
        }
        if (model %in%  c("5","6"))
        {
          for (i in 1:N)
          {
            W_d_04[i,w[[i]]]= -exp(-gamma*W[[i]])
          }
          B=diag(N)- param$a * W_d_04
          diag(B) <-1
          #samples[["Diff"]][t-1,]=B
          B_inv=solve(B)
          M=diag(B_inv)
          M=sqrt(diag(1/M))
          param$PSI = M %*% B_inv %*% M 
        }
      
      if (t > burn_in)
	  {
	  samples[["PSI"]][t-burn_in-1,]=as.vector(param$PSI) 
	  for (s in (0:(N-1)))
      {
       samples[["y_rep"]][t-burn_in-1,(J*s+1):(J *(s+1))]= param$mu +
       param$lambda[,1] * param$delta[(s+1)]
      }
	  }
  }
  cat("\nDONE\n","Metropolis Hastings Acceptance rate", count/(n_iter))
  return(samples)
}


