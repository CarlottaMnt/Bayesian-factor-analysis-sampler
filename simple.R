#install.packages("matlab")
#install.packages("Matrix")
#install.packages("mvtnorm")
#install.packages("invgamma")
#install.packages("matrixcalc")
library(matlab)
library(Matrix)
library(mvtnorm)
library(invgamma)
library(tidyverse)
library(matrixcalc)

#Function
simple <- function(n_iter = 1000,
                   y,
                   intis,
                   g=0,G=1000,
                   alpha=1/1000,beta=1/1000,
                   V_mu,
                   mu_a=0,V_a=1000,
                   R=NULL,
                   W=NULL,
                   O=NULL,
                   gamma=40,
                   w=NULL,
                   tuning=1,
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
  samples <- list(lambda=matrix(NA,n_iter,J*L),
                  delta=matrix(NA,n_iter,N*L),
                  a=rep(NA,n_iter),
                  y_rep = matrix(NA,n_iter,N*J),
                  d=matrix(NA,n_iter,N*L),
                  D=matrix(NA,n_iter,(N*L)*(N*L)),
                  Diff=matrix(NA,n_iter,(N*L)*(N*L)),
                  PSI=matrix(NA,n_iter,(N*L)*(N*L)),
                  Sigma=matrix(NA,n_iter,J*J)
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
                                         (y[,j]-ones(N)[,j]*param$mu[j]))/
                                  param$Sigma[j,j])
      param$lambda[j,L] = rnorm(1,lambda.hat[j],sqrt(V.hat[j]))
      if (j==1)
      {
        while (param$lambda[j,L]<=0)
        {
          param$lambda[j,L] = rnorm(1,lambda.hat[j],sqrt(V.hat[j]))
        }
      }
      
      #Step 2
      V.mu[j]= param$Sigma[j,j]*V_mu/(param$Sigma[j,j]+V_mu*N)
      mu.hat[j]=V.mu[j]*t(ones(N)[,j])%*%(y[,j]-param$delta*param$lambda[j,L])/
        param$Sigma[j,j]
      param$mu[j]= rnorm(1,mu.hat[j],sqrt(V.mu[j]))
      samples[["lambda"]][t-1,j]=param$lambda[j,L]
    }
    
      #Step 3
      cat("Start step 3")
      y_1=as.vector(t(y))
      D= solve(param$PSI + diag(sum((param$lambda[,1])^2 /
                                      diag(param$Sigma)),N))
      #D = solve(param$PSI+t(Lambda)%*%(SSigma)%*%(Lambda))
      for (j in (0:(N-1)))
      {
        prova[j+1]= sum((param$lambda[,1]/diag(param$Sigma)) * 
                          (y_1[(J*j+1):(J *(j+1))]-param$mu))
      }
      d = D %*% prova
      cat("End step 3")
      param$delta= t(rmvnorm(1,d,D))
      samples[["delta"]][t-1,]=t((param$delta))
      samples[["d"]][t-1,]=d
      samples[["D"]][t-1,]=as.vector(D)
      #Step 4
      for (j in 1:J)
      {
        alpha.hat=(alpha + N )/2
        beta.hat[j]=(t(y[,j] - ones(N)[,j]* param$mu[j]- (param$delta)
                       * param$lambda[j,L])%*%
                       ((y[,j] - ones(N)[,j]*param$mu[j]-
                                    param$delta * param$lambda[j,L]))+ beta)/2
        param$Sigma[j,j]=rinvgamma(1,alpha.hat,beta.hat)
      }
      samples[["Sigma"]][t-1,]=as.vector(diag(param$Sigma))
      #Step 5 Metropolis Hastings
      if (model!="0")
      {
        proposed_a = rnorm(1,param$a,sqrt(tuning))
      #Model 2
      
        if (model == "1")
        {
          target_1 = function(x,mu_x,W)
            {
              if (x <= 0)
              {
                print("Proposal less equal 0")
                return(0)
              }
              PSI=exp(-x* W)
              PSI[lower.tri(PSI)] <- t(PSI)[lower.tri(PSI)]
              return(dmvnorm(t(param$delta),rep(0,N),PSI)*
                       dnorm(x,mu_a,V_a)*dnorm(mu_x,x,sqrt(tuning)))
            }
          A =  (target_1(proposed_a,param$a,W=PSI_04))/
            (target_1(param$a,proposed_a,W=PSI_04))
          if (A == "NaN")
          {
            cat("A",A)
            A=param$a
          }
            
        }
      #Model 3A
      
        if (model== "2")
        { 
          target_2 <- function(x,mu_x,R,O)
            {
              Xi <- eigen(R)
              xi <- sort(Xi$values)
              xi_inv <- 1/(xi)
              if (xi_inv[1] > x  | x > xi_inv[N])
              {
                print("Proposal out of support")
                return(0)
              }
              B=O-x*R
              B_inv=solve(B)
              M=diag(B_inv)
              M_1=sqrt(diag(1/M))
              PSI=M_1 %*% B_inv %*% M_1
              #PSI[lower.tri(PSI)] <- t(PSI)[lower.tri(PSI)]
              return(dmvnorm(t(param$delta),rep(0,N),PSI)*
                       dnorm(x,mu_a,V_a)* dnorm(mu_x,x,sqrt(tuning)))
            }   
          A =  (target_2(proposed_a,param$a,R,O))/
            (target_2(param$a,proposed_a,R,O))
        }
      #Model 3B
        if (model=="3")
        { 
          target_3 <- function(x,mu_x,R,W,O)
          {
            Xi <- eigen(R)
            xi <- sort(Xi$values)
            xi_inv <- 1/(xi)
            if (xi_inv[1] > x  | x > xi_inv[N])
            {
              print("Proposal out of support")
              return(0)
            }
            B=O-x*W
            diag(B) <- diag(O)
            B_inv=solve(B)
            M=diag(B_inv)
            M_1=sqrt(diag(1/M))
            PSI=M_1 %*% B_inv %*% M_1
            return(dmvnorm(t(param$delta),rep(0,N),PSI)*
                     dnorm(x,mu_a,V_a)* dnorm(mu_x,x,sqrt(tuning)))
          }   
          A =  (target_3(proposed_a, param$a,R,W,O))/
            (target_3(param$a,proposed_a,R,W,O))
        }
      #Model 3C
        if (model =="4")
        {
          target_4 <- function(x,mu_x,W)
            {
              if (x <= 0)
              {
                print("Proposal less equal 0")
                return(0)
              }
              for (i in 1:N){
                o_c[i]=sum(exp(-x*W[i,]))
              }
              O_c=diag(o_c)
              for (i in 1:dim(y)[1])
              {
                W_c[i,]= exp(-x*W[i,])
              }
              
              Diff= O_c-W_c
              diag(Diff) <- 1
              PSI=solve(Diff)
              return(dmvnorm(t(param$delta),rep(0,N),PSI)*
                       dnorm(x,mu_a,V_a)* dnorm(mu_x,x,sqrt(tuning)))
            }
          A =  (target_4(proposed_a, param$a,W))/
            (target_4(param$a,proposed_a,W))
          if (A == "NaN")
          {
            cat("A",A)
            A=param$a
          }
        }
        #Model 3D
        if (model=="5")
        {
          target_5 <- function(x,mu_x,R,W,w,gamma)
          { 
            Xi <- eigen(R)
            xi <- sort(Xi$values)
            xi_inv <- 1/(xi)
            if (xi_inv[1] > x  | x > xi_inv[N])
            {
              print("Proposal out of support")
              return(0)
            }
            for (i in 1:N)
            {
              W_c[i,w[[i]]]= exp(-W[[i]] * gamma)
            }
            B=diag(N)- x * W_c
            diag(B) <- 1
            B_inv=solve(B)
            A=diag(B_inv)
            A=sqrt(diag(1/A))
            PSI = A %*% B_inv %*% A
            return(dmvnorm(t(param$delta),rep(0,N),PSI)*
                     dnorm(x,mu_a,V_a)* dnorm(mu_x,x,sqrt(tuning)))
          }
          A =  (target_5(proposed_a,param$a,W,R=R,gamma=1,w))/
            (target_5(param$a,proposed_a,W,R=R,gamma=1,w))
          if (A == "NaN")
          {
            cat("A",A)
            A=param$a
          }
        }
        #Model 3E
        if(model=="6")
        {
          target_6 <- function(x,mu_x,W,gamma,w)
          {
            if (-1 > x  | x > 1)
            {
              print("Proposal out of support")
              return(0)
            }
            for (i in 1:N){
              o_c[i]=sum(exp(-gamma*W[[i]]))
            }
            O_c=diag(o_c)
            for (i in 1:N)
            {
              W_c[i,w[[i]]]=exp(-gamma*W[[i]])
            }
            B = O_c - x*W_c
            B_inv = solve(B)
            A=diag(B_inv)
            A=sqrt(diag(1/A))
            param$PSI = A %*% B_inv %*% A
            return(dmvnorm(t(param$delta),rep(0,N),param$PSI)*
                     dnorm(x,mu_a,V_a)* dnorm(mu_x,x,sqrt(tuning)))
          }
          A =  (target_6(proposed_a, param$a,W,gamma,w))/
            (target_6(param$a,proposed_a,W))  
        }
        cat("Proposal",A)
        if(runif(1)<= A)
        {
          param$a = proposed_a
          count = count + 1
        }
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
          diag(B) <- diag(O)
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
          O_c=diag(o_c)
          for (i in 1:N[1])
          {
            W_c[i,]= -exp(-param$a*W[i,])
          }
          
          Diff= O_c-W_c
          diag(Diff) <- 1
          samples[["Diff"]][t-1,]=Diff
          param$PSI=solve(Diff)
        }
        if (model =="5")
        {
          for (i in 1:N)
          {
            W_d_04[i,w[[i]]]= -exp(-gamma*W[[i]])
          }
          B=diag(N)- param$a * W_d_04
          diag(B)<-1
          samples[["Diff"]][t-1,]=B
          B_inv=solve(B)
          M=diag(B_inv)
          M=sqrt(diag(1/M))
          param$PSI = M %*% B_inv %*% M 
        }
      samples[["a"]][t-1]=param$a
      samples[["PSI"]][t-1,]=as.vector(param$PSI)
      }
      for (j in (0:(N-1)))
      {
        samples[["y_rep"]][t-1,(J*j+1):(J *(j+1))]= param$mu +
          param$lambda * param$delta[(j+1)]
      }
  }
  cat("\nDONE\n","Metropolis Hastings Acceptance rate", count/n_iter )
  return(samples)
}

