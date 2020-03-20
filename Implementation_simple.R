#Import packages
library(tidyverse)
library(sf)
library(bestNormalize)
library(spdep)
#Import function 
source("simple.R")
#Import data
data <- read.csv("Imputed_dataset.csv")

data_2014 <- data %>%
  as.data.frame() %>% 
  filter(Anno == "Anno_2014", comune!="OSMATE")
data_2015 <- data %>% 
  as.data.frame() %>% 
  filter(Anno =="Anno_2015", comune!="OSMATE")
longlat_2014 <- data_2014 %>% 
  dplyr::select(long,lat) %>% as.matrix()
longlat_2015 <- data_2015 %>% 
  dplyr::select(long,lat) %>% as.matrix()
#Maps---------------------------------------------------------------------------
lombgeo <-  sf::st_read("LOMBARDIA/Comuni_2020_poligonali.shp")
vargeo <- lombgeo[lombgeo$NOME_PRO=="VARESE",]

vargeo$NOME_COM <- as.character(vargeo$NOME_COM)

vargeo <- vargeo %>% 
  mutate(NOME_COM=case_when(
    grepl("BRISSAGO - VALTRAVAGLIA",NOME_COM)~ "BRISSAGO",
    grepl("CADEGLIANO - VICONAGO",NOME_COM)~"CADEGLIANO",
    grepl("CADREZZATE CON OSMATE" ,NOME_COM)~ "CADREZZATE",
    grepl("COCQUIO - TREVISAGO",NOME_COM)~"COCQUIO",
    grepl("CUGLIATE - FABIASCO",NOME_COM)~"CUGLIATE", 
    grepl("LAVENO - MOMBELLO",NOME_COM)~"LAVENO",
    grepl("TRAVEDONA - MONATE",NOME_COM)~ "TRAVEDONA",
    grepl("VIGGIU`",NOME_COM)~"VIGGI?",
    TRUE ~ as.character(NOME_COM)
  ))
vargeo <- vargeo %>% 
  dplyr:: select(-c("CLASSREF","BELFIORE","SIG_PRO","NOME_REG","ANNO","ISTAT","LEGGE_ISTI",
                    "ATS_COD", "ATS_DEC","CMETR_COD","CMETR_DEC")) %>% 
  mutate(comune=NOME_COM)

#Create outcome variables-------------------------------------------------------
#Normalize all the dataset
data_norm_2014 <- apply(data_2014[,3:46],2,bestNormalize) 
data_norm_2015 <- apply(data_2015[,3:46],2,bestNormalize) 

#2014
y_normalize_2014 <- matrix(NA,138,44)
for (i in 1:44)
{
  y_normalize_2014[,i]=data_norm_2014[[i]]$x.t
}

colnames(y_normalize_2014) <- colnames(data_2014[,3:46])

#2015
y_normalize_2015 <- matrix(NA,138,44)
for (i in 1:44)
{
  y_normalize_2015[,i]=data_norm_2015[[i]]$x.t
}
colnames(y_normalize_2015) <- colnames(data_2015[,3:46])

J=dim(y_normalize_2014)[2]
N=dim(y_normalize_2014)[1]
L=1
#Spatial Matrices--------------------------------------------------------
#Model 1: marginal spatial correlation
PSI_2014 <- as.matrix(dist(longlat_2014))#Euclidean Distance
PSI_2015 <- as.matrix(dist(longlat_2015))

diag(PSI_2014) <- 0
diag(PSI_2015) <- 0

#Create matrix of spatial correlation for model 3A
w <- poly2nb(vargeo,row.names=as.character(vargeo$NOME_COM),queen=TRUE)#List of neighbouring areas
R <- nb2mat(w, style='B',zero.policy=TRUE) #Adjagency matrix

#Create the matrix of spatial correlation for model 3B: O and W
o_b<- rep(NA,N)
for (i in 1:N){
  o_b[i]=length(w[[i]])
}
O_b=diag(o_b)

W_b=matrix(0,N,N)
for (i in 1:N)
{
  for (j in 1:N)
  {
    if (j %in% w[[i]])
      W_b[i,j] = sqrt(o_b[i]*o_b[j])
    else
      W_b[i,j]= 0
  }
}
#Create the matrices of spatial correlation for model 3C. 
#Gamma parameter is fixed and can take any values

o_c <- rep(NA,N)
gamma <- 1
for (i in 1:N){
  o_c[i]=sum(exp(-gamma*PSI_2014[i,]))
}
O_c=diag(o_c)
W_c=matrix(0,N,N)
for (i in 1:N)
{
  W_c[i,]=-exp(-gamma*PSI_2014[i,])
}
#Create the matrices of spatial correlation for Model 3D and 3E
PSI_2014_d <- nbdists(w, coords= longlat_2014,longlat = TRUE) %>% as.list()#Create the list of neighbors

#Function implementation
J=dim(y)[2]
N=dim(y)[1]
L=1
rep.row <- function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
inits <- list(lambda=matrix(0,J,L),
              mu=matrix(0,J,1),
              Sigma=diag(J),
              PSI=diag(N),
              delta=matrix(0,nrow=N,ncol=1),
              a = 0
)
#Model 0: No spatial correlation
simple_0_2014 <- simple(1000,y=y_normalize_2014, inits,0,10000,1/1000,1/1000,V_mu=1000,model="0")
simple_0_15 <- simple(1000,y=y_normalize_2015, inits,0,10000,1/1000,1/1000,V_mu=1000,model="0")

# Model 1: Marginal parametrization.
#For this model the initial value for "a" MUST be equal to 1. 
simple_1_14 <- simple(n_iter=100,
                      y=y_normalize_2014,
                      inits,
                      0,10000,
                      1/1000,1/1000,
                      W=PSI_2014,
                      mu_a=0,
                      V_mu=1000,
                      tuning=0.02,
                      model="1")

simple_1_15 <- simple(n_iter=100,
                      y=y_normalize_2015,
                      inits,
                      0,10000,
                      1/1000,1/1000,
                      W=PSI_2015,
                      mu_a=0,
                      V_mu=1000,
                      tuning=0.02,
                      model="1")
#Model 2: CAR model 3A.
simple_2_14 <- simple(100,
                      y=y_normalize_2014,
                      inits,
                      0,10000,
                      1/1000,1/1000,
                      R=R,O=diag(N),
                      V_mu=1000,
                      tuning=0.02,
                      model="2")

simple_2_15 <- simple(1000,
                      y=y_normalize_2015,
                      inits,
                      0,10000,
                      1/1000,1/1000,
                      R=R,O=diag(N),
                      V_mu=1000,
                      tuning=0.002,
                      model="2")
#Model 3: CAR Model 3B
simple_3_14 <- simple(1000,
                      y=y_normalize_2014,
                      inits,
                      0,10000,
                      1/1000,1/1000,
                      R=R,
                      W=W_b,
                      O=O_b,
                      V_mu=1000,
                      tuning=0.02,
                      model="3")

simple_3_15 <- simple(1000,
                      y=y_normalize_2015,
                      inits,
                      0,10000,
                      1/1000,1/1000,
                      R=R,
                      W=W_b,
                      O=O_b,
                      V_mu=1000,
                      tuning=0.02,
                      model="3")

#Model 4: CAR model 3C.
#For this model the starting value for "a" should be equal to 1. 
simple_4_14 <- simple(n_iter=1000,
                      y=y_normalize_2014,
                      inits,
                      0,10000,
                      1/1000,1/1000,
                      W=PSI_2014,
                      V_mu=1000,
                      tuning=0.2,
                      model="4")

simple_4_15 <- simple(n_iter=1000,
                      y=y_normalize_2015,
                      inits,
                      0,10000,
                      1/1000,1/1000,
                      W=PSI_2015,
                      V_mu=1000,
                      tuning=0.2,
                      model="4")
#Model 5: CAR Model 3D
#For this model starting value for "a" must be equal to 0. 
simple_5_14 <- simple(1000,
                      y=y_normalize_2014,
                      inits,
                      0,10000,
                      1/1000,1/1000,
                      R=R,
                      W=PSI_2014_c,
                      w=w,
                      gamma=1,
                      V_mu=1000,
                      tuning=0.02,
                      model="5")

#Now I change the value of the fixed parameter gamma to 40
simple_5_14_b <- simple(1000,
                        y=y_normalize_2014,
                        inits,
                        0,10000,
                        1/1000,1/1000,
                        R=R,
                        W=PSI_2014_c,
                        w=w,
                        gamma=40,#Changed from before
                        V_mu=1000,
                        tuning=0.02,
                        model="5")

simple_5_15 <- simple(1000,
                      y=y_normalize_2015,
                      inits,
                      0,10000,
                      1/1000,1/1000,
                      R=R,
                      W=PSI_2015_c,
                      w=w,
                      gamma=1,
                      V_mu=1000,
                      tuning=0.02,
                      model="5")
simple_5_15_b <- simple(100,
                        y=y_normalize_2015,
                        inits,
                        0,10000,
                        1/1000,1/1000,
                        R=R,
                        W=PSI_2015_c,
                        w=w,
                        gamma=40,#Changed from before
                        V_mu=1000,
                        tuning=0.02,
                        model="5")
