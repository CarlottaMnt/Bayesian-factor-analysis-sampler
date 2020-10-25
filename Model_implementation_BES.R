setwd("C:/Users/carlo/Dropbox/Casa/Otta/ISTA_project/ISTA_project")
setwd("C:/Users/carlo/Documents/ISTAT_project/ISTA_project")
# install.packages("spdep")
# install.packages("knitr")
# install.packages("bestNormalize")
# install.packages("dlm")
# install.packages("xtable")
#install.packages("rstan")

library(ggmap)
library(tmaptools)
library(tmap)
library(tidyverse)
library(tidybayes)
library(spdep)
library(knitr)
library(bestNormalize)
library(dlm)
library(xtable)
library(plyr)
library(gtable)
library(grid)
library(gridExtra)
library(rstan)
#devtools::install_github("karthik/wesanderson")
#library(wesanderson)
rm(list=ls())
source("simple.R")
#Read_in the data
data<- read.csv("BES_data/Imputed_BES_dataset.csv")
data <- data %>% 
  mutate(anno=sapply(strsplit(Anno, split='_', fixed=TRUE), function(x) (x[2]))) %>% 
  select(-Anno)
location <- readxl::read_xlsx("BES_data/Location_BES.xlsx")
data <- cbind(data,location)
data_norm <- read.csv("BES_data/Imputed_BES_normalized.csv")

#factor_eco <- read.csv("factor_econ.csv")
italgeo <-  sf::st_read("BES_data/italgeo.shp")
italgeo$DEN_UTS <- as.character(italgeo$DEN_UTS)
#English names of the elementary indicato
names <- vector("list",length=4)
names(names) <- c("Overall","Social","Economic","Environmental")
names[[1]] <-  c("Employees in cultural business",
                 "Prison density", "Other reported crimes","Youth (<40 years old) political representation in municipalities",
                 "Women's political representation in municipalities","Collection capacity: provincial governments",
                 "Children who benefited of early childhood services","Collection capacity: municipalities",
                 "Landfill of urban waste","Widespread crimes reported","Density of historical green areas",
                 "Dissemination of holidays farm","Availability of urban green areas","Regional health services outflows (hospital admittances)",
                 "Energy from renewable sources","Working days of paid (employee)",
                 "People not in education, employment, or training (Neet)","Average yearly per capita pension income",
                 "Irregular electricity services","People having completed tertiary education (25-34 years)",
                 "Graduates mobility (25-39 years)","Infant mortality rate",
                 "Age-standardised mortality rate for dementia and nervous diseases",
                 "Roads accidents mortality rate (15-34 years old)","Age-standardised cancer mortality rate",
                 "Mortality rate in extra-urban road accidents","Homicide rate","Participation in long-life learning",
                 "Participation to childhood school",
                 "Electoral participation (european elctions)",
                 "Electoral participation (regional elections)",
                 "Per capita asset","Pensioners with low-pension",
                 "People with at least upper secondary education level (25-64 years)",
                 "Public transport network","Separate collection of municipal waste",
                 "Average disposable per capita income",
                 "Average yearly earnings of employee",
                 "Life expectancy at birth","Mortal accidents and inabilities rate",
                 "Rate of bank's non-performing loans to households","Non-participation rate",
                 "Youth non-participation rate (15-29 years)","Employment rate(20-64 years)",
                 "Youth employment rate (15-29 years)"
)
names[[2]] <- c("Infant mortality rate","Life expectancy at birth",
                   "Roads accidents mortality rate (15-34 years old)",
                   "Employees in cultural business",
                   "Other reported crimes",
                   "Age-standardised cancer mortality rate",
                   "Youth (<40 years old) political representation in municipalities",
                   "Age-standardised mortality rate for dementia and nervous diseases",
                   "Participation in long-life learning",
                   "People with at least upper secondary education level (25-64 years)",
                   "People not in education, employment, or training (Neet)",
                   "People having completed tertiary education (25-34 years)",
                   "Participation to childhood school",
                   "Electoral participation (european elctions)",
                   "Electoral participation (regional elections)",
                   "Women's political representation in municipalities",
                   "Prison density",
                   "Graduates mobility (25-39 years)",
                   "Collection capacity: municipalities",
                   "Collection capacity: provincial governments",
                   "Homicide rate",
                   "Widespread crimes reported",
                   "Mortality rate in extra-urban road accidents",
                   "Children who benefited of early childhood services",
                   "Irregular electricity services",
                   "Public transport network",
                   "Regional health services outflows (hospital admittances)",
                   "Dissemination of holidays farm",
                   "Density of historical green areas")
names[[3]] <- c("Average disposable per capita income",
                   "Average yearly earnings of employee",
                   "Average yearly per capita pension income",
                   "Pensioners with low-pension",
                   "Per capita asset",
                   "Rate of bank's non-performing loans to households",
                   "Employment rate(20-64 years)",
                   "Non-participation rate",
                   "Youth non-participation rate (15-29 years)",
                   "Youth employment rate (15-29 years)",
                   "Working days of paid (employee)",
                   "Mortal accidents and inabilities rate"
                   
)
names[[4]] <- c("Landfill of urban waste",
                  "Availability of urban green areas",
                  "Energy from renewable sources",
                  "Separate collection of municipal waste"
)

#Assign the names to the data
colnames(data)[2:46] <- names[[1]]
colnames(data_norm)[1:45] <- names[[1]]
#Remove absent provinces and uniformy the names
data <- data %>% filter(!(comune %in% c("Olbia-Tempio","Medio Campidano","Ogliastra")))
data <- data %>% 
  mutate(comune= case_when(
    grepl("ForlÃ¬-Cesena",comune) ~ "Forli'-Cesena",
    grepl("Carbonia-Iglesias",comune)~ "Sud Sardegna",
    TRUE ~ as.character(comune)))
data_norm <- data_norm %>% 
  mutate(anno = data$anno, comune=data$comune)
#Italy geo locations
italgeo <- italgeo %>% 
  mutate(comune=DEN_UTS)
#Remove unnecessary columns
italgeo <- italgeo %>% 
  dplyr:: select(-c("COD_CM","TIPO_UTS","DEN_PROV","SIGLA","COD_RIP"))
#Replace names
italgeo <- italgeo %>% 
  mutate(comune=case_when(
    grepl("Massa Carrara",comune)~ "Massa-Carrara",
    grepl("Bolzano",comune)~ "Bolzano Bozen",
    TRUE ~ as.character(comune)
  ))
#Spatial correlation matrices---------------------------------------------------
#Model 2
longlat04 <- data %>% 
  filter(anno=="2004") %>% 
  dplyr::select(`data$long`,`data$lat`) %>% as.matrix()
PSI_04 <- as.matrix(dist(longlat04))  #Euclidean Distance
PSI_04.inv <- 1/PSI_04  #Inverse of the distances:weights
diag(PSI_04.inv) <- 0

#Model CAR 2-3

w <- poly2nb(italgeo,row.names=as.character(italgeo$comune),queen=TRUE)
class(w)
summary(w)
R <- nb2mat(w, style='B',zero.policy=TRUE)
N <- dim(R)[1]

#Create the matrix for model 3B
o_b<- rep(NA,N)
Xi <- eigen(R)
xi <- sort(Xi$values)
xi_inv <- 1/(xi)
#if (xi_inv[1] > x  | x > xi_inv[N])
#  x=0.12
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
#Model 3D
PSI_04_d <- nbdists(w, coords= longlat04,longlat = TRUE) #Crate the list of neighbors

#Model implementation-------------------------------------------------------

models <-  c("0","1","2","3","4","5","6")
years <- c("2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017")
domains <- c("Social","Economic","Environmental")
Seed=123
comune <-  data %>%
       gather(Indicatore, Value,-c(anno,comune)) %>% 
       filter(Indicatore %in% names[[1]]) %>% 
       spread(Indicatore,Value) %>% 
       dplyr::filter(anno=="2004") %>% 
       select(comune)
province <- as.vector(comune$comune)

spmatrices <- vector("list",length = 7)
names(spmatrices) <- models
  for (m in models)
  {
    spmatrices[[m]] <- vector("list",length=4)
    names(spmatrices[[m]]) <- c("R","W","O","w")
  }
spmatrices[["1"]][["W"]] <- PSI_04
spmatrices[["2"]][["R"]] <- R
spmatrices[["2"]][["O"]] <- diag(N)
spmatrices[["3"]][["R"]] <- R
spmatrices[["3"]][["W"]] <- W_b
spmatrices[["3"]][["O"]] <- O_b
spmatrices[["4"]][["W"]] <- PSI_04
spmatrices[["5"]][["R"]] <- R
spmatrices[["5"]][["W"]] <- PSI_04_d
spmatrices[["5"]][["w"]] <- w
spmatrices[["6"]][["R"]] <- R
spmatrices[["6"]][["W"]] <- PSI_04_d
spmatrices[["6"]][["w"]] <- w

#OUTCOME VARIABLE
y <- vector("list",length=4) #Basic
y_norm <- vector("list",length=4) #Normalized
y_scale <- vector("list",length=4) #Scaled

names(y) <- domains
names(y_norm) <- domains
names(y_scale) <- domains
		     
for (d in domains)
{
  y[[d]] <- vector("list",length=length(years))
  names(y[[d]]) <- years 
  y_norm[[d]] <- vector("list",length=length(years))
  names(y_norm[[d]]) <- years 
  for (t in years)
  {
     y[[d]][[t]] <- data %>%
       tidyr::gather(Indicatore, Value,-c(anno,comune,`data$long`,`data$lat`)) %>% 
       filter(Indicatore %in% names[[d]]) %>% 
       tidyr::spread(Indicatore,Value) %>% 
       dplyr::filter(anno==t) %>% 
       select(-c(comune,anno,`data$long`,`data$lat`)) %>% 
       as.matrix()
     y_norm[[d]][[t]] <- data_norm %>%
       tidyr::gather(Indicatore, Value,-c(anno,comune)) %>% 
       filter(Indicatore %in% names[[d]]) %>% 
       tidyr::spread(Indicatore,Value) %>% 
       dplyr::filter(anno==t) %>% 
       select(-c(anno,comune)) %>% 
       as.matrix()
     y_scale[[d]][[t]] <- data %>%
       tidyr::gather(Indicatore, Value,-c(anno,comune)) %>% 
       filter(Indicatore %in% names[[d]]) %>% 
       tidyr::spread(Indicatore,Value) %>% 
       dplyr::filter(anno==t) %>% 
       select(-c(anno,comune)) %>% 
       as.matrix() %>% 
       scale(.,center = FALSE,scale=apply(., 2, sd, na.rm = TRUE))
  }
}

#Implementaiton of model "0","1","2","3"
		     
for (d in domains)
{
  for (m in models[models %in% c("0","1","2","3")])
  {
    for (t in years)
    {
      N=dim(y_scale[[d]][[t]])[1]
      J=dim(y_scale[[d]][[t]])[2]
      L=1
      inits <- list(lambda=matrix(1,J,L),
                    mu=matrix(1,J,1),
                    Sigma=diag(10,J),
                    PSI=diag(N),
                    delta=matrix(1,nrow=N,ncol=1),
                    a= ifelse (m %in% c("1","4"),1,0)
      )
      set.seed(Seed + as.numeric(m) + ifelse(m %in% c("3","5"),2,1))
	  xx <- sprintf("Data/output2_%s_%s_%s.RData",d,m,t)
	  readline(paste(d, m, t))
      output <- simple(n_iter = 23500,burn_in = 16000,
                                      y = y_scale[[d]][[t]],
                                      inits,
                                      0,10000,
                                      1/1000,1/1000,
                                      R=spmatrices[[m]][["R"]],
                                      J=J,
                                      W=spmatrices[[m]][["W"]],
                                      O=spmatrices[[m]][["O"]],
                                      w=spmatrices[[m]][["w"]],
                                      V_mu=1000,
                                      gamma = ifelse(m == "5",1,40), 
                                      mu_a = ifelse (m %in% c("1","4"),1,0),
                                      tuning=ifelse(m == "1",sqrt(1000), sqrt(0.4)),
                                      model = m)
	save(output,file=xx);
   
    }
  }
}


dev.off()



##Model selection----------------------------------------------------------------
y_rep <- vector("list",length=length(domains))
sd <-  vector("list",length=length(domains))
G <- vector("list",length=length(domains))
p <- vector("list",length=length(domains))
C <- vector("list",length=length(domains))
G_model <- vector("list",length =length(domains))
C_model <- vector("list",length =length(domains))
P_model <- vector("list",length =length(domains))
names(sd) <- domains
names(y_rep) <- domains
names(G)<- domains
names(p) <- domains
names(C) <- domains
names(G_model) <- domains
names(C_model) <- domains
names(P_model) <- domains

for (d in domains)
{
  y_rep[[d]] <- vector("list",length=length(models))
  names(y_rep[[d]]) <- models
  sd[[d]] <- vector("list",length=length(models))
  names(sd[[d]]) <- models
  G[[d]] <- vector("list",length=length(models))
  names(G[[d]]) <- models
  p[[d]] <- vector("list",length=length(models))
  names(p[[d]]) <- models
  C[[d]] <- vector("list",length=length(models))
  names(C[[d]]) <- models
  for (m in models)
  {
    y_rep[[d]][[m]] <- vector("list",length=length(years))
    names(y_rep[[d]][[m]]) <- years
    sd[[d]][[m]] <- vector("list",length=length(years))
    names(sd[[d]][[m]]) <- years
    G[[d]][[m]] <- vector("list",length=length(years))
    names(G[[d]][[m]]) <- years
    p[[d]][[m]] <- vector("list",length=length(years))
    names(p[[d]][[m]]) <- years
    C[[d]][[m]] <- vector("list",length=length(years))
    names(C[[d]][[m]]) <- years
  }
}
for (d in domains)
  {
  for (m in models[models %in% c("0","1","2","3")])
    {
    for (t in years)
      {
          xx <- sprintf("Data/output2_%s_%s_%s.RData",d,m,t)
          load(file=xx);
          
          N=dim(y[[d]][[t]])[1]
          J=ifelse((d == "Social" & t == "2017"),
          dim(y_center[[d]][[t]])[2]-2,
          dim(y_center[[d]][[t]])[2])
          y_rep[[d]][[m]][[t]] <- apply(output$y_rep,2,mean) %>% matrix(N,J,byrow=TRUE)
          sd[[d]][[m]][[t]] <- apply(output$y_rep,2,var) %>% matrix(N,J,byrow=TRUE)
          G[[d]][[m]][[t]] <- sum(rowSums((y_rep[[d]][[m]][[t]]- ifelse((d == "Social" & t == "2017"),y_center[[d]][[t]][,-c(1,16)],y_center[[d]][[t]]))^2))
          p[[d]][[m]][[t]] <- sum(rowSums(sd[[d]][[m]][[t]]))
          C[[d]][[m]][[t]] <- G[[d]][[m]][[t]] + p[[d]][[m]][[t]]
      }
    }
}

G_model[[d]][[t]] <-rep(NA,length(models))
P_model[[d]][[t]] <- rep(NA,length(models))
C_model[[d]][[t]] <- rep(NA,length(models))
for (d in domains[1])
{
  for (t in years)
  {
    for (m in models)
    {
      G_model[[d]][[t]]= c(G[[d]][[m]][[t]])
      P_model[[d]][[t]]=c(p[[d]][[m]][[t]])
      C_model[[d]][[t]]=c(C[[d]][[m]][[t]])
    }
  }
}
  

#Factor loadings----------------------------------------------------------------
lambda <- vector("list",length=length(domains))
mu <- vector("list",length=length(domains))
summary <- vector("list",length=length(domains))
mu_summary <- vector("list",length=length(domains))
names(lambda) <- domains
names(mu) <- domains
names(summary) <- domains
names(mu_summary) <- domains
for (d in domains)
{
  lambda[[d]] <- vector("list",length=length(models))
  names(lambda[[d]]) <- models
  mu[[d]] <- vector("list",length=length(models))
  names(mu[[d]]) <- models
  summary[[d]] <- vector("list",length=length(models))
  names(summary[[d]]) <- models
  mu_summary[[d]] <- vector("list",length=length(models))
  names(mu_summary[[d]]) <- models
  for (m in models)
  {
    lambda[[d]][[m]] <- vector("list",length=length(years))
    names(lambda[[d]][[m]]) <- years
    mu[[d]][[m]] <- vector("list",length=length(years))
    names(mu[[d]][[m]]) <- years
    summary[[d]][[m]] <- vector("list",length=length(years))
    names(summary[[d]][[m]]) <- years
    mu_summary[[d]][[m]] <- vector("list",length=length(years))
    names(mu_summary[[d]][[m]]) <- years
  }
}
for (d in domains)
{
  for (m in models[models %in% c("0","1","2","3")])
  {
    for (t in years[years %in% c("2017")])
      {
          xx <- sprintf("Data/output2_%s_%s_%s.RData",d,m,t)
          load(file=xx);
          #cat(xx)
      lambda <- output$lambda %>% as.data.frame()
      mu <- output$mu %>% as.data.frame()
      summary[[d]][[m]][[t]] <- rbind(apply(lambda, 2,mean),apply(lambda,2, quantile,c(.025, .50, .975)),apply(lambda,2, IQR))
      mu_summary[[d]][[m]][[t]] <- rbind(apply(mu, 2,mean),apply(mu,2, quantile,c(.025, .50, .975)),apply(mu,2, IQR))
      if (d == "Social" & t =="2017")
      {
      colnames(summary[[d]][[m]][[t]]) <- names[[d]][-c(1,16)]
      colnames(mu_summary[[d]][[m]][[t]]) <- names[[d]][-c(1,16)]
      } else 
      {
      colnames(summary[[d]][[m]][[t]]) <- names[[d]]
      colnames(mu_summary[[d]][[m]][[t]]) <- names[[d]]    
      }      
    }
  }
}
#variances------------------------------------------------------------
sigma <- vector("list",length=length(domains))
f <- vector("list",length=length(domains))
names(sigma) <- domains
names(f) <- domains
for (d in domains)
{
  sigma[[d]] <- vector("list",length=length(models))
  names(sigma[[d]]) <- models
  f[[d]] <- vector("list",length=length(models))
  names(f[[d]]) <- models
  for (m in models)
  {
    sigma[[d]][[m]] <- vector("list",length=length(years))
    names(sigma[[d]][[m]]) <- years
    f[[d]][[m]] <- vector("list",length=length(years))
    names(f[[d]][[m]]) <- years
    for (t in years)
    { 
      J= ifelse((d == "Social" & t =="2017"),dim(y_center[[d]][[t]])[2]-2,dim(y_center[[d]][[t]])[2])
      f[[d]][[m]][[t]] <- matrix(NA,5,J)
      if (d == "Social" & t =="2017")
      {
      colnames(f[[d]][[m]][[t]]) <- names[[d]][-c(1,16)]
      } else colnames(f[[d]][[m]][[t]]) <- names[[d]]
    }
  }
}

for (d in domains)
{
  for (m in models[models %in% c("0","1","2","3")])
  {
    for (t in years[years %in% c("2017")])
    { 
      xx <- sprintf("Data/output2_%s_%s_%s.RData",d,m,t)
      load(file=xx)
      cat(xx)
      J=ifelse((d == "Social" & t == "2017"),dim(output$lambda)[2]-2,dim(output$lambda)[2])
      sigma[[d]][[m]][[t]] <- rbind(apply(output$Sigma,2,mean),apply(output$Sigma,2, quantile,c(.05, .50, .95)),apply(output$Sigma,2,IQR))
      if (d == "Social" & t =="2017")
      {
      colnames(sigma[[d]][[m]][[t]]) <- names[[d]][-c(1,16)]
      }else colnames(sigma[[d]][[m]][[t]]) <- names[[d]]
    }
  }
}

for (d in domains)
{
  for (m in models[models %in% c("0","1","2","3")])
  {
    for (t in years)
    { 
      J= dim(output$lambda)[2]
      for (i in if(d == "Social" & t == "2017")
      {names[[d]][-c(1,16)]}else names[[d]])
      {
        f[[d]][[m]][[t]][,i] <- (summary[[d]][[m]][[t]][,i]^2)/(summary[[d]][[m]][[t]][,i]^2+ sigma[[d]][[m]][[t]][,i])
        row.names(f[[d]][[m]][[t]]) <- c("mean", "2.5%","50%","97.5%","IQR")
      }
    }
  }
}

#Spatial parameter--------------------------------------------------------------
a <- vector("list",length=length(domains))
names(a) <- domains
for (d in domains)
{
  a[[d]] <- vector("list",length=length(models[-1]))
  names(a[[d]]) <- models[-c(1,5,6,7)]
  for (m in models[-1])
  {
    a[[d]][[m]] <- vector("list",length=length(years))
    names(a[[d]][[m]]) <- years
  }
}
for (d in domains)
{
  for (m in models[models %in% c("1","2","3")])
  {
    for (t in years)
    {
      xx <- sprintf("Data/output2_%s_%s_%s.RData",d,m,t)
      load(file=xx)
      a[[d]][[m]][[t]] <- output$a %>% 
        as.data.frame() %>% 
        dplyr::summarise(min=min(.),first_qu = quantile(.,0.25),
                         mean=mean(.),median=median(.),third_qu=quantile(.,0.75),
                         max=max(.),sd=sd(.),IQR=IQR(.),year=t,model= m,domain = d) %>% 
        as.list()
    }
  }
}
a_df <- vector("list",length=length(domains))
names(a_df) <- domains
for (d in domains)
{
  a_df[[d]] <- vector("list",length=length(models[c(2,3,4)]))
  names(a_df[[d]]) <-  models[c(2,3,4)]
}
for (d in domains)
{
  for (m in models[c(2,3,4)])
  {
    a_df[[d]][[m]] <- do.call(rbind.data.frame, a[[d]][[m]])
  }
}
a_df_d <- vector("list",length=length(domains))
names(a_df_d) <- domains
for (d in domains)
{
  a_df_d[[d]] <- do.call(rbind.data.frame, a_df[[d]])
}
a_df_d_tot <- do.call(rbind.data.frame, a_df_d)
a_df_d_tot$model <- plyr::revalue(x = a_df_d_tot$model, 
c("1" = "Marginal correlation", "2" = "CAR 3A", "3" = "CAR 3B"))
a1 <- a_df_d_tot %>% filter(domain == "Social", model == "Marginal correlation") %>% select(9,1:8)
a2 <- a_df_d_tot %>% filter(domain == "Social", model == "CAR 3A") %>% select(1:8)
a3 <- a_df_d_tot %>% filter(domain == "Social", model == "CAR 3B") %>% select(1:8)
row.names(a1) <- NULL
row.names(a2) <- NULL
row.names(a3) <- NULL
print(xtable(cbind(a1,a2,a3), digits = 2),format.arg = list(big.mark = " ", decimal.mark = "."))

#Comparison tra domini preso lo stesso modello
for (m in models)
{
  a_df_m[[m]] <- do.call(rbind.data.frame, a_df[[m]])
}

pd <- position_dodge(0.1) 

p <- ggplot(a_df_d_tot %>% filter(model %in% c("Marginal correlation","CAR 3A","CAR 3B")) ,aes(y=mean,x=year,col = domain))+
    geom_errorbar(aes(ymin=first_qu, ymax=third_qu), width=.4,lty=1, position=pd) +
    geom_point(size=1)+
    facet_wrap(. ~ model, scales="free_y",ncol= 1)+
    geom_line(aes(colour = domain, linetype = domain,group = domain)) +
    geom_hline(yintercept=0, linetype="dashed")+
    #scale_y_discrete(expand = expansion(mult = c(0, .02)))+
    #guides(fill=FALSE, color=FALSE)+
    ylab(expression(omega))+
    theme_bw() +
    theme(strip.background = element_rect(fill = NA, color = "black"),
          axis.text.x = element_text(size = 8, angle = 45,vjust=0.5),
          plot.caption = element_text(hjust = 1))
png("plotspatial.png")
print(p)     
dev.off() 
egg::ggarrange(p1, p2,p3,p4,p5 ncol = 1,heights = c(1,10))
p[["Overall"]]
p[["Social"]]
p[["Economic"]]
p[["Environmental"]]

#Delta (composite indicator)----------------------------------------------------
multi.sapply <- function(...) {
  arglist <- match.call(expand.dots = FALSE)$...
  var.names <- sapply(arglist, deparse)
  has.name <- (names(arglist) != "")
  var.names[has.name] <- names(arglist)[has.name]
  arglist <- lapply(arglist, eval.parent, n = 2)
  x <- arglist[[1]]
  arglist[[1]] <- NULL
  result <- sapply(arglist, function (FUN, x) sapply(x, FUN), x)
  colnames(result) <- var.names[-1]
  return(result)
}
delta <- vector("list",length=length(domains))
d_summary <- vector("list",length=length(domains))
names(delta) <- domains
names(d_summary) <- domains
for (d in domains)
{
  delta[[d]] <- vector("list",length=length(models[-c(5,6,7)]))
  names(delta[[d]]) <- models[-c(5,6,7)]
  d_summary[[d]] <- vector("list",length=length(models[-c(5,6,7)]))
  names(d_summary[[d]]) <- models[-c(5,6,7)]
  for (m in models[-c(5,6,7)])
  {
    delta[[d]][[m]] <- vector("list",length=length(years))
    names(delta[[d]][[m]]) <- years
    d_summary[[d]][[m]] <- vector("list",length=length(years))
    names(d_summary[[d]][[m]]) <- years
  }
}
for (d in domains)
{
  for (m in models[-c(5,6,7)])
  {
    for (t in years[years %in% c("2017")])
    {
      xx <- sprintf("Data/output2_%s_%s_%s.RData",d,m,t)
      load(file=xx)
      delta[[d]][[m]][[t]] <- output$d %>% as.data.frame()
      d_summary[[d]][[m]][[t]] <- multi.sapply(delta[[d]][[m]][[t]], mean, median , third_qu=function(x) quantile(x, 0.75),firt_qu= function(x) quantile(x,0.25), IQR = IQR) %>%
       as.data.frame()%>%
       mutate(comune = province 
       , anno=t , model = m) %>% 
        as.list()
    }
  }
}

#Plot
plot_delta <- vector("list",length=length(domains))
names(plot_delta) <- domains
delta_df <- vector("list",length=length(domains))
names(delta_df) <- domains
for (d in domains)
{
  plot_delta[[d]] <- vector("list",length=length(models[-c(5,6,7)]))
  names(plot_delta[[d]]) <- models[-c(5,6,7)]
  delta_df[[d]] <- vector("list",length=length(models[-c(5,6,7)]))
  names(delta_df[[d]]) <- models[-c(5,6,7)]
  for (m in models[-c(5,6,7)])
  {
    plot_delta[[d]][[m]] <- vector("list",length=length(years))
    names(plot_delta[[d]][[m]]) <- years
  }
}

for (d in domains)
{
  for (m in models[-c(5,6,7)])
  {
   delta_df[[d]][[m]] <- do.call(rbind.data.frame, d_summary[[d]][[m]])
  }
}
delta_df_d <- vector("list",length=length(domains))
names(delta_df_d) <- domains
for (d in domains)
{
  delta_df_d[[d]] <- do.call(rbind.data.frame, delta_df[[d]])
}



social0 <- delta_df_d[["Social"]] %>% filter (anno == "2017", model %in% c(0), comune %in% c(province[1:34])) %>%  
        dplyr::select(c(6,1:5))
social1 <- delta_df_d[["Social"]] %>% filter (anno == "2017", model %in% c(1), comune %in% province[1:34]) %>% dplyr::select(c(1:5))

social01 <- delta_df_d[["Social"]] %>% filter (anno == "2017", model %in% c(0), comune %in% province[35:77]) %>%  
        dplyr::select(c(6,1:5))
social11 <- delta_df_d[["Social"]] %>% filter (anno == "2017", model %in% c(1), comune %in% province[35:77]) %>% dplyr::select(c(1:5))
row.names(social0) <- NULL
row.names(social1) <- NULL
row.names(social01) <- NULL
row.names(social11) <- NULL
kable(cbind(cbind(social0,social1),cbind(social01,social11)),format="latex", digit=2)


economic0 <- delta_df_d[["Economic"]] %>% filter (anno == "2017", model %in% c(0), comune %in% c(province[104:107])) %>% dplyr::select(c(6,1:5))
economic1 <- delta_df_d[["Economic"]] %>% filter (anno == "2017", model %in% c(3),comune %in% c(province[73:88]))%>% dplyr::select(c(1:5))
economic01 <- delta_df_d[["Economic"]] %>% filter (anno == "2017", model %in% c(0), comune %in% c(province[89:104])) %>% dplyr::select(c(6,1:5))
economic11 <- delta_df_d[["Economic"]] %>% filter (anno == "2017", model %in% c(3),comune %in% c(province[89:104])) %>% dplyr::select(c(1:5))
row.names(economic0) <- NULL
row.names(economic1 ) <- NULL
row.names(economic01) <- NULL
row.names(economic11) <- NULL
kable(cbind(cbind(economic0,economic1),cbind(economic01,economic11)),format="latex", digit=2)

environ0 <- delta_df_d[["Environmental"]] %>% filter (anno == "2017", model %in% c(0), comune %in% c(province[104:107])) %>% dplyr::select(c(6,1:5))
environ1 <- delta_df_d[["Environmental"]] %>% filter (anno == "2010", model %in% c(2), comune %in% c(province[73:88])) %>% dplyr::select(c(1:5))
environ01 <- delta_df_d[["Environmental"]] %>% filter (anno == "2017", model %in% c(0), comune %in% c(province[89:104])) %>% dplyr::select(c(6,1:5))
environ11 <- delta_df_d[["Environmental"]] %>% filter (anno == "2010", model %in% c(2), comune %in% c(province[89:104])) %>% dplyr::select(c(1:5))
row.names(environ0) <- NULL
row.names(environ1) <- NULL
row.names(environ01) <- NULL
row.names(environ11) <- NULL
kable(cbind(cbind(environ0,environ1),cbind(environ01,environ11)) ,format="latex", digit=2)


 #2006
 delta_df_d[[3]]$mean <- ifelse(delta_df_d[[3]]$anno %in% c("2006"),-delta_df_d[[3]]$mean, delta_df_d[[3]]$mean)
 delta_df_d[[3]]$firt_qu <- ifelse(delta_df_d[[3]]$anno %in% c("2006"),-delta_df_d[[3]]$firt_qu, delta_df_d[[3]]$firt_qu)
delta_df_d[[3]]$third_qu <- ifelse(delta_df_d[[3]]$anno %in% c("2006"),-delta_df_d[[3]]$third_qu, delta_df_d[[3]]$third_qu)
 #2005
 delta_df_d[[3]]$mean <- ifelse(delta_df_d[[3]]$anno %in% c("2005"),-delta_df_d[[3]]$mean, delta_df_d[[3]]$mean)
 delta_df_d[[3]]$firt_qu <- ifelse(delta_df_d[[3]]$anno %in% c("2005"),-delta_df_d[[3]]$firt_qu, delta_df_d[[3]]$firt_qu)
 delta_df_d[[3]]$third_qu <- ifelse(delta_df_d[[3]]$anno %in% c("2005"),-delta_df_d[[3]]$third_qu, delta_df_d[[3]]$third_qu)
#2007
delta_df_d[[3]]$mean <- ifelse((delta_df_d[[3]]$anno %in% c("2007") & delta_df_d[[3]]$model %in% c("0")),-delta_df_d[[3]]$mean, delta_df_d[[3]]$mean)
delta_df_d[[3]]$firt_qu <- ifelse((delta_df_d[[3]]$anno %in% c("2007") & delta_df_d[[3]]$model %in% c("0")),-delta_df_d[[3]]$firt_qu, delta_df_d[[3]]$firt_qu)
delta_df_d[[3]]$third_qu <- ifelse((delta_df_d[[3]]$anno %in% c("2007") & delta_df_d[[3]]$model %in% c("0")),-delta_df_d[[3]]$third_qu, delta_df_d[[3]]$third_qu)
# #2008
 delta_df_d[[3]]$mean <- ifelse((delta_df_d[[3]]$anno %in% c("2008") & delta_df_d[[3]]$model %in% c("1","3")),-delta_df_d[[3]]$mean, delta_df_d[[3]]$mean)
 delta_df_d[[3]]$firt_qu <- ifelse((delta_df_d[[3]]$anno %in% c("2008") & delta_df_d[[3]]$model %in% c("1","3")),-delta_df_d[[3]]$firt_qu, delta_df_d[[3]]$firt_qu)
 delta_df_d[[3]]$third_qu <- ifelse((delta_df_d[[3]]$anno %in% c("2008") & delta_df_d[[3]]$model %in% c("1","3")),-delta_df_d[[3]]$third_qu, delta_df_d[[3]]$third_qu)
# #2009
# delta_df_d[[3]]$mean <- ifelse(delta_df_d[[3]]$anno %in% c("2009"),-delta_df_d[[3]]$mean, delta_df_d[[3]]$mean)
# delta_df_d[[3]]$firt_qu <- ifelse(delta_df_d[[3]]$anno %in% c("2009"),-delta_df_d[[3]]$firt_qu, delta_df_d[[3]]$firt_qu)
# delta_df_d[[3]]$third_qu <- ifelse(delta_df_d[[3]]$anno %in% c("2009"),-delta_df_d[[3]]$third_qu, delta_df_d[[3]]$third_qu)
# #2010
# delta_df_d[[3]]$mean <- ifelse((delta_df_d[[3]]$anno %in% c("2010")& delta_df_d[[3]]$model %in% c("1")),-delta_df_d[[3]]$mean, delta_df_d[[3]]$mean)
# delta_df_d[[3]]$firt_qu <- ifelse((delta_df_d[[3]]$anno %in% c("2010")& delta_df_d[[3]]$model %in% c("1")),-delta_df_d[[3]]$firt_qu, delta_df_d[[3]]$firt_qu)
# delta_df_d[[3]]$third_qu <- ifelse((delta_df_d[[3]]$anno %in% c("2010")& delta_df_d[[3]]$model %in% c("1")),-delta_df_d[[3]]$third_qu, delta_df_d[[3]]$third_qu)
#2011
delta_df_d[[3]]$mean <- ifelse((delta_df_d[[3]]$anno %in% c("2011")& delta_df_d[[3]]$model %in% c("0","1")),-delta_df_d[[3]]$mean, delta_df_d[[3]]$mean)
delta_df_d[[3]]$firt_qu <- ifelse((delta_df_d[[3]]$anno %in% c("2011")& delta_df_d[[3]]$model %in% c("0","1")),-delta_df_d[[3]]$firt_qu, delta_df_d[[3]]$firt_qu)
delta_df_d[[3]]$third_qu <- ifelse((delta_df_d[[3]]$anno %in% c("2011")& delta_df_d[[3]]$model %in% c("0","1")),-delta_df_d[[3]]$third_qu, delta_df_d[[3]]$third_qu)
#2012
delta_df_d[[3]]$mean <- ifelse((delta_df_d[[3]]$anno %in% c("2012")& delta_df_d[[3]]$model %in% c("0")),-delta_df_d[[3]]$mean, delta_df_d[[3]]$mean)
delta_df_d[[3]]$firt_qu <- ifelse((delta_df_d[[3]]$anno %in% c("2012")& delta_df_d[[3]]$model %in% c("0")),-delta_df_d[[3]]$firt_qu, delta_df_d[[3]]$firt_qu)
delta_df_d[[3]]$third_qu <- ifelse((delta_df_d[[3]]$anno %in% c("2012")& delta_df_d[[3]]$model %in% c("0")),-delta_df_d[[3]]$third_qu, delta_df_d[[3]]$third_qu)
#2015
delta_df_d[[3]]$mean <- ifelse((delta_df_d[[3]]$anno %in% c("2015")& delta_df_d[[3]]$model %in% c("1")),-delta_df_d[[3]]$mean, delta_df_d[[3]]$mean)
delta_df_d[[3]]$firt_qu <- ifelse((delta_df_d[[3]]$anno %in% c("2015")& delta_df_d[[3]]$model %in% c("1")),-delta_df_d[[3]]$firt_qu, delta_df_d[[3]]$firt_qu)
delta_df_d[[3]]$third_qu <- ifelse((delta_df_d[[3]]$anno %in% c("2015")& delta_df_d[[3]]$model %in% c("1")),-delta_df_d[[3]]$third_qu, delta_df_d[[3]]$third_qu)
#2016
delta_df_d[[3]]$mean <- ifelse((delta_df_d[[3]]$anno %in% c("2016")& delta_df_d[[3]]$model %in% c("1")),-delta_df_d[[3]]$mean, delta_df_d[[3]]$mean)
delta_df_d[[3]]$firt_qu <- ifelse((delta_df_d[[3]]$anno %in% c("2016")& delta_df_d[[3]]$model %in% c("1")),-delta_df_d[[3]]$firt_qu, delta_df_d[[3]]$firt_qu)
delta_df_d[[3]]$third_qu <- ifelse((delta_df_d[[3]]$anno %in% c("2016")& delta_df_d[[3]]$model %in% c("1")),-delta_df_d[[3]]$third_qu, delta_df_d[[3]]$third_qu)


for (d in domains[3])
{
  {
   for (t in years)
   {
    plot_delta[[d]][[t]] <- delta_df_d[[d]]%>%     
       filter(anno == t, comune %in% c("Milano","Roma","Torino","Genova","Palermo","Cagliari",
                                               "Trento","Venezia","Perugia","Potenza","Bari",
                                              "Napoli","Bologna","Firenze","Aosta","Reggio di Calabria",
                                             "Ancona","Trieste","Campobasso","Catanzaro"))%>%
      #tidyr::gather(comune, Wellbeing) %>% 
      ggplot(aes(y= reorder(comune, mean),x = mean, fill = model))+
      geom_bar(stat = "identity",alpha = 0.5, position = position_dodge())+
      geom_point(position = position_dodge(0.9), alpha = 0.8, size = 2 )+
      #facet_wrap(vars(model),scales='free')+
      geom_errorbar(aes(xmin=firt_qu, xmax=third_qu), width=.2,lty=1,alpha = 0.5, position=position_dodge(1))+
      # geom_text(aes(label= ifelse(mean %in% min(mean)
                                 # ,as.character(comune),'')),col="black",size = 5,hjust = -.05,vjust=0)+
      # geom_text(aes(label= ifelse(mean %in% c(max(mean), median(mean))
                                  
                                  # ,as.character(comune),'')),col="black",size =5, hjust = 1.08,vjust=0)+
      geom_vline(xintercept=0, linetype="dashed")+
      #scale_y_discrete(expand = expansion(mult = c(0, .02)))+
      #guides(fill=FALSE, color=FALSE)+
      xlab(paste0(t,"\nComposite indicator"))+
      ylab("")+
      theme_bw() +
      scale_fill_discrete(name  = "Models",
                           breaks = c("0", "1","2","3"),
                           labels = c("Spatial Independence","Marginal correlation", "CAR 3A", "CAR 3B"),
	            guide = guide_legend()) +
      scale_shape_discrete(name  = "Models",
                           breaks = c("0", "1","2","3"),
                           labels = c("Spatial Independence","Marginal correlation", "CAR 3A", "CAR 3B"))+
      theme(strip.background = element_rect(fill = NA, color = "black"),
            plot.caption = element_text(hjust = 1),
            text = element_text(size=rel(5.5)),
            legend.position="bottom",
            axis.text.y=element_text(size = 25),
            axis.text.x=element_text(size = 25),
            axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
            legend.title=element_text(size=rel(10.5)), 
            legend.text= element_text(size=rel(10.5)),
            axis.ticks.y=element_blank())
  }
 }
}


p1 <- plot_delta[["Economic"]][["2005"]]
p2 <- plot_delta[["Environmental"]][["2006"]]
p3 <- plot_delta[["Environmental"]][["2014"]]
p4 <- plot_delta[["Environmental"]][["2017"]]
legend = gtable_filter(ggplotGrob(plot_delta[["Economic"]][["2017"]]),"guide-box")     # Make sure the legend has been extracted
jpeg("Delta_Environ.jpeg", width=2500, height=(2500))
grid.arrange(arrangeGrob(
 p1 + theme(legend.position="none"),
 p2 + theme(legend.position="none"), 
 p3 + theme(legend.position="none"),
 p4 + theme(legend.position="none"), nrow =1, ncol = 4), legend,  nrow=2, heights=c(20, 5)) 
dev.off()
##MAPS
#VARmap <- merge(italgeo,delta_df_d[[2]] %>% filter(model == "0",anno =="2004"),by = "comune")
VARmap <- vector("list",length=length(domains))
names(VARmap) <- domains
map <- vector("list",length=length(domains))
names(map) <- domains
for (d in domains)
{
  VARmap[[d]] <- vector("list",length=length(models[-c(5,6,7)]))
  names(VARmap[[d]]) <- models[-c(5,6,7)]
  map[[d]] <- vector("list",length=length(models[-c(5,6,7)]))
  names(map[[d]]) <- models[-c(5,6,7)]
  for (m in models[-c(5,6,7)])
  {
    VARmap[[d]][[m]] <- vector("list",length=length(years))
    names(VARmap[[d]][[m]]) <- years
    map[[d]][[m]] <- vector("list",length=length(years))
    names(map[[d]][[m]]) <- years
  }
}

for (d in domains[2])
{ 
  for (m in models[models %in% c("0","3")])
  {
  for (t in years[years %in% c("2004","2006","2014","2017")])
  {
  VARmap[[d]][[m]][[t]] <- merge(italgeo,delta_df_d[[d]] %>% filter(model == m ,anno == t, comune != "Crotone"),by = "comune")
  #xm <- sprintf("Figures/map2_%s_%s_%s.jpeg",d,m,t)
  map[[d]][[t]][[m]] <- tm_shape(VARmap[[d]][[m]][[t]])+
  tm_fill("mean",title= t, style = "quantile", midpoint = NA)+
  tm_borders(alpha=.5)+
  tm_shape(VARmap[[d]][[m]][[t]] %>%  filter(comune %in% c("Milano","Roma","Torino","Genova","Palermo","Cagliari",
                                                 "Trento","Venezia","Perugia","Potenza","Bari",
                                                 "Napoli","Bologna","Firenze","Aosta","Reggio di Calabria",
                                                 "Ancona","Trieste","Campobasso")))+
  tm_dots("comune",size=0.4, col="black", shape= 21, legend.show = FALSE)+
  tm_layout(legend.title.size = 2,
          legend.text.size = 1.6,
          legend.position = c("right","top"),
          legend.bg.color = "white",
          legend.bg.alpha = 1)
  #readline(paste(d,m,t))
  #tmap_save(map, xm, height=7)
  }
 }
}

tmap_save(tmap_arrange(map[[d]][["2004"]][["3"]], map[[d]][["2006"]][["3"]], map[[d]][["2014"]][["3"]],map[[d]][["2017"]][["3"]], nrow = 2, ncol = 2),"ECONOMIC3.jpeg")
save(delta_df_d,file="dataoverall.RData")

########Hierarchical model#########################################################################################
N_O <- c("Savona","Imperia","La Spezia","Genova",
         "Cuneo","Biella","Vercelli","Alessandria",
         "Asti","Verbano-Cusio-Ossola","Novara","Torino",
         "Bergamo","Brescia","Como","Lecco","Cremona","Lodi",
         "Mantova","Milano","Monza e della Brianza","Pavia","Sondrio","Varese",
         "Aosta")
N_E <- c("Bologna","Forli'-Cesena","Ferrara","Modena",
         "Parma","Piacenza","Ravenna","Reggio nell'Emilia",
         "Rimini","Gorizia","Pordenone","Trieste",
         "Udine","Bolzano Bozen","Trento","Belluno","Padova","Rovigo",
         "Treviso","Venezia","Verona","Vicenza"
)
CE <- c("Frosinone","Latina","Rieti","Roma",
        "Viterbo","Arezzo","Firenze","Grosseto",
        "Livorno","Lucca","Massa-Carrara","Pisa",
        "Pistoia","Prato","Siena","Perugia","Terni","Ancona",
        "Ascoli Piceno","Fermo","Macerata","Pesaro e Urbino"
)
SUD <- c("Cosenza","Crotone","Reggio di Calabria",
  "Vibo Valentia","Avellino","Benevento","Caserta","Napoli",
  "Salerno","Campobasso","Isernia","L'Aquila","Pescara","Teramo","Bari",
  "Barletta-Andria-Trani","Chieti", "Foggia","Lecce","Matera","Potenza" ,
  "Taranto","Brindisi")
ISOLE <- c("Cagliari","Nuoro","Oristano","Sassari",
           "Agrigento","Sud Sardegna","Caltanissetta","Catania","Enna",
           "Messina","Palermo","Ragusa","Siracusa",
           "Trapani","Catanzaro")

mean_common_summary <- list()        
for(d in domains[3])
{
 for (t in years)
 {
 for (m in models[models %in% c("0","1","2","3")])
 {
  mean_common_summary[[d]][[t]][[m]] <-  delta_df_d[[d]]%>% as.data.frame() %>% filter(model == m, anno == t) %>% 
    mutate(area = ifelse(comune %in% N_O, "1",
                                                     ifelse(comune %in% N_E, "2",
                                                            ifelse(comune %in% CE, "3",
                                                                   ifelse(comune %in% SUD, "4",
                                                                          ifelse(comune %in% ISOLE, "5",
                                                                                 NA))))))
}
}
}
model <- list()
for(d in domains[3])
{
for (t in years[years %in% c("2008")])
{
for (m in models[models %in% c("1","2","3")])
{
xh <- sprintf("Data/hierarchical_%s_%s_%s.RData",d,m,t)
readline(paste(d, m, t))
y <- mean_common_summary[[d]][[t]][[m]]
   stan_data <- list(
    N=dim(y)[1],
    y= y$mean,
    jj=as.numeric(y$area),
    J=5
  )
model_code <- "data{
    int<lower=1> J;//number of groups
    int<lower=1> N;
    real y[N];
    int<lower=1,upper=J> jj[N];
  }
  parameters {
  real mu[J];
  real<lower=0> Sigma[J];
  real eta;
  real<lower=0>nu;
  }
  transformed parameters{
  real lognu;
  lognu =log(nu);
  }
  model {
  for (i in 1:J) 
  {
  mu[i]~ normal(0, 1);
  Sigma[i] ~ cauchy(eta, nu);
  }
  for (n in 1:N)
  y[n] ~ normal(mu[jj[n]],Sigma[jj[n]]);
  }
  generated quantities{
  real y_pred[N];
  real log_lik[N];
  for (n in 1:N)
  y_pred[n]=normal_rng(mu[jj[n]],1);//fitted values
  for (n in 1:N)
  log_lik[n]=normal_lpdf(y[n]|mu[jj[n]],1);
  }
"
#stan_model <- stan_model(model_code=model_code)

model <- sampling(stan_model,
                      data= stan_data,
                      chains=4,
                      iter=2000,
                      control=list(adapt_delta=0.95,
                                   max_treedepth = 15
                      ))

samples <- extract(model, c("mu","Sigma"))	
save(samples,file=xh);   
}    
}
}

mean_area_overall <- vector("list",length=length(domains))
names(mean_area_overall) <- domains
for (d in domains)
{
  mean_area_overall[[d]] <- vector("list",length=length(models[-c(5,6,7)]))
  names(mean_area_overall[[d]]) <- models[-c(5,6,7)]
  for (m in models[-c(5,6,7)])
  {
    mean_area_overall[[d]][[m]] <- vector("list",length=length(years))
    names(mean_area_overall[[d]][[m]]) <- years
  }
}


for (d in domains)
{
  for (m in models[-c(1,5,6,7)])
  {
    for (t in years)
    {
      xh <- sprintf("Data/hierarchical_%s_%s_%s.RData",d,m,t)
      load(file=xh)
      mean_area <- samples[["mu"]] %>% 
      as.data.frame() 
      mean_area_overall[[d]][[m]][[t]] <- multi.sapply(mean_area, mean, median , third_qu=function(x) quantile(x, 0.75),firt_qu= function(x) quantile(x,0.25), IQR = IQR) %>%
      as.data.frame() %>%
      mutate(area = c("North-West","North-East","Center","South","Islands"),
      anno=t, model = m) %>% 
      as.list()   
    }
  }
}

deltarea_df <- vector("list",length=length(domains))
names(deltarea_df) <- domains
graphdf <- vector("list",length=length(domains))
names(graphdf) <- domains
for (d in domains)
{
  graphdf[[d]] <- vector("list",length=length(years))
  names(graphdf[[d]]) <- years
}
  deltarea_df[[d]] <- vector("list",length=length(models[-c(5,6,7)]))
  names(deltarea_df[[d]]) <- models[-c(5,6,7)]
  for (m in models[-c(5,6,7)])
  {
    graphdf[[d]][[m]] <- vector("list",length=length(years))
    names(graphdf[[d]][[m]]) <- years
    graphdf[[d]][[m]] <- vector("list",length=length(years))
    names(graphdf[[d]][[m]]) <- years
  }
}

for (d in domains)
{
  for (m in models[-c(1,5,6,7)])
  {
   deltarea_df[[d]][[t]] <- do.call(rbind.data.frame, mean_area_overall[[d]][[m]])
  }
}
deltarea_df1 <- list()
for ( d in domains)

deltarea_df1[[d]] <- do.call(rbind.data.frame, deltarea_df[[d]])


for (d in domains[3])
{
  graphdf[[d]]<- deltarea_df1[[d]] %>%
  as.data.frame() %>% 
  filter(model != "0", !(anno %in% c("2009","2010"))) %>%
  ggplot(aes(x = anno, y = mean ,colour= area, group = interaction(model, area))) + 
  geom_errorbar(aes(x=anno,ymin=firt_qu, ymax=third_qu), width=.2,alpha = 0.5)+
  geom_line(aes(linetype= model)) + geom_point()+
  expand_limits(y=0) + 
  scale_x_discrete()+
  theme_bw() +
  geom_hline(yintercept=0, linetype="dashed")+
  xlab("Years")+
  ylab("Composite indicator")+
  scale_color_discrete(name  = "Macro Area",
	            guide = guide_legend()) +
  scale_linetype(name  = "Models",
                           breaks = c("1","2","3"),
                           labels = c("Marginal correlation", "CAR 3A", "CAR 3B"),
	            guide = guide_legend()) +
  theme(strip.background = element_rect(fill = NA, color = "black"),
        plot.caption = element_text(hjust = 1),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))+     
  ggtitle(paste0(d," Well-being") ,subtitle= "By macro area")
}

jpeg("Common_ENVIRONMENTAL.jpeg")
print(graphdf[["Environmental"]])
dev.off()
for (d in domains)
{
for (t in years)
   {
    graphdf[[d]][[t]] <- deltarea_df1[[d]]%>%     
      filter(anno == t)%>%
      #tidyr::gather(comune, Wellbeing) %>% 
      ggplot(aes(y= reorder(area, mean),x = mean, fill = model))+
      geom_bar(stat = "identity",alpha = 0.5, position = position_dodge())+
      geom_point(position = position_dodge(0.9), alpha = 0.8, size = 2 )+
      #facet_wrap(vars(model),scales='free')+
      geom_errorbar(aes(xmin=firt_qu, xmax=third_qu), width=.2,lty=1,alpha = 0.5, position=position_dodge(1))+
      # geom_text(aes(label= ifelse(mean %in% min(mean)
                                 # ,as.character(comune),'')),col="black",size = 5,hjust = -.05,vjust=0)+
      # geom_text(aes(label= ifelse(mean %in% c(max(mean), median(mean))
                                  
                                  # ,as.character(comune),'')),col="black",size =5, hjust = 1.08,vjust=0)+
      geom_vline(xintercept=0, linetype="dashed")+
      #scale_y_discrete(expand = expansion(mult = c(0, .02)))+
      #guides(fill=FALSE, color=FALSE)+
      xlab(paste0(t))+
      ylab("")+
      theme_bw() +
      scale_fill_discrete(name  = "Models",
                           breaks = c( "1","2","3"),
                           labels = c("Marginal correlation", "CAR 3A", "CAR 3B"),
	            guide = guide_legend()) +
      scale_shape_discrete(name  = "Models",
                           breaks = c( "1","2","3"),
                           labels = c("Marginal correlation", "CAR 3A", "CAR 3B"))+
      theme(strip.background = element_rect(fill = NA, color = "black"),
            plot.caption = element_text(hjust = 1),
            text = element_text(size=rel(5.5)),
            legend.position="bottom",
            axis.text.y=element_text(size = 10.5),
            axis.text.x=element_text(size = 10.5),
            axis.title.x = element_text(margin = unit(c(1, 0, 0, 0), "mm"),size = 14.5),
            legend.title=element_text(size=rel(4)), 
            legend.text= element_text(size=rel(4)),
            axis.ticks.y=element_blank())
  }
}

p1 <- graphdf[["Environmental"]][["2004"]]
p2 <- graphdf[["Environmental"]][["2006"]]
p3 <- graphdf[["Environmental"]][["2014"]]
p4 <- graphdf[["Environmental"]][["2017"]]
legend = gtable_filter(ggplotGrob(graphdf[["Environmental"]][["2017"]]),"guide-box")     # Make sure the legend has been extracted
jpeg("Common_ENVIRONMENTAL.jpeg",width=2500, height = 1000)
grid.arrange(arrangeGrob(
 p1 + theme(legend.position="none"),
 p2 + theme(legend.position="none"), 
 p3 + theme(legend.position="none"),
 p4 + theme(legend.position="none"), nrow =1, ncol = 4), legend,  nrow=2, heights=c(20, 5)) 
dev.off()
