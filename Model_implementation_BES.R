#import packages
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


