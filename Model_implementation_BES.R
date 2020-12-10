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
#English names
names <- vector("list",length=4)
names(names) <- c("Overall","Social","Economic","Environmental")
names[[1]] <-  c("Employees in cultural business",
                 "Prison density", "Other reported crimes","Youth (<40 years old) political representation in municipalities",
                 "Women' s political representation in municipalities","Collection capacity: provincial governments",
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
                 "Rate of bank' s non-performing loans to households","Non-participation rate",
                 "Youth non-participation rate (15-29 years)","Employment rate(20-64 years)",
                 "Youth employment rate (15-29 years)"
)%>% janitor::make_clean_names()
names[[2]] <- c(  "Life expectancy at birth","Infant mortality rate",
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
                   "Women' s political representation in municipalities",
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
                   "Density of historical green areas")%>% janitor::make_clean_names()
names[[3]] <- c("Average disposable per capita income",
                   "Average yearly earnings of employee",
                   "Average yearly per capita pension income",
                   "Pensioners with low-pension",
                   "Per capita asset",
                   "Rate of bank' s non-performing loans to households",
                   "Employment rate(20-64 years)",
                   "Non-participation rate",
                   "Youth non-participation rate (15-29 years)",
                   "Youth employment rate (15-29 years)",
                   "Working days of paid (employee)",
                   "Mortal accidents and inabilities rate"
                   
) %>% janitor::make_clean_names()
names[[4]] <- c("Landfill of urban waste",
                  "Availability of urban green areas",
                  "Energy from renewable sources",
                  "Separate collection of municipal waste"
)%>%  janitor::make_clean_names()


#Non normalized data
colnames(data)[2:46] <- names[[1]]
colnames(data_norm)[1:45] <- names[[1]]
#Remove absent provinces
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
italgeo <- italgeo[order(italgeo$comune),]

#Spatial correlation matrices---------------------------------------------------
#Model 2
longlat04 <- data %>% 
  filter(anno=="2004") %>% 
  dplyr::select(`data$long`,`data$lat`) %>% as.matrix()
PSI_04 <- as.matrix(dist(longlat04))#Euclidean Distance
PSI_04.inv <- 1/PSI_04
diag(PSI_04.inv) <- 0

#Model CAR 2-3
knn1 <- knearneigh(longlat04)
k1 <- knn2nb(knn1)
critical.threshold <- max(unlist(nbdists(k1,longlat04)))
nb.dist.band <- dnearneigh(longlat04, 0, critical.threshold)
summary(nb.dist.band)
plot(nb.dist.band, longlat04, lwd=.2, col="blue", cex = .5)
R <- nb2mat(nb.dist.band, style='B',zero.policy=TRUE)
N <- dim(R)[1]

#Create the matrix for model 3B
o_b<- rep(NA,N)
Xi <- eigen(R)
xi <- sort(Xi$values)
xi_inv <- 1/(xi)
#if (xi_inv[1] > x  | x > xi_inv[N])
#  x=0.12
for (i in 1:N){
  o_b[i]=length(nb.dist.band[[i]])
}
O_b=diag(o_b)

W_b=matrix(0,N,N)
for (i in 1:N)
{
  for (j in 1:N)
  {
    if (j %in% nb.dist.band[[i]])
      W_b[i,j] = sqrt(o_b[i]*o_b[j])
    else
      W_b[i,j]= 0
  }
}
#Model 3D (Not implemented)
		     
PSI_04_d <- nbdists(w, coords= longlat04,longlat = TRUE)  #Create the list of neighbors

#Model implementation-------------------------------------------------------
models <-  c("0","1","2","3")
years <- c("2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015","2016","2017")
domains <- c("Social","Economic","Environmental")
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

load("ord_name.RData")
y <- vector("list",length=4)
y_norm <- vector("list",length=4)
y_center <- vector("list",length=4)
names(y) <- domains
names(y_norm) <- domains
names(y_center) <- domains
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
       dplyr::select(-c(comune,anno,`data$long`,`data$lat`)) %>% 
       as.matrix()
     y_norm[[d]][[t]] <-data_norm %>%
       tidyr::gather(Indicatore, Value,-c(anno,comune)) %>% 
       filter(Indicatore %in% names[[d]]) %>% 
       tidyr::spread(Indicatore,Value) %>%     
       dplyr::filter(anno==t) %>% 
       dplyr::select(-c(anno,comune)) %>% 
       setcolorder(.,as.character(ord_name[[d]][["2017"]])) %>%    
       as.matrix()
     y_center[[d]][[t]] <- data %>%
       tidyr::gather(Indicatore, Value,-c(anno,comune,`data$long`,`data$lat`)) %>% 
       filter(Indicatore %in% names[[d]]) %>% 
       tidyr::spread(Indicatore,Value) %>% 
       dplyr::filter(anno == t) %>% 
       dplyr::select(-c(anno,comune,`data$long`,`data$lat`)) %>% 
       as.matrix() %>% 
       scale(.,center = TRUE,scale=apply(., 2, sd, na.rm = TRUE))

  }
}

source("Sampler_debug.R")
Seed=123
for (d in domains[c(1,2,3)])
{
  for (m in models[models %in% c("0","1","2","3")])
  {
    for (t in years[years %in% c("2008","2011")])
    {
      N=dim(y_norm[[d]][[t]])[1]
      J=dim(y_norm[[d]][[t]])[2]
      L=1
      inits <- list(lambda=matrix(10,J,L),
                    mu=matrix(0,J,1),
                    Sigma=diag(1,J),
                    PSI=diag(N),
                    delta= rmvnorm(N,rep(0),diag(1)),
                    a= ifelse (m %in% c("1","4"),1,0)
      )
      set.seed(Seed + as.numeric(m) + ifelse(m %in% c("3","5"),2,1))
	  xx <- sprintf("Data/output4_%s_%s_%s.RData",d,m,t)
	  print(paste(d, m, t))
      output <- simple(n_iter = 30000,burn_in = 20000,
                                      y =y_norm[[d]][[t]],
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
                                      tuning=ifelse(m == "1",sqrt(1000), sqrt(10/N)),
		        tuning_lambda= ifelse(d == "Social",20,ifelse(d == "Economic",10,5)),
                                      model = m)
	save(output,file=xx);
   
    }
  }
}
}
