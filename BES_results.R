# install.packages("spdep")
# install.packages("knitr")
# install.packages("bestNormalize")
# install.packages("dlm")
# install.packages("xtable")
#install.packages("rstan")
# install.packages("tilingArray")
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
library(data.table)
library(grid)
library(gridExtra)
library(janitor)
library(rstan)
library("tilingArray")
rm(list=ls())

#Import the data.

#This data must include: the vector of observations, normalized or not normalized,
#the geo location of your obervations


load("data.RData")

#Diagnostic convergence (Year 2017)--------------------------------------------------------------

#Factor loading ( First factor )
dev.off()
par(mfrow=c(2,1))
par(mar=c(1,1,1,1))

for (d in domains[c(1,2,3)])
{ 
  for (m in models[models %in% c("1")])
  {
    for (t in years[years %in% c("2017")])
    {
      xx <- sprintf("Data/output4_%s_%s_%s.RData",d,m,t)
      load(file=xx);
      for (l in c(1))
      {  

 
        png(xf, width = 700, height = 350) 
        plot(output$lambda[,l], type="l",xlab = "iterations", ylab=  paste(d, m, t, l))
        plot(cumsum((output$lambda[,l]))/(1:10000), type="l", xlab = "iterations", ylab=  paste(d, m, t, l))
        dev.off()
      }
   }
 } 
}

#Latent variable (Year 2014)
dev.off()       
par(mfrow=c(2,1))
par(mar=c(1,1,1,1))
for (d in domains)
{
  for (m in models[models %in% c("0","1","2","3")])
  {
    for (t in years[years %in% c("2014")])
    {
	   xx <- sprintf("Data/output4_%s_%s_%s.RData",d,m,t)
	   load(file=xx)
      for (l in c(1,2,3))
      { 
        plot(output$delta[,l], type="l", xlab="iterations", ylab= paste(d, m, t, l))
        plot(cumsum(output$delta[,l])/(1:10000), type="l", xlab="iterations", ylab= paste(d, m, t, l)) 
        readline(paste(d, m, t, l))           
      }      
    } 
  } 
}

#Mu
dev.off()

par(mfrow=c(2,1))
par("mar")
par(mar=c(1,1,1,1))

for (d in domains[1])
{
  for (m in models[models %in% c("0","1","2","3")])
  {
    for (t in years[years %in% c("2004","2006","2014","2017")])
    {
	   xx <- sprintf("Data/output3_%s_%s_%s.RData",d,m,t)
	   load(file=xx);
      for (l in 1:length(names[[d]]))
      { 
        plot(output$mu[,l], type="l")
        plot(cumsum((output$mu[,l]))/(1:10000), type="l")
        readline(paste(d, m, t, l))
      }
      
    } 
  } 
}

#Sigma

dev.off()       
par(mfrow=c(2,1))
par(mar=c(1,1,1,1))
for (d in domains[3])

{
  for (m in models[models %in% c("3")])
  {
    for (t in years[years %in% c("2004")])
    {
	   xx <- sprintf("Data/output3_%s_%s_%s.RData",d,m,t)
	   load(file=xx)
      for (l in 1:length(names[[d]]))
      { 
        plot(output$Sigma[,l], type="l", xlab="iterations", ylab= paste(d, m, t, l))
        plot(cumsum(output$Sigma[,l])/(1:5000), type="l", xlab="iterations", ylab= paste(d, m, t, l)) 
        readline(paste(d, m, t, l))           
      }      
    } 
  } 
}
#Creation of empty objects
##Model selection----------------------------------------------------------------

y_rep <- vector("list",length=length(domains))
sd <-  vector("list",length=length(domains))
G <- vector("list",length=length(domains))
p <- vector("list",length=length(domains))
C <- vector("list",length=length(domains))
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

###Factor loadings----------------------------------------------------------------
summary <- vector("list",length=length(domains))
sigma <- vector("list",length=length(domains))
f <- vector("list",length=length(domains))

#Spatial parameter
a <- vector("list",length=length(domains))

#Composite indicator
delta <- vector("list",length=length(domains))
d_summary <- vector("list",length=length(domains))
#Naming
names(summary) <- domains
names(sd) <- domains
names(y_rep) <- domains
names(G)<- domains
names(p) <- domains
names(C) <- domains
names(sigma) <- domains
names(f) <- domains
names(delta) <- domains
names(d_summary) <- domains


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
  summary[[d]] <- vector("list",length=length(models))
  names(summary[[d]]) <- models
  sigma[[d]] <- vector("list",length=length(models))
  names(sigma[[d]]) <- models
  f[[d]] <- vector("list",length=length(models))
  names(f[[d]]) <- models
  a[[d]] <- vector("list",length=length(models))
  names(a[[d]]) <- models
  delta[[d]] <- vector("list",length=length(models[-c(5,6,7)]))
  names(delta[[d]]) <- models
  d_summary[[d]] <- vector("list",length=length(models[-c(5,6,7)]))
  names(d_summary[[d]]) <- models
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
     summary[[d]][[m]] <- vector("list",length=length(years))
    names(summary[[d]][[m]]) <- years
    sigma[[d]][[m]] <- vector("list",length=length(years))
    names(sigma[[d]][[m]]) <- years
    f[[d]][[m]] <- vector("list",length=length(years))
    names(f[[d]][[m]]) <- years
    a[[d]][[m]] <- vector("list",length=length(years))
    names(a[[d]][[m]]) <- years
    delta[[d]][[m]] <- vector("list",length=length(years))
    names(delta[[d]][[m]]) <- years
    d_summary[[d]][[m]] <- vector("list",length=length(years))
    names(d_summary[[d]][[m]]) <- years
  }
}
for (d in domains[c(1,2,3)])
{
  for (m in models[models %in% c("0","1","2","3")])
  {
    for (t in years[years %in% c("2004","2005","2006","2008","2010","2011","2012","2014","2016","2017")])
      {
          xx <- sprintf("Data/output4_%s_%s_%s.RData",d,m,t)
          load(file=xx)
          N=dim(y_norm[[d]][[t]])[1]
          J=dim(y_norm[[d]][[t]])[2]
          y_rep[[d]][[m]][[t]] <- apply(output$y_rep[-c(1:2000),],2,mean) %>% matrix(N,J,byrow=TRUE)
          sd[[d]][[m]][[t]] <- apply(output$y_rep[-c(1:2000),],2,var) %>% matrix(N,J,byrow=TRUE)
          G[[d]][[m]][[t]] <- sum(rowSums((y_rep[[d]][[m]][[t]]- y_norm[[d]][[t]])^2))
          p[[d]][[m]][[t]] <- sum(rowSums(sd[[d]][[m]][[t]]))
          C[[d]][[m]][[t]] <- G[[d]][[m]][[t]] + p[[d]][[m]][[t]]
          Sigma <- output$Sigma[-c(1:2000),] %>% as.data.frame()   
          lambda <- output$lambda[-c(1:2000),] %>% as.data.frame()
          foutput <- (lambda ^ 2) /(lambda^2+Sigma) %>% as.data.frame()
         f[[d]][[m]][[t]] <- rbind(apply(foutput, 2,mean),apply(foutput,2, quantile,c(.025, .50, .975)),apply(foutput,2, IQR))
         colnames(f[[d]][[m]][[t]]) <- ord_name[[d]][["2017"]]
         f[[d]][[m]][[t]] <- f[[d]][[m]][[t]] %>%
                                             as.data.frame() %>%
			mutate(
			anno=t, 
			model = m,
			row_names = c("mean", "first","median","third","IQR")) %>%       
			as.list()         
          summary[[d]][[m]][[t]] <- rbind(apply(lambda, 2,mean),apply(lambda,2, quantile,c(.025, .50, .975)),apply(lambda,2, IQR)) 
          colnames(summary[[d]][[m]][[t]]) <- ord_name[[d]][["2017"]]      
          summary[[d]][[m]][[t]] <-  summary[[d]][[m]][[t]] %>%
                                             as.data.frame() %>%
			mutate(
			anno=t, 
			model = m,
			row_names = c("mean", "first","median","third","IQR")) %>%       
			as.list()  
          
          sigma[[d]][[m]][[t]] <- rbind(apply(Sigma, 2,mean),apply(Sigma,2, quantile,c(.025, .50, .975)),apply(Sigma,2, IQR))
          colnames(sigma[[d]][[m]][[t]]) <-  ord_name[[d]][["2017"]]
          a[[d]][[m]][[t]] <- output$a[-c(1:2000)] %>% 
			as.data.frame() %>% 
			dplyr::summarise(min=min(.),first_qu = quantile(.,0.25),
			mean=mean(.),median=median(.),third_qu=quantile(.,0.75),
			max=max(.),sd=sd(.),IQR=IQR(.),year=t,model= m,domain = d) %>% 
			as.list()	          
          delta[[d]][[m]][[t]] <- output$delta[-c(1:2000),] %>% as.data.frame()
          d_summary[[d]][[m]][[t]] <- multi.sapply(delta[[d]][[m]][[t]], mean, median , third_qu=function(x) quantile(x, 0.75),firt_qu= function(x) quantile(x,0.25), IQR = IQR) %>%
			as.data.frame()%>%
			mutate(
			comune = province,
			anno=t ,
			model = m
			) %>% 
			as.list()
          
     }                  
  }
}


#Model selection
kable(G[["Social"]][["0"]][["2004"]], format="latex", digit=2)
kable(p[["Social"]][["0"]][["2004"]], format="latex", digit=2)
kable(C[["Environmental"]][["0"]][["2004"]], format="latex", digit=2)
r <- matrix(c(G[[1]][["1"]][["2004"]],G[[1]][["1"]][["2006"]],G[[1]][["1"]][["2014"]],G[[1]][["1"]][["2017"]],
               p[[1]][["1"]][["2004"]],p[[1]][["1"]][["2006"]],p[[1]][["1"]][["2014"]],p[[1]][["1"]][["2017"]],
	C[[1]][["1"]][["2004"]],C[[1]][["1"]][["2006"]],C[[1]][["1"]][["2014"]],C[[1]][["1"]][["2017"]]
	),nrow=4, ncol=3)
colnames(r) <- c("G(m)","p(m)","C(m)")
rownames(r) <- c("2004","2006","2014","2017")
kable(t(r), format="latex", digit=2)

#FACTOR LOADINGS TABLE

kable(t(summary[["Economic"]][["0"]][["2005"]]), format="latex",digit=2)
summary[["Social"]][["1"]][["2005"]]
summary[["Economic"]][["1"]][["2004"]]
summary[["Environmental"]][["1"]][["2006"]]


#PLOT FACTOR LOADINGS
plot_summary <- vector("list",length=length(domains))
names(plot_summary) <- domains
summary_df <- vector("list",length=length(domains))
names(summary_df) <- domains
for (d in domains)
{
  plot_summary[[d]] <- vector("list",length=length(models[-c(5,6,7)]))
  names(plot_summary[[d]]) <- models[-c(5,6,7)]
  summary_df[[d]] <- vector("list",length=length(models[-c(5,6,7)]))
  names(summary_df[[d]]) <- models[-c(5,6,7)]
  for (m in models[-c(5,6,7)])
  {
    plot_summary[[d]][[m]] <- vector("list",length=length(years))
    names(plot_summary[[d]][[m]]) <- years
  }
}

for (d in domains)
{
  for (m in models[-c(5,6,7)])
  {
   summary_df[[d]][[m]] <- do.call(rbind.data.frame, summary[[d]][[m]])
  }
}
summary_df_d <- vector("list",length=length(domains))
names(f_df_d) <- domains

for (d in domains)
{
  summary_df_d[[d]] <- do.call(rbind.data.frame, summary_df[[d]])
}
width_scale <- 12
for (d in domains)
{
  {
   for (t in years[years %in% c("2004","2005","2006","2008","2010","2011","2012","2014","2016","2017")])
   {
    plot_summary[[d]][[t]] <- summary_df_d[[d]]%>% 
    clean_names() %>%     
    filter(anno == t)%>%
    tidyr::gather(Indicatore, value,-c(anno,model,row_names)) %>%
    tidyr::spread(row_names, value) %>%
    mutate(xmin = pmin(first,third),xmax = pmax(first,third)) %>%
      ggplot(aes(y= reorder(Indicatore, mean),x = mean, group = model))+
      geom_bar(aes(fill = model), stat = "identity",alpha = 0.5,width=0.75, position ='dodge')+
      geom_point(position = position_dodge(0.75), alpha = 0.8, size = 2 )+
      geom_errorbar(aes(xmin=xmin, xmax= xmax), lty=1, alpha = 0.5,width = 0.2, position=position_dodge(0.75))+
      # geom_text(aes(label= ifelse(mean %in% min(mean)
                                 # ,as.character(comune),'')),col="black",size = 5,hjust = -.05,vjust=0)+
      # geom_text(aes(label= ifelse(mean %in% c(max(mean), median(mean))
                                  
                                  # ,as.character(comune),'')),col="black",size =5, hjust = 1.08,vjust=0)+
      geom_vline(xintercept=0, linetype="dashed")+
      #scale_y_discrete(expand = expansion(mult = c(0, .02)))+
      #guides(fill=FALSE, color=FALSE)+
      xlab(paste0(t,"\n Squared correlations \n",d))+
      ylab("")+
      theme_light() +
      scale_fill_discrete(name  = "Models",
                           breaks = c("0", "1","2","3"),
                           labels = c("Spatial Independence","Marginal correlation", "CAR 3A", "CAR 3B"),
	            guide = guide_legend()) +
      scale_shape_discrete(name  = "Models",
                           breaks = c("0", "1","2","3"),
                           labels = c("Spatial Independence","Marginal correlation", "CAR 3A", "CAR 3B"))+
     theme(strip.background = element_rect(fill = NA, color = "black"),
            plot.caption = element_text(hjust = 1),
            text = element_text(size=rel(6.5)),
            legend.position="bottom",
            axis.text.y=element_text(size = 35),
            axis.text.x=element_text(size = 25),
            axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
            axis.ticks.y=element_blank(),
            legend.title=element_text(size=rel(10.5)), 
            legend.text= element_text(size=rel(10.5)))
   }
 }
}

#SQUARED CORRELATIONS---------------------------------------------------------------
#Plot square correlations
plot_f <- vector("list",length=length(domains))
names(plot_f) <- domains
f_df <- vector("list",length=length(domains))
names(f_df) <- domains
for (d in domains)
{
  plot_f[[d]] <- vector("list",length=length(models[-c(5,6,7)]))
  names(plot_f[[d]]) <- models[-c(5,6,7)]
  f_df[[d]] <- vector("list",length=length(models[-c(5,6,7)]))
  names(f_df[[d]]) <- models[-c(5,6,7)]
  for (m in models[-c(5,6,7)])
  {
    plot_f[[d]][[m]] <- vector("list",length=length(years))
    names(plot_f[[d]][[m]]) <- years
  }
}

for (d in domains)
{
  for (m in models[-c(5,6,7)])
  {
   f_df[[d]][[m]] <- do.call(rbind.data.frame, f[[d]][[m]])
  }
}
f_df_d <- vector("list",length=length(domains))
names(f_df_d) <- domains

for (d in domains)
{
  f_df_d[[d]] <- do.call(rbind.data.frame, f_df[[d]])
}
width_scale <- 12
for (d in domains)
{
  {
   for (t in years[years %in% c("2004","2005","2006","2008","2010","2011","2012","2014","2016","2017")])
   {
    plot_f[[d]][[t]] <- f_df_d[[d]]%>% 
    clean_names() %>%     
    filter(anno == t)%>%
    tidyr::gather(Indicatore, value,-c(anno,model,row_names)) %>%
    tidyr::spread(row_names, value) %>%
    mutate(xmin = pmin(first,third),xmax = pmax(first,third)) %>%
      ggplot(aes(y= reorder(Indicatore, mean),x = mean, group = model))+
      geom_bar(aes(fill = model), stat = "identity",alpha = 0.5,width=0.75, position ='dodge')+
      geom_point(position = position_dodge(0.75), alpha = 0.8, size = 2 )+
      geom_errorbar(aes(xmin=xmin, xmax= xmax), lty=1, alpha = 0.5,width = 0.2, position=position_dodge(0.75))+
      # geom_text(aes(label= ifelse(mean %in% min(mean)
                                 # ,as.character(comune),'')),col="black",size = 5,hjust = -.05,vjust=0)+
      # geom_text(aes(label= ifelse(mean %in% c(max(mean), median(mean))
                                  
                                  # ,as.character(comune),'')),col="black",size =5, hjust = 1.08,vjust=0)+
      geom_vline(xintercept=0, linetype="dashed")+
      #scale_y_discrete(expand = expansion(mult = c(0, .02)))+
      #guides(fill=FALSE, color=FALSE)+
      xlab(paste0(t,"\n Squared correlations \n",d))+
      ylab("")+
      theme_light() +
      scale_fill_discrete(name  = "Models",
                           breaks = c("0", "1","2","3"),
                           labels = c("Spatial Independence","Marginal correlation", "CAR 3A", "CAR 3B"),
	            guide = guide_legend()) +
      scale_shape_discrete(name  = "Models",
                           breaks = c("0", "1","2","3"),
                           labels = c("Spatial Independence","Marginal correlation", "CAR 3A", "CAR 3B"))+
     theme(strip.background = element_rect(fill = NA, color = "black"),
            plot.caption = element_text(hjust = 1),
            text = element_text(size=rel(6.5)),
            legend.position="bottom",
            axis.text.y=element_text(size = 35),
            axis.text.x=element_text(size = 25),
            axis.title.x = element_text(margin = unit(c(4, 0, 0, 0), "mm")),
            axis.ticks.y=element_blank(),
            legend.title=element_text(size=rel(10.5)), 
            legend.text= element_text(size=rel(10.5)))
   }
 }
}

#SPATIAL PARAMETER--------------------------------------------------------------
#DATA FRAME SPATIAL PARAMETER
a_df <- vector("list",length=length(domains))
a_df_d <- vector("list",length=length(domains))
names(a_df) <- domains
names(a_df_d) <- domains
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
for (d in domains)
{
  a_df_d[[d]] <- do.call(rbind.data.frame, a_df[[d]])
}
#PLOT
a_df_d_tot <- do.call(rbind.data.frame, a_df_d)
a_df_d_tot$model <- plyr::revalue(x = a_df_d_tot$model, 
c("1" = "Marginal correlation", "2" = "CAR 3A", "3" = "CAR 3B"))

pd <- position_dodge(0.1) 

p <- ggplot(a_df_d_tot %>% filter(model %in% c("Marginal correlation","CAR 3A","CAR 3B")) ,aes(y=mean,x=year,col = domain))+
    geom_errorbar(aes(ymin=first_qu, ymax=third_qu), width=.4,lty=1, position=pd) +
    geom_point(size=1)+
    facet_wrap(. ~ model, scales="free_y",ncol= 1)+
    geom_line(aes(colour = domain, linetype = domain,group = domain)) +
    geom_hline(yintercept=0, linetype="dashed")+
    ylab(expression(omega))+
    theme_bw() +
    theme(strip.background = element_rect(fill = NA, color = "black"),
          axis.text.x = element_text(size = 8, angle = 45,vjust=0.5),
          plot.caption = element_text(hjust = 1))
png("plotspatial2.png")
print(p)     
dev.off() 
#COMPOSITE INDICATORS----------------------------------------------------------------

#PLOT--------------------------------------------------------------------------------

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

for (d in domains)
{
  {
   for (t in years[years %in% c("2004","2005","2006","2008","2010","2011","2012","2014","2016","2017")])
   {
    plot_delta[[d]][[t]] <- delta_df_d[[d]]%>%     
    filter(anno == t, comune %in% c("Milano","Roma","Torino","Genova","Palermo","Cagliari",
                                                "Trento","Venezia","Perugia","Potenza","Bari",
                                              "Napoli","Bologna","Firenze","Aosta","Reggio di Calabria",
                                              "Ancona","Trieste","Campobasso","Catanzaro"))%>%
      ggplot(aes(y= reorder(comune, mean),x = mean, fill = model))+
      geom_bar(stat = "identity",alpha = 0.5, position = position_dodge())+
      geom_point(position = position_dodge(0.9), alpha = 0.8, size = 2 )+
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
##MAPS----------------------------------------------------------------------------------------------------------
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
  }
}

for (d in domains)
{ 
  for (m in models[models %in% c("0","1","2","3")])
  {
  for (t in years[years %in% c("2004","2005","2006","2008","2010","2011","2012","2014","2016","2017")])
  {
  VARmap[[d]][[m]][[t]] <- merge(italgeo,delta_df_d[[d]] %>% filter(model == m ,anno == t),by = "comune")
  }
 }
}
for (d in domains)
{
 for (m in models[models %in% c("0","1","2","3")])
  {
  VARmap[[d]][[m]] <- rbind(VARmap[[d]][[m]][["2004"]],VARmap[[d]][[m]][["2006"]],VARmap[[d]][[m]][["2008"]],VARmap[[d]][[m]][["2011"]],VARmap[[d]][[m]][["2014"]],VARmap[[d]][[m]][["2017"]])
  map[[d]][[m]] <- tm_shape(VARmap[[d]][[m]])+
  tm_fill("mean",title= paste0(d,m), style = "quantile", midpoint = NA)+
  tm_borders(alpha=.5)+
  tm_facets(by=c("anno"), ncol  =2, showNA = FALSE)+
  tm_shape(VARmap[[d]][[m]] %>%  filter(comune %in% c("Milano","Roma","Torino","Genova","Palermo","Cagliari",
                                                 "Trento","Venezia","Perugia","Potenza","Bari",
                                                 "Napoli","Bologna","Firenze","Aosta","Reggio di Calabria",
                                                 "Ancona","Trieste","Campobasso")))+
  tm_dots("comune",size=0.2, col="black",alpha=0.5, shape= 21, legend.show = FALSE)
  }
}


