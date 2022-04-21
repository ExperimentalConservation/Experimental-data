rm(list=ls())
library(tidyverse)
library(data.table)
library(openxlsx)
setwd(getwd())

#Set all time to hours, remove unnecessary columns and data we don't need
#load growth and decline curve data
d_growth<-read.csv("Prey growth curves.csv")
d_decline <-read.csv("Prey decline curves.csv")

#convert date into hours separately for each temperature
split_temps <- split(d_growth, f = d_growth$Temperature)
for(i in 1:length(split_temps)){
  split_temps[[i]][["Time"]] = as.numeric(difftime(split_temps[[i]][["Date"]], dplyr::first(split_temps[[i]][["Date"]]), units = "hours"))
}
#recombine
d_growth <- unsplit(split_temps, d_growth$Temperature)

#remove Observer and Date
d_growth <- d_growth %>%
  select(-c(Observer, Date)) %>%
  na.omit()

d_decline <- d_decline %>%
  #remove Bleph, Colp, Spi (duplicate curves in d_growth)
  filter(!Species == "Blepharisma"&!Species == "Colpidium"&!Species == "Spirostomum") %>%
  na.omit() %>%
  #species names lowercase to match d_growth
  mutate_all(.funs = tolower)

#change the time for Spirostomum 20 and 30 manually as they stared at a different time
#minus 1512hrs from Spi 20
d_growth$Time[c(which(d_growth$Time > 1511))] <- d_growth$Time[c(which(d_growth$Time > 1511))] - 1512
#minus 168hrs from Spi 30
d_growth$Time[c(which(d_growth$Species == "spirostomum" & d_growth$Temperature == 30))] <- d_growth$Time[c(which(d_growth$Species == "spirostomum" & d_growth$Temperature == 30))] - 168

#join d_growth and d_decline
dd <- rbind(d_growth, d_decline)
dd <- dd %>%
  mutate_at(c(1,3:5), as.numeric)

fin <- rbindlist(lapply(split(dd, list(dd$Species, dd$Temperature, dd$Replicate)), function(x){
  
  #x<-split(dd, list(dd$Species, dd$Temperature, dd$Replicate))[[20]]
  if(nrow(x)>0){

    x<-x[which(x$Density!=0),]
    
    dens1<-x$Density[1:nrow(x)-1]
    dens2<-x$Density[-1]

    pgr<-log(dens2/dens1)/diff(x$Time)
  
    res<-data.frame(Speces=x$Species[1],
                    Temperature=x$Temperature[1],
                    Replicate=x$Replicate[1],
                    pgr=pgr,
                    dens = dens1)

    mod1<-lm(pgr~dens, data=res)

    r <- coef(mod1)[1]
 
    pred.model<-data.frame(dens=seq(0, max(res$dens), length.out=100), 
                           pgr=predict(mod1, newdata=data.frame(dens=seq(0, max(res$dens), length.out=100))))

    if(max(pred.model$pgr)>0){

      k <- pred.model$dens[which(abs(pred.model$pgr-0)==min(abs(pred.model$pgr-0)))]
    }else{

      k<-0
    }

    return(data.frame(Species=x$Species[1],
                      Temperature=x$Temperature[1],
                      Replicate=x$Replicate[1],
                      r=r[[1]],
                      k=k))
    
  }}))
ggplot(fin, aes(x=Temperature, y=r))+geom_point(pch = 1)+facet_wrap(~Species, scales = 'free')

#remove outliers 
#bleph 35 4; Lepa 37.5 4; Colp 30 4; Paur 35; P caud 35 2; Phi 25 4
data <- fin[-c(179,184,169, 43, 90,135,181,91,166), ]

#### Fitting thermal performance curves for r and K using GAMs ####
library(mgcv)
library(ggpubr)

#### trying to run the above in a lapply loop ####
#remove spirostomum
data <- data %>%
  filter(!Species == "spirostomum")

#split data for each species
protist_preds <- rbindlist(lapply(split(data, list(data$Species)), function(x){
  #GAM for r
  gam_r <- gam(r ~ s(Temperature, k = 5, bs = "cs"), data = x, method = "REML")
  #GAM for k
  gam_k <- gam(k ~ s(Temperature, k = 5, bs = "cs"), data = x, method = "REML")
  #data of 100 Temperature vals to generate model predictions
  new_data <- data.frame(Temperature = seq(min(x$Temperature), max(x$Temperature),
                                           length.out = 100))
  #generate and return a data frame with r and k preds, temperature, species names
  return(data.frame(r = predict(gam_r, newdata = new_data),
                    k = predict(gam_k, newdata = new_data),
                    temperature = new_data,
                    species = x$Species[1]))
}))

#pivot the data frames for rs and ks
protist_preds <- pivot_longer(protist_preds, names_to = "parameter", cols = c(r, k), values_to = "parm_pred")
protist_preds <- protist_preds %>%
  rename(temperature = Temperature) #lowercase so both data frames match

plot_data <- pivot_longer(data, names_to = "parameter", cols = c(r, k), values_to = "parm_raw")
plot_data <- plot_data %>%
  rename(species = Species,
         temperature = Temperature) %>% #lowercase so both match
  select(-Replicate) %>% 
  filter(!species == 'spirostomum')

#plot r and k separately, faceted by species

#facet labels
species_names <- c("Blepharisma japonicum", "Colpidium striatum", "Lepadella sp.", "Paramecium aurelia", "Paramecium caudatum", "Philodina sp.")
names(species_names) <- c("blepharisma", "colpidium", "lepadella", "p_aurelia", "p_caudatum", "philodina")

#carrying capacity plot
K_plot <- 
ggplot(data = subset(protist_preds, parameter == "k")) +
  geom_line(aes(x = temperature, y = parm_pred), colour = '#00695C') +
  geom_point(data = subset(plot_data, parameter == "k"), aes(x = temperature, y = parm_raw), 
             pch = 1, colour = '#37474F') +
  facet_wrap(vars(species), scales = "free_y", #to vary y axis
             labeller = labeller(species = species_names),
             nrow = 1) + 
  scale_x_continuous(breaks = c(10, 15, 20, 25, 30, 35))+
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        complete = F,
        text = element_text(size = 10),
        strip.text = element_text(face = "italic"),
        axis.title.x = element_blank()) +
  ylab("Carrying capacity")

#r plot 
r_plot <- 
ggplot(data = subset(protist_preds, parameter == "r")) +
  geom_line(aes(x = temperature, y = parm_pred), colour = '#C0CA33') +
  geom_point(data = subset(plot_data, parameter == "r"), aes(x = temperature, y = parm_raw), 
             pch = 1, colour = '#37474F') +
  facet_wrap(vars(species), scales = "free_y",
             labeller = labeller(species = species_names),
             nrow = 1) + 
  scale_x_continuous(breaks = c(10, 15, 20, 25, 30, 35))+
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        complete = F,
        text = element_text(size = 10),
        strip.text = element_blank()) +
  ylab("r") +
  xlab("Temperature (ºC)") 

r_K_plot <- ggarrange(K_plot, NULL, r_plot, nrow = 3,
          heights = c(1, -0.5, 1), align = "v")

