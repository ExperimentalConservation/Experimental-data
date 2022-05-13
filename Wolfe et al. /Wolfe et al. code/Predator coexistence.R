##############################
#### predator coexistence ####
##############################

#### load packages ####
library(tidyverse)
library(MuMIn)
library(pROC)
library(ResourceSelection)
library(dplyr)

#### preparing the data frame ####
#read in patch-level abundance data
abundance_data <- read.csv("final_count_patch_data.csv", header = T)

#generate the response variables
abundance_data <- abundance_data %>%
  mutate(generalist_predator = replace(stentor_coeruleus, stentor_coeruleus > 0, 1),#change Stentor abundances to presence/absence
         specialist_predators = rowSums(abundance_data[,c(13,15)]), #change specialist predator combined abundances to presence/absence
         specialist_predators = replace(specialist_predators, specialist_predators > 0, 1)) 

#create the model (specialist predator response, generalist predator predictor, but can go other way too)

glm_fit <- glm(specialist_predators ~ generalist_predator, data = abundance_data, family=binomial)
summary(glm_fit)

#model diagnostics
#add preds to dataframe
abundance_data$pred_binomial <- predict(glm_fit, type = "response")

#Hosmer-Lemeshow goodness of fit test
hl <- hoslem.test(abundance_data$specialist_predators, abundance_data$pred_binomial, g = 7)
hl #p value non-significant so model is not a bad fit 

#### plotting model output ####
final_predator_plot <- abundance_data %>%
  mutate(prob = ifelse(specialist_predators == "1", 1, 0)) %>%
  ggplot(aes(generalist_predator, specialist_predators))+
  geom_jitter(height= 0.03, width = 0.03, colour = "#2D708EFF", pch = 1, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), colour = "#2D708EFF")+
  labs(x = "Generalist predator presence",
       y = "Specialist predator presence")+
  theme_classic()+
  theme(text = element_text(size = 12),
        legend.position = "top",
        aspect.ratio = 1,
        complete = F)  

ggsave("Figure_5.tiff", final_predator_plot, device = "tiff", units = "in", width = 5, height = 5)
