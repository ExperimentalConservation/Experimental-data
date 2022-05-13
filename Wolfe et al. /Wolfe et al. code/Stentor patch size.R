###############################################################
#### effects of patch size on generalist predator presence ####
###############################################################

#load required packages
library(tidyverse)
library(ResourceSelection) 
library(pROC)
library(dplyr)

stentor_data <- read.csv("stentor_patch_size.csv", header = T)

stentor_data <- stentor_data %>%
  mutate(stentor_presence = as.numeric(replace(stentor_abundance, stentor_abundance > 0, 1))) %>%
  filter(patch_size <= 16)

#run the model
model1 <- glm(stentor_presence ~ patch_size,
              family="binomial", 
              data = stentor_data)
summary(model1)

#significant positive effect of patch size on Stentor
#return predicted values from model 
stentor_data$pred_binomial <- predict(model1, type = "response")

#model diagnostics
hoslem.test(stentor_data$stentor_presence, stentor_data$pred_binomial,
            g = 15)
#p-value non significant so no evidence of a bad fit

#area under the ROC curve
#generate a roc object with response and fitted values
roccurve <- roc(stentor_data$stentor_presence, stentor_data$pred_binomial)
#estimate area under the ROC curve
auc(roccurve)# = 0.6814 so an okay fit

subset_plot <- ggplot(data = stentor_data) +
  geom_jitter(mapping = aes(x = patch_size, y = stentor_presence), 
              height = 0.03, width = 0.8, colour = "#2D708EFF", pch = 1, alpha = 0.5) +
  geom_smooth(mapping = aes(x = patch_size, y = stentor_presence),
              method = "glm",
              method.args = list(family = "binomial"),
              colour = "#2D708EFF") +
  xlab("Patch size (mL)")+
  ylab("Probability of generalist predator presence")+
  scale_x_continuous(breaks = c(4, 8, 12, 16), limits = c(0, 20))+
  theme(strip.background = element_rect(colour = "white", fill = "white",
                                        linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        text = element_text(size = 12),
        aspect.ratio = 1,
        complete = F)


