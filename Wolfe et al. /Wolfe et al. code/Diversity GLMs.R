########################
#### Diversity GLMs ####
########################

#load required packages
library(dplyr)
library(tidyverse)
library(MuMIn)
library(scales)

#### Read the data in ####
microcosm_data <- read.csv("diversity_data.csv", header = T)
microcosm_data$X <- NULL
#Re-order the factors so the 'control' comes first (not alphabetical)
microcosm_data <- microcosm_data %>%
  mutate_if(is.character, as.factor)%>%
  mutate(matrix_dispersal = relevel(matrix_dispersal, "none"),
         corridor_dispersal = relevel(corridor_dispersal, "none"),
         heterogeneity = relevel(heterogeneity, "homogeneous"),
         #add ordered matrix and corridor dispersal for plotting facets later
         ordered_matrix_dispersal = factor(matrix_dispersal, levels = c("none", "low", "high")),
         ordered_corridor_dispersal = factor(corridor_dispersal, levels = c("none", "low", "high")))

beta_data <- microcosm_data %>%
  filter(!microcosm_type == "SL") %>% #remove SL and 6HoLC3 microcosms
  filter(!microcosm_ID == "6HoLC3") %>%
  dplyr::select(-c(alpha_diversity, gamma_diversity)) %>%
  mutate(beta_diversity = as.numeric(as.character(beta_diversity)))

#### Alpha Diversity ####
#create saturated model
alpha_glm <- glm(alpha_diversity ~ (heterogeneity + patch_number + corridor_dispersal + matrix_dispersal)^2,
                 data = microcosm_data, family = "gaussian", na.action="na.fail")

summary(alpha_glm)
#visual check of model fit
par(mfrow = c(2, 2))
plot(alpha_glm)

#model simplification using AIC
dredge_alpha_glm <- dredge(alpha_glm, evaluate = TRUE)

best_alpha_models <- subset(dredge_alpha_glm, delta < 2)
#just one top model so don't need to average
best_alpha_glm <- glm(alpha_diversity ~ heterogeneity + patch_number + corridor_dispersal
                      + matrix_dispersal + matrix_dispersal:patch_number, family = "gaussian",
                      data = microcosm_data, na.action = "na.fail")
plot(best_alpha_glm)
summary(best_alpha_glm) 

#extract se to make a final plot
alpha_ilink <- family(best_alpha_glm)$linkinv
                                          
#create a dataframe with unique variables
new_alpha_data <- microcosm_data[,c(3,5,6,7)]

#add fit and se.fit on the link scale
new_alpha_data <- bind_cols(new_alpha_data,
                            setNames(as_tibble(predict(best_alpha_glm, new_alpha_data, se.fit = TRUE)[1:2]),
                                     c('alpha_fit_link', 'alpha_se_link')))
#create the interval and backtransform
new_alpha_data <- new_alpha_data %>%
  mutate(alpha_fit_resp = alpha_ilink(alpha_fit_link),
         alpha_right_upr = alpha_ilink(alpha_fit_link + (1.96 * alpha_se_link)),
         alpha_right_lwr = alpha_ilink(alpha_fit_link - (1.96 * alpha_se_link)))

#add upr and lwr to microcosm data
microcosm_data <- microcosm_data %>%
  mutate(alpha_preds = new_alpha_data$alpha_fit_resp,
         alpha_right_upr = new_alpha_data$alpha_right_upr,
         alpha_right_lwr = new_alpha_data$alpha_right_lwr)

#change the facet labels
matrix.labs <- c("No dispersal", "Low mortality", "High mortality")
names(matrix.labs) <- c("none", "low" ,"high")
corridor.labs <- c("No dispersal", "Low frequency", "High frequency")
names(corridor.labs) <-c("none", "low", "high")

#plot it all together
alpha_final_plot <- ggplot(data = microcosm_data)+
  geom_line(mapping = aes(x = patch_number, y = alpha_preds, colour = heterogeneity))+
  geom_ribbon(mapping = aes(x = patch_number,
                            ymin = alpha_right_lwr,
                            ymax = alpha_right_upr,
                            group = heterogeneity),
              alpha = 0.1) +
  geom_jitter(aes(x = patch_number, y = alpha_diversity, colour = heterogeneity),
              width = 0.2, height = 0.2)+
  xlab("Number of patches") +
  ylab("Alpha diversity") +
  labs(colour = "Heterogeneity") +
  scale_colour_viridis_d()+
  facet_grid(ordered_matrix_dispersal~ordered_corridor_dispersal,
             labeller=ggplot2::labeller(ordered_matrix_dispersal = matrix.labs,
                                        ordered_corridor_dispersal = corridor.labs))+
  scale_x_continuous(breaks = c(1, 4, 6))+
  scale_colour_manual(values = c("homogeneous" = "#B8DE29FF", "heterogeneous" = "#2D708EFF"))+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        aspect.ratio = 1,
        text = element_text(size = 12),
        legend.position = "top")+
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 1)+
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 1 )


#### Beta diversity ####
#create saturated model
beta_glm <- glm(beta_diversity ~ (heterogeneity + patch_number + corridor_dispersal + matrix_dispersal)^2,
                data = beta_data, family = "gaussian", na.action="na.fail")

summary(beta_glm)
#visually inspect model fit
par(mfrow = c(2, 2))
plot(beta_glm)

#model selection using AIC
dredge_beta_glm <- dredge(beta_glm, evaluate = TRUE)
best_beta_models <- subset(dredge_beta_glm, delta < 2)

# model averaging 
beta_averages <- model.avg(best_beta_models, fit = TRUE)
#extract new model preds and add to microcosm_data
summary(beta_averages)

#create a dataframe with unique variables
new_beta_data <- beta_data[,c(3,5,6,7)]

#beta_data$avg_beta_preds <- predict(beta_averages, new_beta_data, full = TRUE, type = "response")

#add fitted values and se.fit on the link scale
new_beta_data <- bind_cols(new_beta_data, 
                           setNames(as_tibble(predict(beta_averages,
                                                      new_beta_data, se.fit = TRUE)[1:2]),
                                    c('average_beta_fit_link', 'average_beta_se_link')))

#extract family and create inverse link
ilink <- family(beta_glm)$linkinv

new_beta_data <- new_beta_data %>%
  mutate(average_beta_fit_resp = ilink(average_beta_fit_link),
         average_beta_right_upr = ilink(average_beta_fit_link + (1.96 * average_beta_se_link)),
         average_beta_right_lwr = ilink(average_beta_fit_link - (1.96 * average_beta_se_link)))

#add right_upr and right_lwr to microcosm data - makes plotting easier
beta_data <- beta_data %>%
  mutate(beta_fit_resp = new_beta_data$average_beta_fit_resp,
         average_beta_right_lwr = new_beta_data$average_beta_right_lwr,
         average_beta_right_upr = new_beta_data$average_beta_right_upr)

#plot it all together
average_beta_final_plot <- ggplot(data = beta_data)+
  geom_line(mapping = aes(x = patch_number, y = beta_fit_resp, colour = heterogeneity))+
  geom_ribbon(mapping = aes(x = patch_number,
                            ymin = average_beta_right_lwr,
                            ymax = average_beta_right_upr,
                            group = heterogeneity),
              alpha = 0.1) +
  geom_jitter(aes(x = patch_number, y = beta_diversity, colour = heterogeneity),
              width = 0.2, height = 0.2)+
  xlab("Number of patches") +
  ylab("Beta diversity") +
  labs(colour = "Heterogeneity") +
  facet_grid(ordered_matrix_dispersal~ordered_corridor_dispersal,
             labeller=ggplot2::labeller(ordered_matrix_dispersal = matrix.labs,
                                        ordered_corridor_dispersal = corridor.labs))+
  scale_x_continuous(breaks = c(1, 4, 6))+
  scale_colour_manual(values = c("homogeneous" = "#B8DE29FF", "heterogeneous" = "#2D708EFF"))+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        aspect.ratio = 1,
        text = element_text(size = 12),
        legend.position = "top")+
  annotate("segment", x = -Inf, xend = Inf, y = -Inf, yend = -Inf, size = 1)+
  annotate("segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf, size = 1)
average_beta_final_plot

#### Gamma diversity ####
#create saturated model
gamma_glm <- glm(gamma_diversity ~ (heterogeneity + patch_number + corridor_dispersal + matrix_dispersal)^2,
                 data = microcosm_data, family = "gaussian", na.action="na.fail")
summary(gamma_glm)
#visual inspection model fit
par(mfrow = c(2, 2))
plot(gamma_glm)

#model simplification by AIC
dredge_gamma_glm <- dredge(gamma_glm, evaluate = TRUE)

best_gamma_models <- subset(dredge_gamma_glm, delta < 2)

#model averaging because there are multiple best models
gamma_averages <- model.avg(best_gamma_models, fit = TRUE)
summary(gamma_averages)

#extract new model preds and add to microcosm_data
new_gamma_average_data <- microcosm_data[,c(3,5,6,7)]

microcosm_data$avg_gamma_preds <- predict(gamma_averages, new_gamma_average_data, full = TRUE, type = "response")

#add fitted values and se.fit on the link scale
new_gamma_average_data <- bind_cols(new_gamma_average_data, 
                                    setNames(as_tibble(predict(gamma_averages,
                                                               new_gamma_average_data, se.fit = TRUE)[1:2]),
                                             c('average_gamma_fit_link', 'average_gamma_se_link')))

#create the interval and backtransform
new_gamma_average_data <- new_gamma_average_data %>%
  mutate(average_gamma_fit_resp = ilink(average_gamma_fit_link),
         average_gamma_right_upr = ilink(average_gamma_fit_link + (1.96 * average_gamma_se_link)),
         average_gamma_right_lwr = ilink(average_gamma_fit_link - (1.96 * average_gamma_se_link)))

#add right_upr and right_lwr to microcosm data - makes plotting easier
microcosm_data <- microcosm_data %>%
  mutate(average_gamma_fit_resp = new_gamma_average_data$average_gamma_fit_resp,
         average_gamma_right_lwr = new_gamma_average_data$average_gamma_right_lwr,
         average_gamma_right_upr = new_gamma_average_data$average_gamma_right_upr)

average_gamma_final_plot <- ggplot(data = microcosm_data)+
  geom_line(mapping = aes(x = patch_number, y = average_gamma_fit_resp, colour = heterogeneity))+
  geom_ribbon(mapping = aes(x = patch_number,
                            ymin = average_gamma_right_lwr,
                            ymax = average_gamma_right_upr,
                            group = heterogeneity),
              alpha = 0.1) +
  geom_jitter(aes(x = patch_number, y = gamma_diversity, colour = heterogeneity),
              width = 0.2, height = 0.2)+
  xlab("Number of patches") +
  ylab("Gamma diversity") +
  labs(colour = "Heterogeneity") +
  facet_grid(ordered_matrix_dispersal~ordered_corridor_dispersal,
             labeller=labeller(ordered_matrix_dispersal = matrix.labs,
                               ordered_corridor_dispersal = corridor.labs))+
  scale_x_continuous(breaks = c(1, 4, 6))+
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5), limits = c(0, 6))+
  scale_colour_manual(values = c("homogeneous" = "#B8DE29FF", "heterogeneous" = "#2D708EFF"))+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        aspect.ratio = 1,
        text = element_text(size = 12),
        legend.position = "top")+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = 1)
average_gamma_final_plot

ggsave("Figure_3.tiff", average_gamma_final_plot, device = "tiff", units = "in", width = 5.25, height = 7.05)
