###########################
#### Predator presence ####
###########################

#load required packages
library(tidyverse)
library(MuMIn)
library(ResourceSelection) 
library(pROC)
library(viridis)
library(ggpubr)
library(dplyr)

#load the data and calculate the response variable
#presence/absence at day 21
microcosm_data <- read.csv("final_presence_microcosm_data.csv", header=T, row.names=1)

microcosm_data <- microcosm_data %>%
  mutate_if(is.character, as.factor)%>%
  mutate(matrix_dispersal = relevel(matrix_dispersal, "none"),
         corridor_dispersal = relevel(corridor_dispersal, "none"),
         heterogeneity = relevel(heterogeneity, "homogeneous"),
         #add ordered matrix and corridor dispersal for plotting facets later
         ordered_matrix_dispersal = factor(matrix_dispersal, levels = c("none", "low", "high")),
         ordered_corridor_dispersal = factor(corridor_dispersal, levels = c("none", "low", "high")),
         #generate predator response variables
         generalist_predator = replace(stentor_coeruleus, stentor_coeruleus > 0, 1),#change Stentor abundances to presence/absence
         specialist_predators = rowSums(microcosm_data[,c(13,15)]), #change specialist predator combined abundances to presence/absence
         specialist_predators = as.numeric(replace(specialist_predators, specialist_predators > 0, 1)))

#### generalist predator ####
#GLM 
#create the saturated model
gen_model1 <- glm(generalist_predator~(heterogeneity+corridor_dispersal+matrix_dispersal+patch_number)^2,
                  family="binomial", 
                  data = microcosm_data,
                  na.action="na.fail")
summary(gen_model1)

#model simplification using AIC
dredge_gen_model1 <- dredge(gen_model1, evaluate = TRUE)
#subset to select best models (within 2 AIC)
subset_gen <- subset(dredge_gen_model1, delta < 2)

# model averaging because there were multiple top models #
gen_averages <- model.avg(subset_gen, fit = TRUE)
summary(gen_averages)
#create new new data with the variables in the top models
new_data_avg_gen <- microcosm_data[,c(3,4,6)]

#add fitted values and se.fit on the link scale
new_data_avg_gen <- bind_cols(new_data_avg_gen, 
                              setNames(as_tibble(predict(gen_averages,
                                                         new_data_avg_gen, se.fit = TRUE)[1:2]),
                                       c('average_gen_fit_link', 'average_gen_se_link')))
#extract model family to get inverse link
ilink <- family(gen_model1)$linkinv

#create the interval and backtransform
new_data_avg_gen <- new_data_avg_gen %>%
  mutate(average_gen_fit_resp = ilink(average_gen_fit_link),
         average_gen_right_upr = ilink(average_gen_fit_link + (1.96 * average_gen_se_link)),
         average_gen_right_lwr = ilink(average_gen_fit_link - (1.96 * average_gen_se_link)))


#add model fit and CIs to microcosm data - makes plotting easier
microcosm_data <- microcosm_data %>%
  mutate(average_gen_fit_resp = new_data_avg_gen$average_gen_fit_resp,
         average_gen_right_lwr = new_data_avg_gen$average_gen_right_lwr,
         average_gen_right_upr <- new_data_avg_gen$average_gen_right_upr)

#facet labels
matrix.labs <- c("No dispersal", "Low mortality", "High mortality")
names(matrix.labs) <- c("none", "low" ,"high")

average_generalist_predator_plot <- 
ggplot(data = microcosm_data) +
  geom_line(aes(x = patch_number, y = average_gen_fit_resp, colour = heterogeneity)) +
  geom_ribbon(aes(x = patch_number, ymin = average_gen_right_lwr, ymax = average_gen_right_upr, group = heterogeneity), alpha = 0.1)+
  geom_jitter(aes(x = patch_number, y = generalist_predator, colour = heterogeneity),
              width = 0.2, height = 0.15)+
  ylab("Probability of generalist predator presence") +
  xlab(NULL) +
  labs(colour = "Heterogeneity") +
  facet_grid(.~ ordered_matrix_dispersal,
             labeller = labeller(ordered_matrix_dispersal = matrix.labs),
             scales = "free")+
  theme_classic()+
  scale_colour_manual(values = c ("homogeneous" = "#B8DE29FF", "heterogeneous" = "#2D708EFF"))+
  scale_x_continuous(breaks = c(1, 4, 6))+
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00))+
  theme(strip.background = element_rect(colour = "white", fill = "white",
                                        linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        text = element_text(size = 12),
        aspect.ratio = 1,
        legend.position = "top",
        axis.text.x = element_blank(),
        complete = F)+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = 1)

#### specialist predators ####
#create saturated model
spec_model1 <- glm(specialist_predator~(heterogeneity+corridor_dispersal+matrix_dispersal+patch_number)^2,
                   family="binomial", 
                   data = microcosm_data,
                   na.action="na.fail")
summary(spec_model1)

#model simplification by AIC
dredge_spec_model1 <- dredge(spec_model1, evaluate = TRUE)
#subset top models within 2 AIC of top model
subset_spec <- subset(dredge_spec_model1, delta < 2)

### model averaging as there are multiple top models ####
spec_averages <- model.avg(subset_spec, fit = TRUE)
summary(spec_averages)
#create new new data with the variables present in the top models
new_data_avg_spec <- microcosm_data[,c(3,4,6)]
#microcosm_data$avg_spec_preds <- predict(spec_averages, new_data_avg_spec, full = TRUE, type = "response")

#add fitted values and se.fit on the link scale
new_data_avg_spec <- bind_cols(new_data_avg_spec, 
                               setNames(as_tibble(predict(spec_averages,
                                                          new_data_avg_spec, se.fit = TRUE)[1:2]),
                                        c('average_spec_fit_link', 'average_spec_se_link')))
#extract model family to get inverse link
ilink <- family(spec_model1)$linkinv

#create the interval and backtransform
new_data_avg_spec <- new_data_avg_spec <-
  mutate(average_spec_fit_resp = ilink(average_spec_fit_link),
         average_spec_right_upr = ilink(average_spec_fit_link + (1.96 * average_spec_se_link)),
         average_spec_right_lwr = ilink(average_spec_fit_link - (1.96 * average_spec_se_link)))

#add right_upr and right_lwr to microcosm data - makes plotting easier
microcosm_data <- microcosm_data %>%
  mutate(average_spec_fit_resp = new_data_avg_spec$average_spec_fit_resp,
         average_spec_right_lwr <- new_data_avg_spec$average_spec_right_lwr,
         average_spec_right_upr <- new_data_avg_spec$average_spec_right_upr)

#now plot it
average_specialist_predator_plot <- ggplot(data = microcosm_data) +
  geom_line(aes(x = patch_number, y = average_spec_fit_resp, colour = heterogeneity)) +
  geom_ribbon(aes(x = patch_number, ymin = average_spec_right_lwr, ymax = average_spec_right_upr, group = heterogeneity), alpha = 0.1)+
  geom_jitter(aes(x = patch_number, y = specialist_predator, colour = heterogeneity),
              height = 0.15, width = 0.2)+
  xlab("Number of patches") + 
  ylab("Probability of specialist predator presence") +
  labs(colour = "Heterogeneity") +
  facet_grid(.~ ordered_matrix_dispersal, 
             labeller = labeller(ordered_matrix_dispersal = matrix.labs))+
  #  ylim(-0.5,1.2) +
  scale_colour_manual(values = c ("homogeneous" = "#B8DE29FF", "heterogeneous" = "#2D708EFF"))+
  theme(strip.background = element_rect(colour = "black",
                                        fill = "white",
                                        linetype = "blank"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 12),
        strip.text.x = element_blank(),
        aspect.ratio = 1,
        legend.position = "none")+
  scale_x_continuous(breaks = c(1, 4, 6))+
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size = 1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size = 1)
average_specialist_predator_plot

leg <- get_legend(average_generalist_predator_plot)

oeco_plot <- ggarrange(leg, average_generalist_predator_plot, NULL, average_specialist_predator_plot, ncol = 1,
                       labels = c("","a)","", "b)"), common.legend = T, align = "hv", label.x = 0.05, label.y = 0.7,
                       heights = c(-0.2,1, -0.55, 1), widths = c(1,1,0.1,1))


ggsave("Figure_2.tiff", oeco_plot, device = "tiff", units = "in", width = 7, height = 9.4)


