#####################
#### Extinctions ####
#####################

#load required packages
library(tidyverse)
library(MuMIn)
library(MASS)
library(data.table)

#load data and calculate the response variable
#load presence/absence data at day 21
microcosm_data <- read.csv("final_presence_microcosm_data.csv", row.names = 1, header = T)

#ensure control comes first and add extinctions variable
microcosm_data <- microcosm_data %>%
  mutate_if(is.character, as.factor)%>%
  mutate(matrix_dispersal = relevel(matrix_dispersal, "none"),
         corridor_dispersal = relevel(corridor_dispersal, "none"),
         heterogeneity = relevel(heterogeneity, "homogeneous"),
         extinctions = (8 - rowSums(microcosm_data[,8:15])))

#--------------------
#### fit GLM ####

#fit saturated model with all variables and two-way interactions
model1 <- glm(extinctions ~ (patch_number + heterogeneity + matrix_dispersal + corridor_dispersal)^2,
              family = "poisson", na.action = "na.fail", data = microcosm_data)

#dredge to select best model based on AIC
summary(get.models(dredge(model1, evaluate = TRUE), 1)[[1]])

best_model1 <- glm(extinctions ~ heterogeneity + patch_number + 1, 
                   family = "poisson", data = microcosm_data, na.action = "na.fail")
#check model fit
par(mfrow = c(2,2))
plot(best_model1)
#underdispersed 

#fit Quasi-Poisson instead because of underdispersion
#fit saturated model with all variables and two-way interactions
quasi_model1 <- glm(extinctions ~ (patch_number + heterogeneity + matrix_dispersal + corridor_dispersal)^2,
              family = "quasipoisson", na.action = "na.fail", data = microcosm_data)

#extract QAIC (AIC for quasi-poisson)
dfun <- function(object) {
  with(object, sum((weights * residuals ^ 2)[weights > 0])/df.residual)
}

x.quasipoisson <- function() {
  res <- quasipoisson()
  res$aic <- poisson()$aic
  res
}
quasi1 <- update(model1, family = "x.quasipoisson",
                 na.action = "na.fail")
gg <- dredge(quasi1, rank = "QAICc", chat = dfun(model1))
#select models within <2 AIC of the top model
subset_quasi1 <-subset(gg, delta < 2)

#### model averaging because there are 2 top models ####
quasi1_averages <- model.avg(subset_quasi1, fit = TRUE)

#extract new model preds and add to microcosm_data
#create a dataframe with variabless present in top models (patch number and heterogeneity)
new_data <- microcosm_data[,3:4]

#add fitted values and se.fit on the link scale
new_data <- bind_cols(new_data, 
                      setNames(as_tibble(predict(quasi1_averages,
                                                 new_data, se.fit = TRUE)[1:2]),
                               c('average_fit_link', 'average_se_link')))

#create the interval and backtransform
ilink <- family(quasi1)$linkinv #extract family and inverse link

#extract model fit, upper and lower confidence intervals
new_data <- new_data %>%
  mutate(average_fit_resp = ilink(average_fit_link),
         average_right_upr = ilink(average_fit_link + (1.96 * average_se_link)),
         average_right_lwr = ilink(average_fit_link - (1.96 * average_se_link)))

#add right_upr and right_lwr to microcosm data for easier plotting 
microcosm_data <- microcosm_data %>%
  mutate(average_right_lwr = new_data$average_right_lwr,
         average_right_upr = new_data$average_right_upr,
         average_fit_resp = new_data$average_fit_resp)

#plot model output 
average_final_plot <- ggplot(data = microcosm_data)+
  geom_line(aes(x = patch_number, y = average_fit_resp, colour = heterogeneity))+
  geom_ribbon(aes(x = patch_number, ymin = average_right_lwr, ymax = average_right_upr, group = heterogeneity),
              alpha = 0.1)+
  geom_jitter(aes(x = patch_number, y = extinctions, colour = heterogeneity),
              height = 0.2, width = 0.2)+
  xlab("Number of patches")+
  ylab("Number of extinctions")+
  labs(colour = "Heterogeneity")+
  scale_colour_manual(values = c ("homogeneous" = "#B8DE29FF", "heterogeneous" = "#2D708EFF"))+
  scale_x_continuous(breaks = c(1, 4, 6))+
  theme_classic()+
  theme(text = element_text(size = 12),
        legend.position = "top",
        aspect.ratio = 1,
        complete = F)

ggsave("average_final_plot.tiff", average_final_plot, device = "tiff", units = "in", width = 5, height = 5)

