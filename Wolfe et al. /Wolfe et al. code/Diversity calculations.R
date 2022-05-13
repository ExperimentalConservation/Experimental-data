######################################
#### Diversity Response Variables ####
######################################

#load required packages
library(entropart)
library(data.table)
library(dplyr)

#load in species abundance at the patch level
protist_data <- read.csv("final_count_patch_data.csv", row.names =1, header = T)
#read in microcosm info data to add diversity variables to later
microcosm_data <- read.csv("final_microcosm_data.csv", header = T)

#add weights (patch sizes)
weightsSL <- as.vector(rep(48, times = 36))
weights4Ho <- as.vector(rep(c(12,12,12,12), times = 35)) #12mL x 4 patches x 35 microcosms
weights4He <- as.vector(rep(c(16,8,16,8), times = 36)) #alternating 16 and 8mL for 36 microcosms
weights6Ho <- as.vector(rep(8, times = 204)) #8mL x 6 patches x 34 microcosms
weights6He <- as.vector(rep(c(12,4,12,4,12,4), times = 35)) #alternating 12 and 4mL for 35 microcosms

#combine them and add to protist_data
protist_data$weights <- c(weightsSL, weights4Ho, weights4He, weights6Ho, weights6He)

#remove all rows where patch is empty (diversity calculations won't work)
protist_data <- subset(protist_data, rowSums(protist_data[,c(7:14)])!=0)

#preparing protist_data_2 dataframe for SL and 6HoCL3 - separate analyses because they are only one patch
protist_data_2 <- protist_data %>%
  filter(microcosm_type == "SL"|
         microcosm_ID == "6HoLC3") %>% 
  select(-c(microcosm_type, replicate, heterogeneity, matrix_dispersal, corridor_dispersal)) #remove unnecessary columns

#preparing protist_data dataframe for all other analysis
protist_data <- protist_data %>%
  filter(!microcosm_type == "SL") %>% #remove SL and 6HoLC3 microcosms
  filter(!microcosm_ID == "6HoLC3") %>%
  select(-c(microcosm_type, replicate, heterogeneity, matrix_dispersal, corridor_dispersal))#remove uneccessary columns  

#split protist abundance and weights data, according to microcosm ID
split_data <- split(protist_data[2:10], protist_data$microcosm_ID)

#for loops to calculate alpha and beta diversity. Have removed 6HoLC3 and SL microcosms because there is only one patch
alpha_diversity <- vector("double", 139)
beta_diversity <- vector("double", 139)
gamma_diversity <- vector("double", 139)

for(i in 1:length(split_data)){
  alpha_diversity[[i]] <- as.numeric(DivPart(MC = (MetaCommunity(Abundances = t(split_data[[i]][1:8]), Weights = split_data[[i]]$weights)))$TotalAlphaDiversity)
  beta_diversity[[i]] <- as.numeric(DivPart(MC = (MetaCommunity(Abundances = t(split_data[[i]][1:8]), Weights = split_data[[i]]$weights)))$TotalBetaDiversity)
  gamma_diversity[[i]] <- as.numeric(DivPart(MC = (MetaCommunity(Abundances = t(split_data[[i]][1:8]), Weights = split_data[[i]]$weights)))$GammaDiversity)
}

#create a data frame and put the data in
diversity_data <- microcosm_data %>%
  rename(microcosm_ID = X) %>%
  subset(microcosm_type != "SL") %>%
  subset(microcosm_ID != "6HoLC3") %>%
  arrange(microcosm_ID) %>% #order alphabetically to match order of diversity measures
  mutate(alpha_diversity = alpha_diversity,
         beta_diversity = beta_diversity,
         gamma_diversity = gamma_diversity)

###########################################
#### gamma diversity for SL and 6HoLC3 ####
###########################################

SL_gamma <- DivPart(MC = (MetaCommunity(Abundances = t(protist_data_2[1:36,2:9]), Weights = protist_data_2[1:36,10])))$CommunityAlphaDiversities
sixHoLC3 <- Diversity(as.numeric(as.vector(protist_data_2[37,2:9])), q = 1)

#combine
diversity <- c(SL_gamma, sixHoLC3)

#dataframe for results
diversity_data_2 <- microcosm_data %>%
  rename(microcosm_ID = X) %>%
  filter(microcosm_type == "SL" |
         microcosm_ID == "6HoLC3") %>%
  mutate(alpha_diversity = diversity,
         beta_diversity = rep("n/a", times = 37),
         gamma_diversity = diversity)

#combine dataframes
diversity_data_fin <- rbind(diversity_data, diversity_data_2)

write.csv(diversity_data_fin, "diversity_data.csv")
