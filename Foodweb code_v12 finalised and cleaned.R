# RestREcol web beta diversity paper   
#
# R version 4.4

# Code used to combine the networks of insect-plant relationships from 60 
# woodlands and 60 grassland sites. The trophic networks of sites are combined
# from 1 to 10 sites randomly picked either from all 60 sites (control) or
# from sites aggregated in groupings based on their site age, site size
# or landscape proximity. Each time sites are combined bipartite web metrics (connectance, 
# generality, and NODF) are derived.   This process of random site selection and
# aggregation is repeated 500 times.  Mean and 95% CI are derived from these
# repeated simulations for each number of aggregated sites.

# set working directory

setwd("P:\\07583 RestREco (Restoring Resilient Ecosystems)\\Data\\Web Beta diversity Paper\\")


# the four required dependent files for the analysis

# DATA_sites_grass.csv
# For the 60 grassland sites described the site size, age and  proximity indexes for each site and the 
# treatment allocation to the relevant groups of low, medium and high values 
# for each of these.

# Data_sites_wood.csv
# For the 60 woodland sites described the site size, age and  proximity indexes for each site and the 
# treatment allocation to the relevant groups of low, medium and high values 
# for each of these.

# Wood_web_21.csv
# describes the feeding relationships between insects and trees each of the 60 woodland sites.

# PolWeb_21.csv
# describes the feeding relationships between insects and trees each of the 60 grassland sites.

#  You will also need to create in the directory the folders
# /results/Wood/Proximity1km
# /results/Wood/Age
# /results/Wood/Size
# /results/Grass/Proximity1km
# /results/Grass/Age
# /results/Grass/Size

# Required libraries
library(dplyr)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(bipartite)
library(circlize)
library(lmtest)


# a function that will calculate your confidence interval according to the t-distribution

confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

################################################################################

#  NUMBER OF ITERATIONS FOR PERMUTATIONS OF WEBs

n_iter<-500
set.seed(111)

################################################################################

# WOODLAND HERBIVORE  WEBS



# Site information
locations_df<- read.csv("DATA_sites_wood.csv")
res_sites_only<-filter(locations_df, between(Site,0,60)) # Remove six target restoration communites not considered in the analysis
res_sites_only$Age<-as.numeric(res_sites_only$Age)   # covert to numeric


################################################################################
################################################################################

# Site age 

res_sites_only<-res_sites_only %>% mutate(age_bin = cut(Age, breaks=c(0, 26, 46, 100)))  # identified woodland site age treatment levels
unique(res_sites_only$age_bin)
res_sites_only$site_order<-sapply(as.character(res_sites_only$age_bin), switch, "(0,26]" = 1, "(26,46]" = 2, "(46,100]"=3,
                                  USE.NAMES = F)    #   Gives a number to each of the site age groups

site_order<-unique(res_sites_only$site_order)    #  We now have sites split into 3 groups 
                                                # site_order is the defining  field for sub-setting the webs in the subsequent loops

# Read in herbivore webs
Pol<-read.csv("Wood_Web_21.csv")  
nrow(Pol)
nrow(res_sites_only)
Pol <- merge(Pol,res_sites_only,by.x="Site", all.x=TRUE)   #  merge Herbivore data frame and site variables including 'site_order'
nrow(Pol)
par(mfrow = c(1, 1))   #  just specify single graph format

number_of_site_groups<-max(Pol$site_order) # number of categories for which we are creating webs from site_order, e.g. sites split into 3 age class groups
max_sites_to_combine<-10 # max number of sites to combine withing each site_order group


# create data matrix to store the means & 95% CIs across the iterations for the combined number of sites combined and the grouping of sites by age (divided into three groups)
Mean<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(Mean) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(Mean) <-c("ages 0-26", "ages 27-46","ages 47+")

CI_low<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_low) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low) <-c("ages 0-26", "ages 27-46","ages 47+")


CI_hig<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_low) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low) <-c("ages 0-26", "ages 27-46","ages 47+")

# create data matrices to store the means & CIs across the control random site picks from all 60 sites
Mean_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(Mean_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(Mean_rand) <-c("rand")

CI_low_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_low_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low_rand) <-c("rand")


CI_hig_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_hig_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low_rand) <-c("rand")


#-----------------------------------------------------------------------
#  Number of plant species

# Calculating web metrics for the treatment levels

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[2]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("ages 0-26", "ages 27-46","ages 47+"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[2]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

write.csv(results_combined, "\\results\\Wood\\Age\\SRplant-500iterations.csv")
results_combined_plant_sr_age<-read.csv("results\\Wood\\Age\\SRplant-500iterations.csv")


# Plots of web metric response with accumulating sites
ggplot(data=results_combined_plant_sr_age, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                               fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Plant species richness") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)




#---------------------------------------------------------------------------------
#  Number of herbivore species

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[1]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("ages 0-26", "ages 27-46","ages 47+"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[1]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

write.csv(results_combined, "results\\Wood\\Age\\SRpol-500iterations.csv")
results_combined_polSR_age<-read.csv( "results\\Wood\\Age\\SRpol-500iterations.csv")




# Plots of web metric response with accumulating sites
ggplot(data=results_combined_polSR_age, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                            fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+ theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Herbivore species richness") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)





#-------Weighted connectance----------------------------------------------------
a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-2
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      web<-round(web,0)  #   round up web values to whole numbers
      value_metric<-networklevel(web, index="weighted connectance")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(nullmodel(web, N=1, method=1)[[1]], index="weighted connectance")  } else {null_model_vec[c]<- NA} # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("ages 0-26", "ages 27-46","ages 47+"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      web<-round(web,0)  #   round up web values to whole numbers
      value_metric<-networklevel(web, index="weighted connectance")  #  calculates actual web value
      null_model_vec <- rep(NA, times=10)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) { 
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(nullmodel(web, N=1, method=1)[[1]], index="weighted connectance")  } else {null_model_vec[c]<- NA}# calculates a null model and the relvent metric
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

write.csv(results_combined, "results\\Wood\\Age\\conect-500iterations.csv")

results_combined_connect_age<-read.csv( "results\\Wood\\Age\\conect-500iterations.csv")

# Plots of web metric response with accumulating sites
ggplot(data=results_combined_connect_age, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                               fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Weighted connectance") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)


#  Exponential Decay Model to understand rate of decline in connectance according to age classes    
Working_Data_trends<-results_combined_connect_age
low <- na.omit(Working_Data_trends %>% filter(groups =="ages 0-26"))
med <- na.omit(Working_Data_trends %>% filter(groups =="ages 27-46"))
high <- na.omit(Working_Data_trends %>% filter(groups =="ages 47+"))



#  young sites exponential decay slope
low_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = low, 
               start = list(a = min(low$mean), b = max(low$mean) - min(low$mean), c = 0.1))   
# a is the plateau (as number_combined_sites → ∞, mean value for network metric response approaches a) 
# b: Initial Drop or the difference between the initial value and the plateau. 
# c: Decay Rate. It controls how quickly the response variable decreases over time. Larger values of c result in faster decay, while smaller values result in slower decay.𝑐
summary(low_fit)
low_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                  data = low, 
                  start = list(a = min(low$mean), c = 0.1))    # null slope only model
summary(low_fit_linear)
lrtest(low_fit, low_fit_linear)  #  likelihood ratio test using library lmtest


params <- coef(low_fit)
a_low <- params["a"]
b_low <- params["b"]
c_low <- params["c"]
low_no_sites99perc_plateau <- -log(0.01) / c_low # Compute the number of sites when the response network metric is within 1% of the plateau


# medium age class
med_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = med, 
               start = list(a = min(med$mean), b = max(med$mean) - min(med$mean), c = 0.1))
med_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                    data = med, 
                    start = list(a = min(med$mean), c = 0.1))    # null slope only model
summary(med_fit_linear)
lrtest(med_fit, med_fit_linear)  #  likelihood ratio test using library lmtest
summary(med_fit)
params <- coef(med_fit)
a_med <- params["a"]
b_med <- params["b"]
c_med <- params["c"]
med_no_sites99perc_plateau <- -log(0.01) / c_med


# old age class
high_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
                data = high, 
                start = list(a = min(high$mean), b = max(high$mean) - min(high$mean), c = 0.1))
high_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                     data = high, 
                     start = list(a = min(high$mean), c = 0.1))    # null slope only model
summary(high_fit_linear)
lrtest(high_fit, high_fit_linear)  #  likelihood ratio test using library lmtest
summary(high_fit)
params <- coef(high_fit)
a_high <- params["a"]
b_high <- params["b"]
c_high <- params["c"]
high_no_sites99perc_plateau <- -log(0.01) / c_high
high_no_sites99perc_plateau


#           NODF
# 
# Another index for nestedness, calling nested nodf. 
# High values indicate nestedness. According to the analysis
# of Almeida-Neto et al. (2008, 2010), NODF is more consistent
# and “better” than usual measures of nestedness.

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="weighted NODF")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="weighted NODF") } else {null_model_vec[c]<- NA}  # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("ages 0-26", "ages 27-46","ages 47+"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="weighted NODF")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="weighted NODF") } else {null_model_vec[c]<- NA}  # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

write.csv(results_combined, "results\\Wood\\Age\\NODf-500iterations.csv")
results_combined_nodf_age<-read.csv("results\\Wood\\Age\\NODf-500iterations.csv")

# Plots of web metric response with accumulating sites
ggplot(data=results_combined_nodf_age, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                            fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Nestedness (weighted NODF)") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)

#  Linear models to assess rate of change in NODF according to age class'
# Responses within the range of the data are linear 
 
Working_Data_trends<-results_combined_nodf_age
low <- na.omit(Working_Data_trends %>% filter(groups =="ages 0-26"))
med <- na.omit(Working_Data_trends %>% filter(groups =="ages 27-46"))
high <- na.omit(Working_Data_trends %>% filter(groups =="ages 47+"))



#  young
low_fit <- lm(mean ~  Number_combined_sites,           data = low)   
summary(low_fit)
low_intercept <- coef(low_fit)[1]
low_slope <- coef(low_fit)[2]

# medium
med_fit <- lm(mean ~  Number_combined_sites,           data = med)   
summary(med_fit)
med_intercept <- coef(med_fit)[1]
med_slope <- coef(med_fit)[2]


# old
high_fit <- lm(mean ~  Number_combined_sites,           data = high)   
summary(high_fit)  #  ns
high_fit <- lm(mean ~  1,           data = high)   


high_intercept <- coef(high_fit)[1]
high_slope <- coef(high_fit)[2]


low_intercept
low_slope

med_intercept
med_slope

high_intercept
high_slope

#---------- generality----------------------------------------------------------
# 
#(Weighted) mean effective number of LL species per HL species (generality; HL
# species per LL species for vulnerability), weighted by their marginal totals (row 
#sums); see Tylianakis et al. (2007) and Bersier et al. (2002). This is identical
#to exp(“partner diversity”, i.e., simply the Jost (2006)-recommended version of
#        diversity.
a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="generality")[1]  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="generality")[1] } else {null_model_vec[c]<- NA} # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("ages 0-26", "ages 27-46","ages 47+"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="generality")[1]  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="generality")[1] } else {null_model_vec[c]<- NA} # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

write.csv(results_combined, "results\\Wood\\Age\\generality-500iterations.csv")
results_combined_gen_age<-read.csv( "results\\Wood\\Age\\generality-500iterations.csv")

# Plots of web metric response with accumulating sites
ggplot(data=results_combined_gen_age, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                           fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+ theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Generality") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)


Working_Data_trends<-results_combined_gen_age
llow <- na.omit(Working_Data_trends %>% filter(groups =="ages 0-26"))
med <- na.omit(Working_Data_trends %>% filter(groups =="ages 27-46"))
high <- na.omit(Working_Data_trends %>% filter(groups =="ages 47+"))


#  young sites exponential decay slope
low_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = low, 
               start = list(a = min(low$mean), b = max(low$mean) - min(low$mean), c = 0.1))   
# a is the plateau (as number_combined_sites → ∞, mean value for network metric response approaches a) 
# b: Initial Drop or the difference between the initial value and the plateau. 
# c: Decay Rate. It controls how quickly the response variable decreases over time. Larger values of c result in faster decay, while smaller values result in slower decay.𝑐
summary(low_fit)
low_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                  data = low, 
                  start = list(a = min(low$mean), c = 0.1))    # null slope only model
summary(low_fit_linear)
lrtest(low_fit, low_fit_linear)  #  likelihood ratio test using library lmtest


params <- coef(low_fit)
a_low <- params["a"]
b_low <- params["b"]
c_low <- params["c"]
low_no_sites99perc_plateau <- -log(0.01) / c_low # Compute the number of sites when the response network metric is within 1% of the plateau


# medium age class
med_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = med, 
               start = list(a = min(med$mean), b = max(med$mean) - min(med$mean), c = 0.1))
med_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                    data = med, 
                    start = list(a = min(med$mean), c = 0.1))    # null slope only model
summary(med_fit_linear)
lrtest(med_fit, med_fit_linear)  #  likelihood ratio test using library lmtest
summary(med_fit)
params <- coef(med_fit)
a_med <- params["a"]
b_med <- params["b"]
c_med <- params["c"]
med_no_sites99perc_plateau <- -log(0.01) / c_med


# old age class
high_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
                data = high, 
                start = list(a = min(high$mean), b = max(high$mean) - min(high$mean), c = 0.1))
high_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                     data = high, 
                     start = list(a = min(high$mean), c = 0.1))    # null slope only model
summary(high_fit_linear)
lrtest(high_fit, high_fit_linear)  #  likelihood ratio test using library lmtest
summary(high_fit)
params <- coef(high_fit)
a_high <- params["a"]
b_high <- params["b"]
c_high <- params["c"]
high_no_sites99perc_plateau <- -log(0.01) / c_high



#################################



# Linear models to assesss rate of chage in genealtiy according to age class
# Responses within the range of the data are linear 
 
Working_Data_trends<-  results_combined_gen_age2
  low <- Working_Data_trends %>% filter(groups =="ages 0-26")
med <- Working_Data_trends %>% filter(groups =="ages 27-46")
high <- Working_Data_trends %>% filter(groups =="ages 47+")


#  Young
low_fit <- lm(mean ~  Number_combined_sites,           data = low)   
summary(low_fit)
low_intercept <- coef(low_fit)[1]
low_slope <- coef(low_fit)[2]

# Medium
med_fit <- lm(mean ~  Number_combined_sites,           data = med)   
summary(med_fit)
med_intercept <- coef(med_fit)[1]
med_slope <- coef(med_fit)[2]


# Old
high_fit <- lm(mean ~  Number_combined_sites,           data = high)   
summary(high_fit)
high_intercept <- coef(high_fit)[1]
high_slope <- coef(high_fit)[2]


low_intercept
low_slope

med_intercept
med_slope

high_intercept
high_slope

################################################################################
###########  Herbivores -  Site Size    ########################################


#  Field_area_order    -  provides three equally split groups labeled 1 (small size), 2 (medium size) and 3 (laerge size)

res_sites_only$site_order<-res_sites_only$Field_area_order  #   specifies the three groups of small, medium and large sites

site_order<-unique(res_sites_only$site_order)    #  We now have sites split into 3 groups (1-3) 
#------------------------------------  site_order is the defining  field for sub-setting the webs in the subsequent loops



# Read in herbivore webs
Pol<-read.csv("Wood_Web_21.csv")
nrow(Pol)
nrow(res_sites_only)
Pol <- merge(Pol,res_sites_only,by.x="Site", all.x=TRUE)   #  merge Herbivore data frame and site variables including 'site_order'
nrow(Pol)

par(mfrow = c(1, 1))   #  just specify single graph format



number_of_site_groups<-max(Pol$site_order) # number of categories for which we are creating webs from site_order, e.g. sites split into 6 age class groups
max_sites_to_combine<-10 # max number of sites to combine withing each site_order group, e.g. cluster of age bracketed sites (e.g. 6  -  a 1/10)


# create data matricies to store the means & CIs across the interactions for each number of sites combined and the grouping of sites by environment (e.g. age divided into 6 groups)
Mean<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(Mean) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(Mean) <-c("Small", "Med.","Large")

CI_low<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_low) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low) <-c("Small", "Med.","Large")


CI_hig<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_hig) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low) <-c("Small", "Med.","Large")

# create data matricies to store the means & CIs across the Random site picks
Mean_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(Mean_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(Mean_rand) <-c("rand")

CI_low_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_low_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low_rand) <-c("rand")


CI_hig_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_hig_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low_rand) <-c("rand")


#-------------------------------------------------------------------------------
#  Number of plant species

# Calculating web metrics for the treatment levels

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[2]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Small", "Med.","Large"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[2]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)




# Plots of web metric response with accumulating sites

levels(results_combined$groups)

results_combined$groups <- factor(results_combined$groups, levels = c('Small', 'Med.', 'Large', 'Random'))

write.csv(results_combined, "results\\Wood\\Size\\SRplant-500iterations.csv")
results_combined_plantSR_size<-read.csv("results\\Wood\\Size\\SRplant-500iterations.csv")
results_combined_plantSR_size$groups <- factor(results_combined_plantSR_size$groups, levels = c('Small', 'Med.', 'Large', 'Random'))  

ggplot(data=results_combined_plantSR_size, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+  theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Plant species richness") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)



#  Number of herbivore species

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[1]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Small", "Med.","Large"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[1]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

results_combined$groups <- factor(results_combined$groups, levels = c('Small', 'Med.', 'Large', 'Random'))

write.csv(results_combined, "results\\Wood\\Size\\SRpol-500iterations.csv")
results_combined_polSR_size<-read.csv("results\\Wood\\Size\\SRpol-500iterations.csv")
results_combined_polSR_size$groups <- factor(results_combined_polSR_size$groups, levels = c('Small', 'Med.', 'Large', 'Random'))  



# Plots of web metric response with accumulating sites
ggplot(data=results_combined_polSR_size, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                             fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+  theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Herbivore species richness") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)





#-------Weighted connectance----------------------------------------------------
a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      web<-round(web,0)  #   round up web values to whole numbers
      value_metric<-networklevel(web, index="weighted connectance")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(nullmodel(web, N=1, method=1)[[1]], index="weighted connectance")  } else {null_model_vec[c]<- NA} # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Small", "Med.","Large"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      web<-round(web,0)  #   round up web values to whole numbers
      value_metric<-networklevel(web, index="weighted connectance")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to store the null model values in
      c<-1
      for (c in 1:1000) { 
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(nullmodel(web, N=1, method=1)[[1]], index="weighted connectance")  } else {null_model_vec[c]<- NA}# calculates a null model and the relvent metric
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

results_combined$groups <- factor(results_combined$groups, levels = c('Small', 'Med.', 'Large', 'Random'))

write.csv(results_combined, "results\\Wood\\Size\\conect-500iterations.csv")
results_combined_conect_size<-read.csv("results\\Wood\\Size\\conect-500iterations.csv")
results_combined_conect_size$groups <- factor(results_combined_conect_size$groups, levels = c('Small', 'Med.', 'Large', 'Random'))  


# Plots of web metric response with accumulating sites
ggplot(data=results_combined_conect_size, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                               fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Weighted connectance") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)


#  Exponential Decay Model to understand rate of decline in connectacne according to site size classes    
Working_Data_trends<-results_combined_conect_size
  Working_Data_trends <- Working_Data_trends %>% drop_na()
low <- Working_Data_trends %>% filter(groups =="Small")
med <- Working_Data_trends %>% filter(groups =="Med.")
high <- Working_Data_trends %>% filter(groups =="Large")


#  small
low_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = low, 
               start = list(a = min(low$mean), b = max(low$mean) - min(low$mean), c = 0.1))   
# a is the plateau (as number_combined_sites → ∞, mean value for network metric response approaches a) 
# b: Initial Drop or the difference between the initial value and the plateau. 
# c: Decay Rate. It controls how quickly the response variable decreases over time. Larger values of c result in faster decay, while smaller values result in slower decay.𝑐

low_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                    data = low, 
                    start = list(a = min(low$mean), c = 0.1))    # null slope only model
summary(low_fit_linear)
lrtest(low_fit, low_fit_linear)  #  liklihood ration test using library lmtest
summary(low_fit)

params <- coef(low_fit)
a_low <- params["a"]
b_low <- params["b"]
c_low <- params["c"]
low_no_sites99perc_plateau <- -log(0.01) / c_low # Compute the number of sites when the response network metric is within 1% of the plateau


# med
med_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = med, 
               start = list(a = min(med$mean), b = max(med$mean) - min(med$mean), c = 0.1))
med_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                    data = med, 
                    start = list(a = min(med$mean), c = 0.1))    # null slope only model
summary(med_fit_linear)
lrtest(med_fit, med_fit_linear)  #  likelihood ratio test using library lmtest
summary(med_fit)
params <- coef(med_fit)
a_med <- params["a"]
b_med <- params["b"]
c_med <- params["c"]
med_no_sites99perc_plateau <- -log(0.01) / c_med


# large
high_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
                data = high, 
                start = list(a = min(high$mean), b = max(high$mean) - min(high$mean), c = 0.1))
high_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                     data = high, 
                     start = list(a = min(high$mean), c = 0.1))    # null slope only model
summary(high_fit_linear)
lrtest(high_fit, high_fit_linear)  #  likelihood ratio test using library lmtest
summary(high_fit)
params <- coef(high_fit)
a_high <- params["a"]
b_high <- params["b"]
c_high <- params["c"]
high_no_sites99perc_plateau <- -log(0.01) / c_high

a_low
b_low
c_low
low_no_sites99perc_plateau

a_med
b_med
c_med
med_no_sites99perc_plateau

a_high
b_high
c_high
high_no_sites99perc_plateau




#           NODF
# 
# Another index for nestedness, calling nested nodf. 
# High values indicate nestedness. According to the analysis
# of Almeida-Neto et al. (2008, 2010), NODF is more consistent
# and “better” than usual measures of nestedness.

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="weighted NODF")  #  calculates actual web value
      null_model_vec <- rep(NA, times=100)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:100) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="weighted NODF") } else {null_model_vec[c]<- NA}  # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Small", "Med.","Large"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="weighted NODF")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="weighted NODF") } else {null_model_vec[c]<- NA}  # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
      
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

results_combined$groups <- factor(results_combined$groups, levels = c('Small', 'Med.', 'Large', 'Random'))

write.csv(results_combined, "results\\Wood\\Size\\NODf-500iterations.csv")
results_combined_nodf_size<-read.csv( "results\\Wood\\Size\\NODf-500iterations.csv")
results_combined_nodf_size$groups <- factor(results_combined_nodf_size$groups, levels = c('Small', 'Med.', 'Large', 'Random'))  



# Plots of web metric response with accumulating sites
ggplot(data=results_combined_nodf_size, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                             fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+ theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Nestedness (weighted NODF)") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)

#  Linear models to assess rate of change in NODF according to site size classes
# Responses within the range of the data are linear 
Working_Data_trends<-results_combined_nodf_size
  Working_Data_trends <- Working_Data_trends %>% drop_na()
low <- Working_Data_trends %>% filter(groups =="Small")
med <- Working_Data_trends %>% filter(groups =="Med.")
high <- Working_Data_trends %>% filter(groups =="Large")



#  small
low_fit <- lm(mean ~  Number_combined_sites,           data = low)   
summary(low_fit)
low_fit <- lm(mean ~  1,           data = low)   
summary(low_fit)

low_intercept <- coef(low_fit)[1]

# medium
med_fit <- lm(mean ~  Number_combined_sites,           data = med)   
summary(med_fit)
med_intercept <- coef(med_fit)[1]
med_slope <- coef(med_fit)[2]


# large
high_fit <- lm(mean ~  Number_combined_sites,           data = high)   
summary(high_fit)
high_fit <- lm(mean ~  1,           data = high)   
summary(high_fit)
high_intercept <- coef(high_fit)[1]

low_intercept
low_slope

med_intercept
med_slope

high_intercept
high_slope

#---------- generality----------------------------------------------------------
# 
#(Weighted) mean effective number of LL species per HL species (generality; HL
# species per LL species for vulnerability), weighted by their marginal totals (row 
#sums); see Tylianakis et al. (2007) and Bersier et al. (2002). This is identical
#to exp(“partner diversity”, i.e., simply the Jost (2006)-recommended version of
#        diversity.
a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="generality")[1]  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to store the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="generality")[1] } else {null_model_vec[c]<- NA} # Calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Small", "Med.","Large"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="generality")[1]  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="generality")[1] } else {null_model_vec[c]<- NA} # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

results_combined$groups <- factor(results_combined$groups, levels = c('Small', 'Med.', 'Large', 'Random'))

write.csv(results_combined, "results\\Wood\\Size\\generality-500iterations.csv")
results_combined_gen_size<-read.csv("results\\Wood\\Size\\generality-500iterations.csv")
results_combined_gen_size2 <- subset(results_combined_gen_size, Number_combined_sites != 1) 
results_combined_gen_size$groups <- factor(results_combined_gen_size$groups, levels = c('Small', 'Med.', 'Large', 'Random'))  



# Plots of web metric response with accumulating sites
ggplot(data=results_combined_gen_size, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                            fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+theme_bw() +
  
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Generality") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)

#  Linear models to assess rate of change in generality according to site size classes
# Responses within the range of the data are linear 
Working_Data_trends<-results_combined_gen_size2
  Working_Data_trends <- Working_Data_trends %>% drop_na()
low <- Working_Data_trends %>% filter(groups =="Small")
med <- Working_Data_trends %>% filter(groups =="Med.")
high <- Working_Data_trends %>% filter(groups =="Large")

# test of asymptote model


#  small
low_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = low, 
               start = list(a = min(low$mean), b = max(low$mean) - min(low$mean), c = 0.1))   
# a is the plateau (as number_combined_sites → ∞, mean value for network metric response approaches a) 
# b: Initial Drop or the difference between the initial value and the plateau. 
# c: Decay Rate. It controls how quickly the response variable decreases over time. Larger values of c result in faster decay, while smaller values result in slower decay.𝑐

low_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                    data = low, 
                    start = list(a = min(low$mean), c = 0.1))    # null slope only model
summary(low_fit_linear)
lrtest(low_fit, low_fit_linear)  #  liklihood ration test using library lmtest
summary(low_fit)

params <- coef(low_fit)
a_low <- params["a"]
b_low <- params["b"]
c_low <- params["c"]
low_no_sites99perc_plateau <- -log(0.01) / c_low # Compute the number of sites when the response network metric is within 1% of the plateau


# med
med_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = med, 
               start = list(a = min(med$mean), b = max(med$mean) - min(med$mean), c = 0.1))
med_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                    data = med, 
                    start = list(a = min(med$mean), c = 0.1))    # null slope only model
summary(med_fit_linear)
lrtest(med_fit, med_fit_linear)  #  likelihood ratio test using library lmtest
summary(med_fit)
params <- coef(med_fit)
a_med <- params["a"]
b_med <- params["b"]
c_med <- params["c"]
med_no_sites99perc_plateau <- -log(0.01) / c_med


# large
high_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
                data = high, 
                start = list(a = min(high$mean), b = max(high$mean) - min(high$mean), c = 0.1))
high_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                     data = high, 
                     start = list(a = min(high$mean), c = 0.1))    # null slope only model
summary(high_fit_linear)
lrtest(high_fit, high_fit_linear)  #  likelihood ratio test using library lmtest
summary(high_fit)
params <- coef(high_fit)
a_high <- params["a"]
b_high <- params["b"]
c_high <- params["c"]
high_no_sites99perc_plateau <- -log(0.01) / c_high

######Linear


#  low
low_fit <- lm(mean ~  Number_combined_sites,           data = low)   
summary(low_fit)
low_intercept <- coef(low_fit)[1]
low_slope <- coef(low_fit)[2]

# med
med_fit <- lm(mean ~  Number_combined_sites,           data = med)   
summary(med_fit)
med_intercept <- coef(med_fit)[1]
med_slope <- coef(med_fit)[2]


# high
high_fit <- lm(mean ~  Number_combined_sites,           data = high)   
summary(high_fit)
high_intercept <- coef(high_fit)[1]
high_slope <- coef(high_fit)[2]


low_intercept
low_slope

med_intercept
med_slope

high_intercept
high_slope



########################################################################################################################
###########   Proximity 1km    #########################################################################################


# Proximity_1km_order    -  provides three equally split groups labeled 1 (low proximity), 2 (medium ) and 3 (high proximity)

res_sites_only$site_order<-res_sites_only$Proximity_1km_order  #   specifies the three groups of low, medium and high proximity index

site_order<-unique(res_sites_only$site_order)    #  We now have sites split into 3 groups (1-3) 
#------------------------------------  site_order is the defining  field for sub-setting the webs in the subsequent loops


# Read in herbivore webs
Pol<-read.csv("Wood_Web_21.csv")
nrow(Pol)
nrow(res_sites_only)
Pol <- merge(Pol,res_sites_only,by.x="Site", all.x=TRUE)   #  merge Herbivore data frame and site variables including 'site_order'
nrow(Pol)

par(mfrow = c(1, 1))   #  just specify single graph format



number_of_site_groups<-max(Pol$site_order) # number of categories for which we are creating webs from site_order, e.g. sites split into 6 age class groups
max_sites_to_combine<-10 # max number of sites to combine withing each site_order group, e.g. cluster of age bracketed sites (e.g. 6  -  a 1/10)

# create data matricies to store the means & CIs across the interactions for each number of sites combined and the grouping of sites by environment (e.g. age divided into 6 groups)
Mean<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(Mean) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(Mean) <-c("Low", "Med.","High")

CI_low<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_low) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low) <-c("Low", "Med.","High")



CI_hig<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_hig) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low) <-c("Low", "Med.","High")

# create data matricies to store the means & CIs across the Random site picks
Mean_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(Mean_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(Mean_rand) <-c("rand")

CI_low_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_low_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low_rand) <-c("rand")


CI_hig_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_hig_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low_rand) <-c("rand")


#-----------------------------------------------------------------------
#  Number of plant species

# Calculating web metrics for the treatment levels

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[2]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Low", "Med.","High"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[2]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)



write.csv(results_combined, "results\\Wood\\Proximity1km\\SRplant-500iterations.csv")
results_combined_plantSR_prox<-read.csv( "results\\Wood\\Proximity1km\\SRplant-500iterations.csv")
results_combined_plantSR_prox$groups <- factor(results_combined_plantSR_prox$groups, levels = c('Low', 'Med.', 'High', 'Random'))   #  gives an order to the catagories for the legend



# Plots of web metric response with accumulating sites
ggplot(data=results_combined_plantSR_prox, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                               fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Plant species richness") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)





#  Number of herbivore species

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[1]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Low", "Med.","High"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[1]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

results_combined$groups <- factor(results_combined$groups, levels = c('Low', 'Med.', 'High', 'Random'))

write.csv(results_combined, "results\\Wood\\Proximity1km\\SRpol-500iterations.csv")
results_combined_polSR_prox<-read.csv("results\\Wood\\Proximity1km\\SRpol-500iterations.csv")
results_combined_polSR_prox$groups <- factor(results_combined_polSR_prox$groups, levels = c('Low', 'Med.', 'High', 'Random'))



# Plots of web metric response with accumulating sites
ggplot(data=results_combined_polSR_prox, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                             fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+ theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Herbivore species richness") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)





#-------Weighted connectance----------------------------------------------------
a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      web<-round(web,0)  #   round up web values to whole numbers
      value_metric<-networklevel(web, index="weighted connectance")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to store the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(nullmodel(web, N=1, method=1)[[1]], index="weighted connectance")  } else {null_model_vec[c]<- NA} # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Low", "Med.","High"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      web<-round(web,0)  #   round up web values to whole numbers
      value_metric<-networklevel(web, index="weighted connectance")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) { 
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(nullmodel(web, N=1, method=1)[[1]], index="weighted connectance")  } else {null_model_vec[c]<- NA}# calculates a null model and the relvent metric
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)



write.csv(results_combined, "results\\Wood\\Proximity1km\\conect-500iterations.csv")
results_combined_connect_prox<-read.csv( "results\\Wood\\Proximity1km\\conect-500iterations.csv")
results_combined_connect_prox$groups <- factor(results_combined_connect_prox$groups, levels = c('Low', 'Med.', 'High', 'Random'))   #  gives an order to the catagories for the legend



# Plots of web metric response with accumulating sites
ggplot(data=results_combined_connect_prox, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                                fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+ theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Weighted connectance") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)


#  Exponential Decay Model to understand rate of decline in connectance with proximity ladnsacpe classes   
Working_Data_trends<-results_combined_connect_prox
  Working_Data_trends <- Working_Data_trends %>% drop_na()
low <- Working_Data_trends %>% filter(groups =="Low")
med <- Working_Data_trends %>% filter(groups =="Med.")
high <- Working_Data_trends %>% filter(groups =="High")


#  low
low_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = low, 
               start = list(a = min(low$mean), b = max(low$mean) - min(low$mean), c = 0.1))   
# a is the plateau (as number_combined_sites → ∞, mean value for network metric response approaches a) 
# b: Initial Drop or the difference between the initial value and the plateau. 
# c: Decay Rate. It controls how quickly the response variable decreases over time. Larger values of c result in faster decay, while smaller values result in slower decay.𝑐

low_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                    data = low, 
                    start = list(a = min(low$mean), c = 0.1))    # null slope only modelsummary(low_fit_null)
lrtest(low_fit, low_fit_linear)  #  liklihood ration test using library lmtest
summary(low_fit)
params <- coef(low_fit)
a_low <- params["a"]
b_low <- params["b"]
c_low <- params["c"]
low_no_sites99perc_plateau <- -log(0.01) / c_low # Compute the number of sites when the response network metric is within 1% of the plateau


# med
med_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = med, 
               start = list(a = min(med$mean), b = max(med$mean) - min(med$mean), c = 0.1))
med_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                    data = med, 
                    start = list(a = min(med$mean), c = 0.1))    # null slope only model
summary(med_fit_linear)
lrtest(med_fit, med_fit_linear)  #  likelihood ratio test using library lmtest
summary(med_fit)
params <- coef(med_fit)
a_med <- params["a"]
b_med <- params["b"]
c_med <- params["c"]
med_no_sites99perc_plateau <- -log(0.01) / c_med


# high
high_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
                data = high, 
                start = list(a = min(high$mean), b = max(high$mean) - min(high$mean), c = 0.1))
high_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                     data = high, 
                     start = list(a = min(high$mean), c = 0.1))    # null slope only model
summary(high_fit_linear)
lrtest(high_fit, high_fit_linear)  #  likelihood ratio test using library lmtest
summary(high_fit)
params <- coef(high_fit)
a_high <- params["a"]
b_high <- params["b"]
c_high <- params["c"]
high_no_sites99perc_plateau <- -log(0.01) / c_high

a_low
b_low
c_low
low_no_sites99perc_plateau

a_med
b_med
c_med
med_no_sites99perc_plateau

a_high
b_high
c_high
high_no_sites99perc_plateau



#           NODF
# 
# Another index for nestedness, calling nested nodf. 
# High values indicate nestedness. According to the analysis
# of Almeida-Neto et al. (2008, 2010), NODF is more consistent
# and “better” than usual measures of nestedness.

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="weighted NODF")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="weighted NODF") } else {null_model_vec[c]<- NA}  # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Low", "Med.","High"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="weighted NODF")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="weighted NODF") } else {null_model_vec[c]<- NA}  # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)



write.csv(results_combined, "results\\Wood\\Proximity1km\\NODf-500iterations.csv")
results_combined_nodf_prox<-read.csv( "results\\Wood\\Proximity1km\\NODf-500iterations.csv")
results_combined_nodf_prox$groups <- factor(results_combined_nodf_prox$groups, levels = c('Low', 'Med.', 'High', 'Random'))   #  gives an order to the catagories for the legend


# Plots of web metric response with accumulating sites
ggplot(data=results_combined_nodf_prox, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                             fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+ theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Nestedness (weighted NODF)") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)

# Linear models to assess rate of change in NODF according to site landscape proximity classess
# Responses within the range of the data are linear 
Working_Data_trends<-results_combined_nodf_prox
  Working_Data_trends <- Working_Data_trends %>% drop_na()
low <- Working_Data_trends %>% filter(groups =="Low")
med <- Working_Data_trends %>% filter(groups =="Med.")
high <- Working_Data_trends %>% filter(groups =="High")

#  low
low_fit <- lm(mean ~  Number_combined_sites,           data = low)   
summary(low_fit)
low_intercept <- coef(low_fit)[1]
low_slope <- coef(low_fit)[2]

# med
med_fit <- lm(mean ~  Number_combined_sites,           data = med)   
summary(med_fit)
med_fit <- lm(mean ~  1,           data = med)   
summary(med_fit)
med_intercept <- coef(med_fit)[1]



# high
high_fit <- lm(mean ~  Number_combined_sites,           data = high)   
summary(high_fit)
high_intercept <- coef(high_fit)[1]
high_slope<- coef(high_fit)[2]

low_intercept
low_slope

med_intercept


high_intercept
high_slope

#---------- generality
# 
#(Weighted) mean effective number of LL species per HL species (generality; HL
# species per LL species for vulnerability), weighted by their marginal totals (row 
#sums); see Tylianakis et al. (2007) and Bersier et al. (2002). This is identical
#to exp(“partner diversity”, i.e., simply the Jost (2006)-recommended version of
#        diversity.
a<-1
for (a in 1:max(Pol$site_order)) {
  print(paste("site catagory:", a))
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    print(paste("number of sites:", j))
    i<-1
    for(i in 1:n_iter) {
      print(paste("Iteration:", i))
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="generality")[1]  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="generality")[1] } else {null_model_vec[c]<- NA} # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Low", "Med.","High"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a Herbivore X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="generality")[1]  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="generality")[1] } else {null_model_vec[c]<- NA} # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

results_combined$groups <- factor(results_combined$groups, levels = c('Low', 'Med.', 'High', 'Random'))

write.csv(results_combined, "results\\Wood\\Proximity1km\\generality-500iterations.csv")
results_combined_gen_prox<-read.csv( "results\\Wood\\Proximity1km\\generality-500iterations.csv")
results_combined_gen_prox$groups <- factor(results_combined_gen_prox$groups, levels = c('Low', 'Med.', 'High', 'Random'))   #  gives an order to the catagories for the legend

# Plots of web metric response with accumulating sites
ggplot(data=results_combined_gen_prox, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                            fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+ theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Generality") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)


# Linear models to assess rate of change in generality according to site landscape proximity classess
# Responses within the range of the data are linear 
Working_Data_trends<-results_combined_gen_prox
  Working_Data_trends <- Working_Data_trends %>% drop_na()
low <- Working_Data_trends %>% filter(groups =="Low")
med <- Working_Data_trends %>% filter(groups =="Med.")
high <- Working_Data_trends %>% filter(groups =="High")

#exp decay


#  low
low_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = low, 
               start = list(a = min(low$mean), b = max(low$mean) - min(low$mean), c = 0.1))   
# a is the plateau (as number_combined_sites → ∞, mean value for network metric response approaches a) 
# b: Initial Drop or the difference between the initial value and the plateau. 
# c: Decay Rate. It controls how quickly the response variable decreases over time. Larger values of c result in faster decay, while smaller values result in slower decay.𝑐

low_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                    data = low, 
                    start = list(a = min(low$mean), c = 0.1))    # null slope only modelsummary(low_fit_null)
lrtest(low_fit, low_fit_linear)  #  liklihood ration test using library lmtest
summary(low_fit)
params <- coef(low_fit)
a_low <- params["a"]
b_low <- params["b"]
c_low <- params["c"]
low_no_sites99perc_plateau <- -log(0.01) / c_low # Compute the number of sites when the response network metric is within 1% of the plateau


# med
med_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = med, 
               start = list(a = min(med$mean), b = max(med$mean) - min(med$mean), c = 0.1))
med_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                    data = med, 
                    start = list(a = min(med$mean), c = 0.1))    # null slope only model
summary(med_fit_linear)
lrtest(med_fit, med_fit_linear)  #  likelihood ratio test using library lmtest
summary(med_fit)
params <- coef(med_fit)
a_med <- params["a"]
b_med <- params["b"]
c_med <- params["c"]
med_no_sites99perc_plateau <- -log(0.01) / c_med


# high
high_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
                data = high, 
                start = list(a = min(high$mean), b = max(high$mean) - min(high$mean), c = 0.1))
high_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                     data = high, 
                     start = list(a = min(high$mean), c = 0.1))    # null slope only model
summary(high_fit_linear)
lrtest(high_fit, high_fit_linear)  #  likelihood ratio test using library lmtest
summary(high_fit)
params <- coef(high_fit)
a_high <- params["a"]
b_high <- params["b"]
c_high <- params["c"]
high_no_sites99perc_plateau <- -log(0.01) / c_high

#linear

#  low
low_fit <- lm(mean ~  Number_combined_sites,           data = low)   
summary(low_fit)
low_intercept <- coef(low_fit)[1]
low_slope <- coef(low_fit)[2]

# med
med_fit <- lm(mean ~  Number_combined_sites,           data = med)   
summary(med_fit)

med_fit <- lm(mean ~  1,           data = med)   
summary(med_fit)
med_intercept <- coef(med_fit)[1]



# high
high_fit <- lm(mean ~  Number_combined_sites,           data = high)   
summary(high_fit)
high_intercept <- coef(high_fit)[1]
high_slope <- coef(high_fit)[2]


low_intercept
low_slope

med_intercept
med_slope

high_intercept
high_slope


################################################################################################################################################################################################################
################################################################################################################################################################################################################

# POLLINATOR WEBS



setwd("P:\\07583 RestREco (Restoring Resilient Ecosystems)\\Data\\Web Beta diversity Paper\\")
# Site information
locations_df<- read.csv("DATA_sites_grass.csv")
locations_df$Field_area<-locations_df$Field_area_checked/10000  #  convert field area to ha
res_sites_only<-filter(locations_df, between(Site,0,60)) # Remove six target restoration communites not considered in the analysis
res_sites_only$Age<-as.numeric(res_sites_only$Age)   # covert to numeric


################################################################################
################################################################################

# Site age 

res_sites_only<-res_sites_only %>% mutate(age_bin = cut(Age, breaks=c(0, 8, 18, 100)))  #  establishes three age classes of young, mid and old based on approximately even splits across the 60 sites
# most bins have 8 or more sites
unique(res_sites_only$age_bin)

res_sites_only$site_order<-sapply(as.character(res_sites_only$age_bin), switch, "(0,8]" = 1, "(8,18]" = 2, "(18,100]"=3,
                                  USE.NAMES = F)

site_order<-unique(res_sites_only$site_order)    #  We now have sites split into 3 groups (1-3) 
#------------------------------------  site_order is the defining  field for sub-setting the webs in the subsequent loops

# Read in pollinator webs
Pol<-read.csv("PolWeb_21.csv")
nrow(Pol)
nrow(res_sites_only)
Pol <- merge(Pol,res_sites_only,by.x="Site", all.x=TRUE)   #  merge pollinator data frame and site variables including 'site_order'
nrow(Pol)

par(mfrow = c(1, 1))   #  just specify single graph format



number_of_site_groups<-max(Pol$site_order) # number of categories for which we are creating webs from site_order, e.g. sites split into 6 age class groups
max_sites_to_combine<-10 # max number of sites to combine withing each site_order group, e.g. cluster of age bracketed sites (e.g. 6  -  a 1/10)


# create data matricies to store the means & CIs across the interactions for each number of sites combined and the grouping of sites by environment (e.g. age divided into 6 groups)
Mean<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(Mean) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(Mean) <-c(" ages 0-8", " ages 9-18","ages 19+")

CI_low<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_low) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low) <-c(" ages 0-8", " ages 9-18","ages 19+")


CI_hig<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_low) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low) <-c(" ages 0-8", " ages 9-18","ages 19+")

# create data matricies to store the means & CIs across the Random site picks
Mean_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(Mean_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(Mean_rand) <-c("rand")

CI_low_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_low_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low_rand) <-c("rand")


CI_hig_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_hig_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low_rand) <-c("rand")


#-----------------------------------------------------------------------
#  Number of plant species

# Calculating web metrics for the treatment levels

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[2]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c(" ages 0-8", " ages 9-18","ages 19+"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[2]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

write.csv(results_combined, "results\\Grass\\Age\\SRplant-500iterations.csv")
results_combined_plant_sr_age<-read.csv("results\\Grass\\Age\\SRplant-500iterations.csv")




# Plots of web metric response with accumulating sites
ggplot(data=results_combined_plant_sr_age, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                               fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Plant species richness") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)





#  Number of pollinator species

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[1]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c(" ages 0-8", " ages 9-18","ages 19+"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[1]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

write.csv(results_combined, "results\\Grass\\Age\\SRpol-500iterations.csv")
results_combined_polSR_age<-read.csv( "results\\Grass\\Age\\SRpol-500iterations.csv")




# Plots of web metric response with accumulating sites
ggplot(data=results_combined_polSR_age, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                            fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+ theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Pollinator species richness") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)





#-------Weighted connectance----------------------------------
a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-2
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      web<-round(web,0)  #   round up web values to whole numbers
      value_metric<-networklevel(web, index="weighted connectance")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(nullmodel(web, N=1, method=1)[[1]], index="weighted connectance")  } else {null_model_vec[c]<- NA} # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c(" ages 0-8", " ages 9-18","ages 19+"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      web<-round(web,0)  #   round up web values to whole numbers
      value_metric<-networklevel(web, index="weighted connectance")  #  calculates actual web value
      null_model_vec <- rep(NA, times=10)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) { 
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(nullmodel(web, N=1, method=1)[[1]], index="weighted connectance")  } else {null_model_vec[c]<- NA}# calculates a null model and the relvent metric
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

write.csv(results_combined, "results\\Grass\\Age\\conect-500iterations.csv")
results_combined_connect_age<-read.csv( "results\\Grass\\Age\\conect-500iterations.csv")

#results_combined_connect_age2 <- subset(results_combined_connect_age, Number_combined_sites != 1) 




# Plots of web metric response with accumulating sites
ggplot(data=results_combined_connect_age, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                               fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Weighted connectance") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)



#  Exponential Decay Model to understand rate of decline according connectance with site age class
Working_Data_trends<-results_combined_connect_age
Working_Data_trends <- Working_Data_trends %>% drop_na()
low <- Working_Data_trends %>% filter(groups ==" ages 0-8")
med <- Working_Data_trends %>% filter(groups ==" ages 9-18")
high <- Working_Data_trends %>% filter(groups =="ages 19+")

Working_Data_trends$groups

#  young
low_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = low, 
               start = list(a = min(low$mean), b = max(low$mean) - min(low$mean), c = 0.1))   
# a is the plateau (as number_combined_sites → ∞, mean value for network metric response approaches a) 
# b: Initial Drop or the difference between the initial value and the plateau. 
# c: Decay Rate. It controls how quickly the response variable decreases over time. Larger values of c result in faster decay, while smaller values result in slower decay.𝑐
summary(low_fit)
low_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                    data = low, 
                    start = list(a = min(low$mean), c = 0.1))    # null slope only model
summary(low_fit_linear)
lrtest(low_fit, low_fit_linear)  #  liklihood ration test using library lmtest


params <- coef(low_fit)
a_low <- params["a"]
b_low <- params["b"]
c_low <- params["c"]
low_no_sites99perc_plateau <- -log(0.01) / c_low # Compute the number of sites when the response network metric is within 1% of the plateau


# medium
med_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = med, 
               start = list(a = min(med$mean), b = max(med$mean) - min(med$mean), c = 0.1))
med_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                    data = med, 
                    start = list(a = min(med$mean), c = 0.1))    # null slope only model
summary(med_fit_linear)
lrtest(med_fit, med_fit_linear)  #  likelihood ratio test using library lmtest
summary(med_fit)
params <- coef(med_fit)
a_med <- params["a"]
b_med <- params["b"]
c_med <- params["c"]
med_no_sites99perc_plateau <- -log(0.01) / c_med


# old
high_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
                data = high, 
                start = list(a = min(high$mean), b = max(high$mean) - min(high$mean), c = 0.1))
high_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                     data = high, 
                     start = list(a = min(high$mean), c = 0.1))    # null slope only model
summary(high_fit_linear)
lrtest(high_fit, high_fit_linear)  #  likelihood ratio test using library lmtest
summary(high_fit)
params <- coef(high_fit)
a_high <- params["a"]
b_high <- params["b"]
c_high <- params["c"]
high_no_sites99perc_plateau <- -log(0.01) / c_high

a_low
b_low
c_low
low_no_sites99perc_plateau

a_med
b_med
c_med
med_no_sites99perc_plateau

a_high
b_high
c_high
high_no_sites99perc_plateau



#           NODF
# 
# Another index for nestedness, calling nested nodf. 
# High values indicate nestedness. According to the analysis
# of Almeida-Neto et al. (2008, 2010), NODF is more consistent
# and “better” than usual measures of nestedness.

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="weighted NODF")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="weighted NODF") } else {null_model_vec[c]<- NA}  # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c(" ages 0-8", " ages 9-18","ages 19+"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="weighted NODF")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="weighted NODF") } else {null_model_vec[c]<- NA}  # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

write.csv(results_combined, "results\\Grass\\Age\\NODf-500iterations.csv")
results_combined_nodf_age<-read.csv("results\\Grass\\Age\\NODf-500iterations.csv")
#results_combined_nodf_age2 <- subset(results_combined_nodf_age, Number_combined_sites != 1) 



# Plots of web metric response with accumulating sites
ggplot(data=results_combined_nodf_age, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                            fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Nestedness (weighted NODF)") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)


# Linear models to assess rate of change in NODF according to site age classes
# Responses within the range of the data are linear 
Working_Data_trends<-results_combined_nodf_age
  Working_Data_trends <- Working_Data_trends %>% drop_na()
low <- Working_Data_trends %>% filter(groups ==" ages 0-8")
med <- Working_Data_trends %>% filter(groups ==" ages 9-18")
high <- Working_Data_trends %>% filter(groups =="ages 19+")


#  Young
low_fit <- lm(mean ~  Number_combined_sites,           data = low)   
summary(low_fit)
low_intercept <- coef(low_fit)[1]
low_slope <- coef(low_fit)[2]

# medium
med_fit <- lm(mean ~  Number_combined_sites,           data = med)   
summary(med_fit)
med_intercept <- coef(med_fit)[1]
med_slope <- coef(med_fit)[2]


# old
high_fit <- lm(mean ~  Number_combined_sites,           data = high)   
summary(high_fit)
high_intercept <- coef(high_fit)[1]
high_slope <- coef(high_fit)[2]


low_intercept
low_slope

med_intercept
med_slope

high_intercept
high_slope


#---------- generality
# 
#(Weighted) mean effective number of LL species per HL species (generality; HL
# species per LL species for vulnerability), weighted by their marginal totals (row 
#sums); see Tylianakis et al. (2007) and Bersier et al. (2002). This is identical
#to exp(“partner diversity”, i.e., simply the Jost (2006)-recommended version of
#        diversity.
a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="generality")[1]  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="generality")[1] } else {null_model_vec[c]<- NA} # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c(" ages 0-8", " ages 9-18","ages 19+"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="generality")[1]  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="generality")[1] } else {null_model_vec[c]<- NA} # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

write.csv(results_combined, "results\\Grass\\Age\\generality-500iterations.csv")
results_combined_gen_age<-read.csv( "results\\Grass\\Age\\generality-500iterations.csv")
#results_combined_gen_age2 <- subset(results_combined_gen_age, Number_combined_sites != 1) 



# Plots of web metric response with accumulating sites
ggplot(data=results_combined_gen_age, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                           fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+ theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Generality") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)


# Linear models to assess rate of change in generality according to site age classess
# Responses within the range of the data are linear 
Working_Data_trends<-results_combined_gen_age
  Working_Data_trends <- Working_Data_trends %>% drop_na()
low <- Working_Data_trends %>% filter(groups ==" ages 0-8")
med <- Working_Data_trends %>% filter(groups ==" ages 9-18")
high <- Working_Data_trends %>% filter(groups =="ages 19+")


#  Young
low_fit <- lm(mean ~  Number_combined_sites,           data = low)   
summary(low_fit)
low_intercept <- coef(low_fit)[1]
low_slope <- coef(low_fit)[2]

# medium
med_fit <- lm(mean ~  Number_combined_sites,           data = med)   
summary(med_fit)
med_intercept <- coef(med_fit)[1]
med_slope <- coef(med_fit)[2]


# old
high_fit <- lm(mean ~  Number_combined_sites,           data = high)   
summary(high_fit)
high_intercept <- coef(high_fit)[1]
high_slope <- coef(high_fit)[2]


low_intercept
low_slope

med_intercept
med_slope

high_intercept
high_slope



################################################################################
###########  Pollinators -  Site Size    #######################################


#  Field_area_order    -  provides three equally split groups labeled 1 (small size), 2 (medium size) and 3 (laerge size)

res_sites_only$site_order<-res_sites_only$Field_area_order  #   specifies the three groups of small, medium and large sites

site_order<-unique(res_sites_only$site_order)    #  We now have sites split into 3 groups (1-3) 
#------------------------------------  site_order is the defining  field for sub-setting the webs in the subsequent loops



# Read in pollinator webs
Pol<-read.csv("PolWeb_21.csv")
nrow(Pol)
nrow(res_sites_only)
Pol <- merge(Pol,res_sites_only,by.x="Site", all.x=TRUE)   #  merge pollinator data frame and site variables including 'site_order'
nrow(Pol)

par(mfrow = c(1, 1))   #  just specify single graph format



number_of_site_groups<-max(Pol$site_order) # number of categories for which we are creating webs from site_order, e.g. sites split into 6 age class groups
max_sites_to_combine<-10 # max number of sites to combine withing each site_order group, e.g. cluster of age bracketed sites (e.g. 6  -  a 1/10)


# create data matricies to store the means & CIs across the interactions for each number of sites combined and the grouping of sites by environment (e.g. age divided into 6 groups)
Mean<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(Mean) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(Mean) <-c("Small", "Med.","Large")

CI_low<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_low) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low) <-c("Small", "Med.","Large")


CI_hig<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_hig) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low) <-c("Small", "Med.","Large")

# create data matricies to store the means & CIs across the Random site picks
Mean_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(Mean_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(Mean_rand) <-c("rand")

CI_low_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_low_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low_rand) <-c("rand")


CI_hig_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_hig_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low_rand) <-c("rand")


#-----------------------------------------------------------------------
#  Number of plant species

# Calculating web metrics for the treatment levels

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[2]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Small", "Med.","Large"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[2]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)




# Plots of web metric response with accumulating sites

levels(results_combined$groups)

results_combined$groups <- factor(results_combined$groups, levels = c('Small', 'Med.', 'Large', 'Random'))

write.csv(results_combined, "results\\Grass\\Size\\SRplant-500iterations.csv")
results_combined_plantSR_size<-read.csv("results\\Grass\\Size\\SRplant-500iterations.csv")
results_combined_plantSR_size$groups <- factor(results_combined_plantSR_size$groups, levels = c('Small', 'Med.', 'Large', 'Random'))  

ggplot(data=results_combined_plantSR_size, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+  theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Plant species richness") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)



#  Number of pollinator species

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[1]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Small", "Med.","Large"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[1]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

results_combined$groups <- factor(results_combined$groups, levels = c('Small', 'Med.', 'Large', 'Random'))

write.csv(results_combined, "results\\Grass\\Size\\SRpol-500iterations.csv")
results_combined_polSR_size<-read.csv("results\\Grass\\Size\\SRpol-500iterations.csv")
results_combined_polSR_size$groups <- factor(results_combined_polSR_size$groups, levels = c('Small', 'Med.', 'Large', 'Random'))  



# Plots of web metric response with accumulating sites
ggplot(data=results_combined_polSR_size, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                             fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+  theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Pollinator species richness") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)





#-------Weighted connectance----------------------------------
a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      web<-round(web,0)  #   round up web values to whole numbers
      value_metric<-networklevel(web, index="weighted connectance")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(nullmodel(web, N=1, method=1)[[1]], index="weighted connectance")  } else {null_model_vec[c]<- NA} # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Small", "Med.","Large"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      web<-round(web,0)  #   round up web values to whole numbers
      value_metric<-networklevel(web, index="weighted connectance")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to store the null model values in
      c<-1
      for (c in 1:1000) { 
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(nullmodel(web, N=1, method=1)[[1]], index="weighted connectance")  } else {null_model_vec[c]<- NA}# calculates a null model and the relvent metric
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

results_combined$groups <- factor(results_combined$groups, levels = c('Small', 'Med.', 'Large', 'Random'))

write.csv(results_combined, "results\\Grass\\Size\\conect-500iterations.csv")
results_combined_conect_size<-read.csv("results\\Grass\\Size\\conect-500iterations.csv")
#results_combined_conect_size2 <- subset(results_combined_conect_size, Number_combined_sites != 1) 
results_combined_conect_size$groups <- factor(results_combined_conect_size$groups, levels = c('Small', 'Med.', 'Large', 'Random'))  




# Plots of web metric response with accumulating sites
ggplot(data=results_combined_conect_size, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                               fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Weighted connectance") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)


#  Exponential Decay Model to understand rate of decline  in connectance with site size classes 
Working_Data_trends<-results_combined_conect_size
  Working_Data_trends <- Working_Data_trends %>% drop_na()
low <- Working_Data_trends %>% filter(groups =="Small")
med <- Working_Data_trends %>% filter(groups =="Med.")
high <- Working_Data_trends %>% filter(groups =="Large")


#  small
low_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = low, 
               start = list(a = min(low$mean), b = max(low$mean) - min(low$mean), c = 0.1))   
# a is the plateau (as number_combined_sites → ∞, mean value for network metric response approaches a) 
# b: Initial Drop or the difference between the initial value and the plateau. 
# c: Decay Rate. It controls how quickly the response variable decreases over time. Larger values of c result in faster decay, while smaller values result in slower decay.𝑐

low_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                    data = low, 
                    start = list(a = min(low$mean), c = 0.1))    # null slope only modelsummary(low_fit_null)
lrtest(low_fit, low_fit_linear)  #  liklihood ration test using library lmtest
summary(low_fit)
params <- coef(low_fit)
a_low <- params["a"]
b_low <- params["b"]
c_low <- params["c"]
low_no_sites99perc_plateau <- -log(0.01) / c_low # Compute the number of sites when the response network metric is within 1% of the plateau


# medium
med_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = med, 
               start = list(a = min(med$mean), b = max(med$mean) - min(med$mean), c = 0.1))
med_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                    data = med, 
                    start = list(a = min(med$mean), c = 0.1))    # null slope only model
summary(med_fit_linear)
lrtest(med_fit, med_fit_linear)  #  likelihood ratio test using library lmtest
summary(med_fit)
params <- coef(med_fit)
a_med <- params["a"]
b_med <- params["b"]
c_med <- params["c"]
med_no_sites99perc_plateau <- -log(0.01) / c_med


# large
high_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
                data = high, 
                start = list(a = min(high$mean), b = max(high$mean) - min(high$mean), c = 0.1))
high_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                     data = high, 
                     start = list(a = min(high$mean), c = 0.1))    # null slope only model
summary(high_fit_linear)
lrtest(high_fit, high_fit_linear)  #  likelihood ratio test using library lmtest
summary(high_fit)
params <- coef(high_fit)
a_high <- params["a"]
b_high <- params["b"]
c_high <- params["c"]
high_no_sites99perc_plateau <- -log(0.01) / c_high

a_low
b_low
c_low
low_no_sites99perc_plateau

a_med
b_med
c_med
med_no_sites99perc_plateau

a_high
b_high
c_high
high_no_sites99perc_plateau





#           NODF
# 
# Another index for nestedness, calling nested nodf. 
# High values indicate nestedness. According to the analysis
# of Almeida-Neto et al. (2008, 2010), NODF is more consistent
# and “better” than usual measures of nestedness.

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="weighted NODF")  #  calculates actual web value
      null_model_vec <- rep(NA, times=100)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:100) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="weighted NODF") } else {null_model_vec[c]<- NA}  # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Small", "Med.","Large"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="weighted NODF")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="weighted NODF") } else {null_model_vec[c]<- NA}  # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
      
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

results_combined$groups <- factor(results_combined$groups, levels = c('Small', 'Med.', 'Large', 'Random'))

write.csv(results_combined, "results\\Grass\\Size\\NODf-500iterations.csv")
results_combined_nodf_size<-read.csv( "results\\Grass\\Size\\NODf-500iterations.csv")
#results_combined_nodf_size2 <- subset(results_combined_nodf_size, Number_combined_sites != 1) 
results_combined_nodf_size$groups <- factor(results_combined_nodf_size$groups, levels = c('Small', 'Med.', 'Large', 'Random'))  



# Plots of web metric response with accumulating sites
ggplot(data=results_combined_nodf_size, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                             fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+ theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Nestedness (weighted NODF)") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)


# Linear models to assess rate of change in NODF according to site size classes
# Responses within the range of the data are linear 
Working_Data_trends<-results_combined_nodf_size
  Working_Data_trends <- Working_Data_trends %>% drop_na()
low <- Working_Data_trends %>% filter(groups =="Small")
med <- Working_Data_trends %>% filter(groups =="Med.")
high <- Working_Data_trends %>% filter(groups =="Large")

#  small
low_fit <- lm(mean ~  Number_combined_sites,           data = low)   
summary(low_fit)
low_intercept <- coef(low_fit)[1]
low_slope <- coef(low_fit)[2]


# medium
med_fit <- lm(mean ~  Number_combined_sites,           data = med)   
summary(med_fit)
med_intercept <- coef(med_fit)[1]
med_slope <- coef(med_fit)[2]


# large
high_fit <- lm(mean ~  Number_combined_sites,           data = high)   
summary(high_fit)
high_fit <- lm(mean ~  1,           data = high)   
summary(high_fit)
high_intercept <- coef(high_fit)[1]

low_intercept


med_intercept
med_slope

high_intercept



#---------- generality
# 
#(Weighted) mean effective number of LL species per HL species (generality; HL
# species per LL species for vulnerability), weighted by their marginal totals (row 
#sums); see Tylianakis et al. (2007) and Bersier et al. (2002). This is identical
#to exp(“partner diversity”, i.e., simply the Jost (2006)-recommended version of
#        diversity.
a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="generality")[1]  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to store the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="generality")[1] } else {null_model_vec[c]<- NA} # Calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Small", "Med.","Large"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="generality")[1]  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="generality")[1] } else {null_model_vec[c]<- NA} # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

results_combined$groups <- factor(results_combined$groups, levels = c('Small', 'Med.', 'Large', 'Random'))

write.csv(results_combined, "results\\Grass\\Size\\generality-500iterations.csv")
results_combined_gen_size<-read.csv("results\\Grass\\Size\\generality-500iterations.csv")
#results_combined_gen_size2 <- subset(results_combined_gen_size, Number_combined_sites != 1) 
results_combined_gen_size$groups <- factor(results_combined_gen_size$groups, levels = c('Small', 'Med.', 'Large', 'Random'))  



# Plots of web metric response with accumulating sites
ggplot(data=results_combined_gen_size, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                            fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+theme_bw() +
  
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Generality") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)


# Linear models to assess rate of change in generality according to site size classes
# Responses within the range of the data are linear 
Working_Data_trends<-results_combined_gen_size
  Working_Data_trends <- Working_Data_trends %>% drop_na()
low <- Working_Data_trends %>% filter(groups =="Small")
med <- Working_Data_trends %>% filter(groups =="Med.")
high <- Working_Data_trends %>% filter(groups =="Large")

#  small
low_fit <- lm(mean ~  Number_combined_sites,           data = low)   
summary(low_fit)
low_intercept <- coef(low_fit)[1]
low_slope <- coef(low_fit)[2]

# medium
med_fit <- lm(mean ~  Number_combined_sites,           data = med)   
summary(med_fit)
med_intercept <- coef(med_fit)[1]
med_slope <- coef(med_fit)[2]


# large
high_fit <- lm(mean ~  Number_combined_sites,           data = high)   
summary(high_fit)
high_intercept <- coef(high_fit)[1]
high_slope <- coef(high_fit)[2]


low_intercept
low_slope

med_intercept
med_slope

high_intercept
high_slope



################################################################################
###########   Proximity 1km    #################################################


# Proximity_1km_order    -  provides three equally split groups labeled 1 (low proximity), 2 (medium ) and 3 (high proximity)

res_sites_only$site_order<-res_sites_only$Proximity_1km_order  #   specifies the three groups of low, medium and high proximity index

site_order<-unique(res_sites_only$site_order)    #  We now have sites split into 3 groups (1-3) 
#------------------------------------  site_order is the defining  field for sub-setting the webs in the subsequent loops


# Read in pollinator webs
Pol<-read.csv("PolWeb_21.csv")
nrow(Pol)
nrow(res_sites_only)
Pol <- merge(Pol,res_sites_only,by.x="Site", all.x=TRUE)   #  merge pollinator data frame and site variables including 'site_order'
nrow(Pol)

par(mfrow = c(1, 1))   #  just specify single graph format



number_of_site_groups<-max(Pol$site_order) # number of categories for which we are creating webs from site_order, e.g. sites split into 6 age class groups
max_sites_to_combine<-10 # max number of sites to combine withing each site_order group, e.g. cluster of age bracketed sites (e.g. 6  -  a 1/10)

# create data matricies to store the means & CIs across the interactions for each number of sites combined and the grouping of sites by environment (e.g. age divided into 6 groups)
Mean<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(Mean) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(Mean) <-c("Low", "Med.","High")

CI_low<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_low) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low) <-c("Low", "Med.","High")



CI_hig<-matrix(ncol = max_sites_to_combine, nrow = number_of_site_groups)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_hig) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low) <-c("Low", "Med.","High")

# create data matricies to store the means & CIs across the Random site picks
Mean_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(Mean_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(Mean_rand) <-c("rand")

CI_low_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_low_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low_rand) <-c("rand")


CI_hig_rand<-matrix(ncol = max_sites_to_combine, nrow = 1)  #  empty data frame to put in the web metrics for each combination of sites and iterations
colnames(CI_hig_rand) <- c("1_site", "2_sites", "3_sites", "4_sites", "5_sites", "6_sites", "7_sites", "8_sites", "9_sites", "10_sites")
rownames(CI_low_rand) <-c("rand")


#-----------------------------------------------------------------------
#  Number of plant species

# Calculating web metrics for the treatment levels

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[2]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Low", "Med.","High"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[2]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)



write.csv(results_combined, "results\\Grass\\Proximity1km\\SRplant-500iterations.csv")
results_combined_plantSR_prox<-read.csv( "results\\Grass\\Proximity1km\\SRplant-500iterations.csv")
results_combined_plantSR_prox$groups <- factor(results_combined_plantSR_prox$groups, levels = c('Low', 'Med.', 'High', 'Random'))


# Plots of web metric response with accumulating sites
ggplot(data=results_combined_plantSR_prox, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                               fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Plant species richness") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)





#  Number of pollinator species

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[1]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Low", "Med.","High"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      repsXsites[i,n_sites]<- networklevel(web, index="number of species")[1]  #  small webs don't work well
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

results_combined$groups <- factor(results_combined$groups, levels = c('Low', 'Med.', 'High', 'Random'))

write.csv(results_combined, "results\\Grass\\Proximity1km\\SRpol-500iterations.csv")
results_combined_polSR_prox<-read.csv("results\\Grass\\Proximity1km\\SRpol-500iterations.csv")
results_combined_polSR_prox$groups <- factor(results_combined_polSR_prox$groups, levels = c('Low', 'Med.', 'High', 'Random'))




# Plots of web metric response with accumulating sites
ggplot(data=results_combined_polSR_prox, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                             fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+ theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Pollinator species richness") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)





#-------Weighted connectance----------------------------------
a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      web<-round(web,0)  #   round up web values to whole numbers
      value_metric<-networklevel(web, index="weighted connectance")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to store the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(nullmodel(web, N=1, method=1)[[1]], index="weighted connectance")  } else {null_model_vec[c]<- NA} # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Low", "Med.","High"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      web<-round(web,0)  #   round up web values to whole numbers
      value_metric<-networklevel(web, index="weighted connectance")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) { 
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(nullmodel(web, N=1, method=1)[[1]], index="weighted connectance")  } else {null_model_vec[c]<- NA}# calculates a null model and the relvent metric
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)



write.csv(results_combined, "results\\Grass\\Proximity1km\\conect-500iterations.csv")
results_combined_connect_prox<-read.csv( "results\\Grass\\Proximity1km\\conect-500iterations.csv")
#results_combined_connect_prox2 <- subset(results_combined_connect_prox, Number_combined_sites != 1) 
results_combined_connect_prox$groups <- factor(results_combined_connect_prox$groups, levels = c('Low', 'Med.', 'High', 'Random'))


# Plots of web metric response with accumulating sites
ggplot(data=results_combined_connect_prox, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                                fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+ theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Weighted connectance") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)


#  Exponential Decay Model to understand rate of decline in connectance with site landscape proximity classs   
Working_Data_trends<-results_combined_connect_prox
  Working_Data_trends <- Working_Data_trends %>% drop_na()
low <- Working_Data_trends %>% filter(groups =="Low")
med <- Working_Data_trends %>% filter(groups =="Med.")
high <- Working_Data_trends %>% filter(groups =="High")


#  low
low_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = low, 
               start = list(a = min(low$mean), b = max(low$mean) - min(low$mean), c = 0.1))   
# a is the plateau (as number_combined_sites → ∞, mean value for network metric response approaches a) 
# b: Initial Drop or the difference between the initial value and the plateau. 
# c: Decay Rate. It controls how quickly the response variable decreases over time. Larger values of c result in faster decay, while smaller values result in slower decay.𝑐

low_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                    data = low, 
                    start = list(a = min(low$mean), c = 0.1))    # null slope only modelsummary(low_fit_null)
lrtest(low_fit, low_fit_linear)  #  liklihood ration test using library lmtest
summary(low_fit)
params <- coef(low_fit)
a_low <- params["a"]
b_low <- params["b"]
c_low <- params["c"]
low_no_sites99perc_plateau <- -log(0.01) / c_low # Compute the number of sites when the response network metric is within 1% of the plateau


# med
med_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
               data = med, 
               start = list(a = min(med$mean), b = max(med$mean) - min(med$mean), c = 0.1))
med_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                    data = med, 
                    start = list(a = min(med$mean), c = 0.1))    # null slope only model
summary(med_fit_linear)
lrtest(med_fit, med_fit_linear)  #  likelihood ratio test using library lmtest
summary(med_fit)
params <- coef(med_fit)
a_med <- params["a"]
b_med <- params["b"]
c_med <- params["c"]
med_no_sites99perc_plateau <- -log(0.01) / c_med


# high
high_fit <- nls(mean ~ a + b * exp(-c * Number_combined_sites), 
                data = high, 
                start = list(a = min(high$mean), b = max(high$mean) - min(high$mean), c = 0.1))
high_fit_linear<-nls(mean ~ a + (c * Number_combined_sites), 
                     data = high, 
                     start = list(a = min(high$mean), c = 0.1))    # null slope only model
summary(high_fit_linear)
lrtest(high_fit, high_fit_linear)  #  likelihood ratio test using library lmtest
summary(high_fit)
params <- coef(high_fit)
a_high <- params["a"]
b_high <- params["b"]
c_high <- params["c"]
high_no_sites99perc_plateau <- -log(0.01) / c_high

a_low
b_low
c_low
low_no_sites99perc_plateau

a_med
b_med
c_med
med_no_sites99perc_plateau

a_high
b_high
c_high
high_no_sites99perc_plateau




#           NODF
# 
# Another index for nestedness, calling nested nodf. 
# High values indicate nestedness. According to the analysis
# of Almeida-Neto et al. (2008, 2010), NODF is more consistent
# and “better” than usual measures of nestedness.

a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="weighted NODF")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="weighted NODF") } else {null_model_vec[c]<- NA}  # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Low", "Med.","High"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="weighted NODF")  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="weighted NODF") } else {null_model_vec[c]<- NA}  # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)



write.csv(results_combined, "results\\Grass\\Proximity1km\\NODf-500iterations.csv")
results_combined_nodf_prox<-read.csv( "results\\Grass\\Proximity1km\\NODf-500iterations.csv")
#results_combined_nodf_prox2 <- subset(results_combined_nodf_prox, Number_combined_sites != 1) 
results_combined_nodf_prox$groups <- factor(results_combined_nodf_prox$groups, levels = c('Low', 'Med.', 'High', 'Random'))



# Plots of web metric response with accumulating sites
ggplot(data=results_combined_nodf_prox, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                             fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+ theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Nestedness (weighted NODF)") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)


# Linear models to assess rate of change in NODF according to site landscape proximity classess
# Responses within the range of the data are linear 
Working_Data_trends<-results_combined_nodf_prox
  Working_Data_trends <- Working_Data_trends %>% drop_na()
low <- Working_Data_trends %>% filter(groups =="Low")
med <- Working_Data_trends %>% filter(groups =="Med.")
high <- Working_Data_trends %>% filter(groups =="High")

#  low
low_fit <- lm(mean ~  Number_combined_sites,           data = low)   
summary(low_fit)
low_intercept <- coef(low_fit)[1]
low_slope <- coef(low_fit)[2]

# med
med_fit <- lm(mean ~  Number_combined_sites,           data = med)   
summary(med_fit)
med_intercept <- coef(med_fit)[1]
med_slope <- coef(med_fit)[2]


# high
high_fit <- lm(mean ~  Number_combined_sites,           data = high)   
summary(high_fit)

high_slope <- coef(med_fit)[2]
high_intercept <- coef(high_fit)[1]

low_intercept
low_slope

med_intercept
med_slope

high_intercept





#---------- generality
# 
#(Weighted) mean effective number of LL species per HL species (generality; HL
# species per LL species for vulnerability), weighted by their marginal totals (row 
#sums); see Tylianakis et al. (2007) and Bersier et al. (2002). This is identical
#to exp(“partner diversity”, i.e., simply the Jost (2006)-recommended version of
#        diversity.
a<-1
for (a in 1:max(Pol$site_order)) {
  Pol_subset<-subset(Pol, site_order==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="generality")[1]  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="generality")[1] } else {null_model_vec[c]<- NA} # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low[a,n_sites]<-CI[1]
    CI_hig[a,n_sites]<-CI[2]
    Mean[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results<-data.frame(c(Mean), c(CI_low) , c(CI_hig))
colnames(results) <- c("mean", "CI_high", "CI_low")
results$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=number_of_site_groups) 
results$groups<-rep(c("Low", "Med.","High"), max_sites_to_combine)

# For the random site picks used as the control

a<-1
for (a in 1:max(Pol$random)) {
  Pol_subset<-subset(Pol, random==a)  #  selects sites that fall in a specific class of the site orders  defined by 'a'
  repsXsites<-matrix(ncol = max_sites_to_combine, nrow = n_iter)  #  empty data frame to put in the web metrics for each combination of sites and itterations
  j<-1
  n_sites<-1
  for (j in 1:max_sites_to_combine) {
    i<-1
    for(i in 1:n_iter) {
      site_codes<-sample(unique(Pol_subset$Site), n_sites)#   identify unique site values # Sample from the vector 'site' 1 element.
      pol_webs_subset<-Pol[Pol$Site %in% site_codes,]   #  selecting only those sites defined in site_codes -  i.e. randomly selected ones
      # reducing this to a single value for a pollinator X site giving mean values for N
      web_linear<-pol_webs_subset %>%
        group_by(across(all_of(c('insect','Plant')))) %>%
        summarise_at(vars(N), list(N_ave = mean))
      web<-acast(web_linear, Plant~insect, value.var="N_ave")
      web<-replace(web, web %in% NA, 0)
      value_metric<-networklevel(web, index="generality")[1]  #  calculates actual web value
      null_model_vec <- rep(NA, times=1000)   #  makes a vector to sotre the null model values in
      c<-1
      for (c in 1:1000) {                              
        if (nrow(web)>1 & ncol(web)>1 ){null_model_vec[c]<-networklevel(vaznull(1, web)[[1]], index="generality")[1] } else {null_model_vec[c]<- NA} # calculates a null model and the relvent metric   
      }
      null_mean<-mean(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      null_sd<-sd(null_model_vec,  na.rm=TRUE)  #  calculates a mean for each of the null model values
      z_score<- (value_metric- null_mean)/null_sd   # z-score for the metric (observed metric-mean metric [null model] divided by standard deviation [null model])
      repsXsites[i,n_sites]<- z_score
    } 
    CI<-confidence_interval(repsXsites[,n_sites], 0.95)    #  95% confidence interval see above function
    CI_low_rand[a,n_sites]<-CI[1]
    CI_hig_rand[a,n_sites]<-CI[2]
    Mean_rand[a,n_sites]<-mean(repsXsites[,n_sites])
    n_sites<-n_sites+1
  }   
}

results_rand<-data.frame(c(Mean_rand), c(CI_low_rand) , c(CI_hig_rand))
colnames(results_rand) <- c("mean", "CI_high", "CI_low")
results_rand$Number_combined_sites<-rep(1:max_sites_to_combine, times=1, each=1) 
results_rand$groups<-rep(c("Random"), max_sites_to_combine)

results_combined<-rbind(results, results_rand)

results_combined$groups <- factor(results_combined$groups, levels = c('Low', 'Med.', 'High', 'Random'))

write.csv(results_combined, "results\\Grass\\Proximity1km\\generality-500iterations.csv")
results_combined_gen_prox<-read.csv( "results\\Grass\\Proximity1km\\generality-500iterations.csv")
#results_combined_gen_prox2 <- subset(results_combined_gen_prox, Number_combined_sites != 1) 
results_combined_gen_prox$groups <- factor(results_combined_gen_prox$groups, levels = c('Low', 'Med.', 'High', 'Random'))  



# Plots of web metric response with accumulating sites
ggplot(data=results_combined_gen_prox, aes(x=as.factor(Number_combined_sites), y=mean, group=groups,
                                            fill=groups)) +
  geom_line(aes(linetype=groups, color=groups))+
  geom_point(aes(color=groups))+ theme_bw() +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=.3, linetype=0) +
  theme(legend.text = element_text(size=10),legend.position="top",  legend.justification = "center", legend.box.just = "center", legend.title=element_blank()) +
  labs(x = "Number of sites combined") +
  labs(y = "Generality") + 
  theme(axis.text=element_text(size=15)) +
  theme(axis.title=element_text(size=15)) +
  expand_limits(y=0)


# Linear models to assess rate of change in generality according to site landscape proximity classess
# Responses within the range of the data are linear 
Working_Data_trends<-results_combined_gen_prox
  Working_Data_trends <- Working_Data_trends %>% drop_na()
low <- Working_Data_trends %>% filter(groups =="Low")
med <- Working_Data_trends %>% filter(groups =="Med.")
high <- Working_Data_trends %>% filter(groups =="High")

#  low
low_fit <- lm(mean ~  Number_combined_sites,           data = low)   
summary(low_fit)
low_intercept <- coef(low_fit)[1]
low_slope <- coef(low_fit)[2]

# med
med_fit <- lm(mean ~  Number_combined_sites,           data = med)   
summary(med_fit)
med_intercept <- coef(med_fit)[1]
med_slope <- coef(low_fit)[2]

# high
high_fit <- lm(mean ~  Number_combined_sites,           data = high)   
summary(high_fit)
high_intercept <- coef(high_fit)[1]
high_slope <- coef(high_fit)[2]


low_intercept
low_slope

med_intercept
med_slope

high_intercept
high_slope
