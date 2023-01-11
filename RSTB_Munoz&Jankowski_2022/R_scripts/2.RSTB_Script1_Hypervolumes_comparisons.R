#######################################################################################
### Manuscript_Flocks in a community context                                         ###
### Hypervolume comparisons                                                                                 ###
### R-code                                                                          ###
### Jenny Munoz                                                                     ###
### Last update: June 02 2022                                                       ###
################################################################################

###################################
## Section 1 Initial set up------------------------------------------------------------------
###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_

#Plots
install.packages("graphics")
install.packages("ggplot2")
install.packages("plotly")
install.packages ("ggpubr")
install.packages ("grid")
install.packages("gridExtra")
install.packages("rgl")


# Modeling
install.packages("vegan")
install.packages("devtools")
install.packages("lme4")
install.packages("stats")
install.packages("nlme")
install.packages("mgcv")
install.packages("Matrix")
install.packages("lattice")
install.packages("betareg")
install.packages("AICcmodavg")
install.packages("MASS")
install.packages("visreg")
install.packages("emmeans")
install.packages("car") #Anova command
#PCA
install.packages("TPD")
install.packages("factoextra")
#traits
install.packages("FD")
install.packages("geometry")
#hypervolumes
install.packages("terra")
install.packages("hypervolume")
#phylogenetis
install.packages("ape")

#overlap 
install.packages("reshape") 
#Data
install.packages("tidyverse") #I added this libraries for easier manipulation of 
install.packages("tibble")
install.packages ("data.table")
install.packages ("extrafont")
install.packages("assertr")
install.packages("convertr")
install.packages("purrr")
install.packages("daisy")
#install.packages("beep")

# To find outliers
install.packages("assertr")
install.packages("stringdist")
install.packages("skimr")

#clusters
install.packages("dendextend")

# plots 

### ### ###
#Libraries 
### ### ### 
#Plots
library(stats) 
library(graphics)
library(ggplot2)
library(plotly)
library(ggpubr)
library(grid)
library(gridExtra)
library(rgl)

#Modeling 
library(vegan)
library(devtools)
library(lme4)
library(stats)
library(nlme)
library(mgcv)
library(Matrix)
library(lattice)
library(eHOF) 
library(betareg)
library(AICcmodavg)
library(MASS)
library(visreg)
library(emmeans)
library(car) #Anova command
#PCA
library(TPD)
library(factoextra)
#hypervolumes 
library(terra)
library(hypervolume)

#traits
library(FD)
library(geometry)

#phylogenetis

library(ape)

# overlap
library(reshape)
#Data
library(tidyverse)   
library(data.table)
library(extrafont)
library(assertr)
library(convertr)
library(purrr)
library(tibble)
library(skimr)

# NOTES ON BANDWITH 

#How do I choose the bandwidth parameter?
#There is no objective way to choose the bandwidth parameter. 
#You can use the provided estimate_bandwidth function to try one possibility that trades off between variance in the data and sample size. 
#However this Silverman estimator is only optimal for univariate normal data and has unknown properties when used elsewhere. 
#In particular, this estimator is not guaranteed to minimize the mean integrated square error.
#Another option is to use a fixed value for the analysis that reflects your understanding of uncertainty in the data. 
#There are two potential caveats. 
#First, you may choose a value so low that all points appear to be disjunct (a value of @DisjunctFactor that approaches one).
#This is especially likely in high dimensionality analyses.
#You probably need to use a fixed value that is higher than you expect. 
#Second, your results may be sensitive to the particular value chosen, especially for analyses of negative features where the appearance of hypervolume 'holes' depends on the padding put around each data point.
#In this case you should repeat your analysis for a range of bandwidth values and determine if the qualitative conclusions are robust to the bandwidth choice.


#Choose the dimensionality of the analysis. Choose as low a value as possible that is consistent with your analytic goals. Subset the data to only these dimensions.
# all axes to a common and comparable scale, e.g. by z- or log- transformation.
#Choose a kernel bandwidth. Use either a fixed value or use estimate_bandwidth to determine a value.
#To ensure results are not sensitive to analytical choices, repeat analysis with other fixed bandwidth values and determine whether conclusions are qualitatively different.
#To ensure results are not sensitive to sample size, repeat analysis in the context of a null model using a dataset with the same number of observations but random values.
#Report dimensionality and bandwidth choice for analysis, as well as quantile threshold and repsperpoint if non-default values were used.


# Part 1 [ Observed hypervolumes] and comparison **Community vs flocks 0.60 fixed/ -----------------------------------
#To ensure results are not sensitive to analytical choices, repeat analysis with other fixed bandwidth values and determine whether conclusions are qualitatively different.
# After readingthe FAQ of Benjamin Bloder (the author of the package) https://www.benjaminblonder.org/hypervolume_faq.html
# Some relevant notes 
#You probably need to use a fixed value that is higher than you expect. 
#Second, your results may be sensitive to the particular value chosen, especially for analyses of negative features where the appearance of hypervolume 'holes' depends on the padding put around each data point.
#In this case you should repeat your analysis for a range of bandwidth values and determine if the qualitative conclusions are robust to the bandwidth choice.

# Also inote when comparing hypervokums
#Some care should also be taken when comparing the volumes of hypervolumes of the same dimensionality if a fixed kernel bandwidth was used to construct them. 
#In this case, the volume of the hypervolume is approximately linearly proportional to the number of observations in the dataset.
#This is because each new data point contributes approximately the same amount of volume, unless it overlaps with a previous data point. 
#The issue is then that the largest hypervolumes will be those constructed from the largest number of data points. 
#This may reflect the true structure of your data, but if it does not, you should instead proceed with a null-modeling procedure
#where you compare the observed hypervolume to that of a distribution of null hypervolumes constructed by resampling an identical number of observations.
#Instead of reporting a raw hypervolume you can report a deviation hypervolume, e.g. a z-score.

### kde 0.6 ( Which is larger than the one automatically choosen by silverman in the commmunity group 5.3 was the max as per suggestion "You probably need to use a fixed value that is higher than you expect. )
#Hypervolume nonsocial
set.seed(17)
obs_community_hypervolume_gaussian_3<-hypervolume_gaussian(bird_traits_z_community, 
                                                           name="Community_2000_samples_060",
                                                           samples.per.point =2000,
                                                           weight = NULL,
                                                           sd.count= 3,
                                                           quantile.requested = 0.95,
                                                           quantile.requested.type = "probability",
                                                           kde.bandwidth= estimate_bandwidth(bird_traits_z_community, method="fixed", value = c(0.60,0.60,0.60,0.60,0.60,0.60,0.60))) 


saveRDS(obs_community_hypervolume_gaussian_3, "data/data_analyses/stats_analyses/obs_community_hypervolume_gaussian_2000s_bw060.rds") #SAVED
obs_community_hypervolume_gaussian_3<-readRDS("data/data_analyses/stats_analyses/obs_community_hypervolume_gaussian_2000s_bw060.rds")


set.seed(17)
obs_nonsocial_hypervolume_gaussian_3<-hypervolume_gaussian(bird_traits_z_nonsocial, 
                                                          name="non social_2000_samples_060",
                                                          samples.per.point =2000,
                                                          weight = NULL,
                                                          sd.count= 3,
                                                          quantile.requested = 0.95,
                                                          quantile.requested.type = "probability",
                                                          kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed",  value = c(0.60,0.60,0.60,0.60,0.60,0.60,0.6))
)
saveRDS(obs_nonsocial_hypervolume_gaussian_3, "data/data_analyses/stats_analyses/obs_nonsocial_hypervolume_gaussian_2000s_bw060.rds") #SAVED
obs_nonsocial_hypervolume_gaussian_3<-readRDS("data/data_analyses/stats_analyses/obs_nonsocial_hypervolume_gaussian_2000s_bw060.rds")


#Hypervolume social
set.seed(17) #
obs_social_hypervolume_gaussian_3<-hypervolume_gaussian(bird_traits_z_social, 
                                                        name="All flocks_2000_samples_060",
                                                        samples.per.point =2000,
                                                        weight = NULL,
                                                        sd.count= 3,
                                                        quantile.requested = 0.95,
                                                        quantile.requested.type = "probability",
                                                        kde.bandwidth= estimate_bandwidth(bird_traits_z_social, method="fixed",  value = c(0.60,0.60,0.60,0.60,0.60,0.60,0.60))) 


saveRDS(obs_social_hypervolume_gaussian_3, "data/data_analyses/stats_analyses/obs_social_hypervolume_gaussian_2000s_bw060.rds") #SAVED
obs_social_hypervolume_gaussian_3<-readRDS("data/data_analyses/stats_analyses/obs_social_hypervolume_gaussian_2000s_bw060.rds")


### 2 #### set raw hypervolumes for comparison PAIRWISE

#The intersection is the points in both hypervolumes, 
#the union those in either hypervolume, and the
#unique components the points in one hypervolume but not the other.

set.seed(17)
social_nonsocial_set<-hypervolume_set(obs_nonsocial_hypervolume_gaussian_3,obs_social_hypervolume_gaussian_3,check.memory = FALSE)
social_nonsocial_overlap<-hypervolume_overlap_statistics(social_nonsocial_set)

summary(social_nonsocial_overlap)
social_nonsocial_set

# I am not sure this ones are neecesary given that they social and community overlap so much 
social_community_set<-hypervolume_set(obs_community_hypervolume_gaussian_3,obs_social_hypervolume_gaussian_3,check.memory = FALSE)
social_community_overlap<-hypervolume_overlap_statistics(social_community_set)
summary(social_community_set)
social_community_set

nonsocial_community_set<-hypervolume_set(obs_community_hypervolume_gaussian_3,obs_nonsocial_hypervolume_gaussian_3,check.memory = FALSE)
nonsocial_community_overlap<-hypervolume_overlap_statistics(nonsocial_community_set)
summary(social_community_set)
social_community_set


a=417.515741/16969.067823
b=10705.620759/16969.067823
c=5845.93/16969.067823
#jaccard 0.3445052 
a
b
c
a+b+c

6863.651980/18175.85 # the overlaping of the three comunity social and non social

# lets do a group comparison to see what is teh overalap of all of them ... shoudl be consistent with teh overlap in pairwise (social no social)
hv_list_community_social_nonsocial=hypervolume_join(obs_community_hypervolume_gaussian_3,obs_social_hypervolume_gaussian_3,obs_nonsocial_hypervolume_gaussian_3)
community_social_nonsocial_intersection<-hypervolume_set_n_intersection(hv_list_community_social_nonsocial, num.points.max = NULL,verbose = TRUE, distance.factor = 1, check.hyperplane = FALSE)
saveRDS(community_social_nonsocial_intersection, "data/data_analyses/stats_analyses/community_social_nonsocial_intersection.rds")
community_social_nonsocial_intersection

# THIS TWO ARE NOT NEEDED
hv_list_community_social=hypervolume_join(obs_community_hypervolume_gaussian_3,obs_social_hypervolume_gaussian_3)
community_social_intersection<-hypervolume_set_n_intersection(hv_list_community_social, num.points.max = NULL,verbose = TRUE, distance.factor = 1, check.hyperplane = FALSE)
community_social_intersection

hv_list_social_nonsocial=hypervolume_join(obs_nonsocial_hypervolume_gaussian_3,obs_social_hypervolume_gaussian_3)
social_nonsocial_intersection<-hypervolume_set_n_intersection(hv_list_social_nonsocial, num.points.max = NULL,verbose = TRUE, distance.factor = 1, check.hyperplane = FALSE)
social_nosocial_intersection

# Part 2 [Randomized hypervolumes] matching sample size **Community vs flocks_bw0.60 ---------------------------
# Hypervolume with matching sample size (minimum diversity in a flocking assemblage 269)
#Commmunity
set.seed(17)
com_hypervolumes_randomized_community_269=hypervolume_resample("269_community_resample_bw", obs_community_hypervolume_gaussian_3,
                                                               method="bootstrap",
                                                               n = 50, 
                                                               points_per_resample = "269", 
                                                               cores = 6,
                                                               verbose = TRUE, 
                                                               weight_func = NULL)

com_hypervolumes_randomized_community_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_community_resample_bw"

#Flocking
set.seed(17)
com_hypervolumes_randomized_flocking_269=hypervolume_resample("269_flocking_resample_bw", obs_social_hypervolume_gaussian_3,
                                                               method="bootstrap",
                                                               n = 50, 
                                                               points_per_resample = "269", 
                                                               cores = 6,
                                                               verbose = TRUE,
                                                               weight_func = NULL)
com_hypervolumes_randomized_flocking_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_flocking_resample_bw"

#No-flocking 
set.seed(18)
com_hypervolumes_randomized_nonflocking_269=hypervolume_resample("269_nonflocking_resample_bw",  obs_nonsocial_hypervolume_gaussian_3,
                                                              method="bootstrap",
                                                              n = 50, 
                                                              points_per_resample = "269", 
                                                              cores = 6,
                                                              verbose = TRUE,
                                                              weight_func = NULL)






com_hypervolumes_randomized_nonflocking_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_nonflocking_resample_bw"


# List of matched assemblages 

com_hypervolumes_randomized_nonflocking_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_nonflocking_resample_bw"
com_hypervolumes_randomized_community_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_community_resample_bw"
com_hypervolumes_randomized_flocking_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_flocking_resample_bw"


# Part 3 [ Observed hypervolumes] ***flocks across the gradient ~fixed bandwidth 0.60 -----------------------------------------------------
set.seed(17)
obs_lowlands_hypervolume_gaussian_3<-hypervolume_gaussian(bird_traits_z_und_lowlands, 
                                                          name="lowlands_2000_s_bw065",
                                                          samples.per.point =2000,
                                                          weight = NULL,
                                                          sd.count= 3,
                                                          quantile.requested = 0.95,
                                                          quantile.requested.type = "probability",
                                                          kde.bandwidth= estimate_bandwidth(bird_traits_z_und_lowlands, method="fixed", value = c(0.60,0.60,0.60,0.60,0.60,0.60,0.60))
                                                          #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.60,0.60,0.60,0.60,0.60,0.60,0.60))
)

saveRDS(obs_lowlands_hypervolume_gaussian_3, "data/data_analyses/stats_analyses/obs_lowlands_hypervolume_gaussian_2000_s_bw065.rds") 
obs_lowlands_hypervolume_gaussian_3<-readRDS("data/data_analyses/stats_analyses/obs_lowlands_hypervolume_gaussian_2000_s_bw065.rds")

# create the hypervolume for canopy flocks ( see than here I a keeping the number of samples per point constant )
set.seed(17)
obs_canopy_hypervolume_gaussian_3<-hypervolume_gaussian(bird_traits_z_canopy, 
                                                        name="canopy_2000_s_bw065",
                                                        samples.per.point=2000,
                                                        weight = NULL,
                                                        sd.count= 3,
                                                        quantile.requested = 0.95,
                                                        quantile.requested.type = "probability",
                                                        kde.bandwidth= estimate_bandwidth(bird_traits_z_canopy, method="fixed", value = c(0.60,0.60,0.60,0.60,0.60,0.60,0.60))
                                                        #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.60,0.60,0.60,0.60,0.60,0.60,0.60))
)

saveRDS(obs_canopy_hypervolume_gaussian_3, "data/data_analyses/stats_analyses/obs_canopy_hypervolume_gaussian_2000_s_bw065.rds") 
obs_canopy_hypervolume_gaussian_3<-readRDS("data/data_analyses/stats_analyses/obs_canopy_hypervolume_gaussian_2000_s_bw065.rds")

set.seed(17)
obs_bamboo_hypervolume_gaussian_3<-hypervolume_gaussian(bird_traits_z_bamboo, 
                                                        name="bamboo_2000_s_bw065",
                                                        samples.per.point =2000,
                                                        weight = NULL,
                                                        sd.count= 3,
                                                        quantile.requested = 0.95,
                                                        quantile.requested.type = "probability",
                                                        kde.bandwidth= estimate_bandwidth(bird_traits_z_bamboo, method="fixed", value = c(0.60,0.60,0.60,0.60,0.60,0.60,0.60))
)

saveRDS(obs_bamboo_hypervolume_gaussian_3, "data/data_analyses/stats_analyses/obs_bamboo_hypervolume_gaussian_2000_s_bw065.rds") 
obs_bamboo_hypervolume_gaussian_3<-readRDS("data/data_analyses/stats_analyses/obs_bamboo_hypervolume_gaussian_2000_s_bw065.rds")

set.seed(17)
obs_montane_hypervolume_gaussian_3<-hypervolume_gaussian(bird_traits_z_andes, 
                                                         name="montane_2000_s_bw065",
                                                         samples.per.point =2000,
                                                         weight = NULL,
                                                         sd.count= 3,
                                                         quantile.requested = 0.95,
                                                         quantile.requested.type = "probability",
                                                         #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5))
                                                         kde.bandwidth= estimate_bandwidth(bird_traits_z_andes, method="fixed", value = c(0.60,0.60,0.60,0.60,0.60,0.60,0.60))
)

saveRDS(obs_montane_hypervolume_gaussian_3, "data/data_analyses/stats_analyses/obs_montane_hypervolume_gaussian_2000_s_bw065.rds") 
obs_montane_hypervolume_gaussian_3<-readRDS("data/data_analyses/stats_analyses/obs_montane_hypervolume_gaussian_2000_s_bw065.rds")

set.seed(17)
obs_h_montane_hypervolume_gaussian_3<-hypervolume_gaussian(bird_traits_z_h_andes, 
                                                           name="h_montane_2000_s_bw065",
                                                           samples.per.point =2000,
                                                           weight = NULL,
                                                           sd.count= 3,
                                                           quantile.requested = 0.95,
                                                           quantile.requested.type = "probability",
                                                           #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5))
                                                           kde.bandwidth= estimate_bandwidth(bird_traits_z_h_andes, method="fixed", value = c(0.60,0.60,0.60,0.60,0.60,0.60,0.60))
)
saveRDS(obs_h_montane_hypervolume_gaussian_3, "data/data_analyses/stats_analyses/obs_h_montane_hypervolume_gaussian_2000_s_bw065.rds") 
obs_h_montane_hypervolume_gaussian_3<-readRDS("data/data_analyses/stats_analyses/obs_h_montane_hypervolume_gaussian_2000_s_bw065.rds") 


###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_

# Hypervolume [Observed]with matching sample***flocks across the gradient size bw 06  --------


#Terrafirme
# make sure you keep the sample size identical to what you want to compare with 
com_hypervolumes_randomized_terrafirme_43_bw=hypervolume_resample("43_terrafirme_resample_bw", obs_lowlands_hypervolume_gaussian_3,
                                                               method="bootstrap",
                                                               n = 50, 
                                                               points_per_resample ="43", 
                                                               cores = 6,
                                                               verbose = TRUE, 
                                                               weight_func = NULL)

com_hypervolumes_randomized_terrafirme_43_bw="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_terrafirme_resample_bw"

#Canopy

com_hypervolumes_randomized_canopy_43_bw=hypervolume_resample("43_canopy_resample_bw",  obs_canopy_hypervolume_gaussian_3,
                                                           method="bootstrap", 
                                                           n = 50, 
                                                           points_per_resample =43,  
                                                           cores = 6,
                                                           verbose = TRUE, 
                                                           weight_func = NULL)



com_hypervolumes_randomized_canopy_43_bw="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_canopy_resample_bw"

#Bamboo

com_hypervolumes_randomized_bamboo_43_bw=hypervolume_resample("43_bamboo_resample_bw", obs_bamboo_hypervolume_gaussian_3,
                                                           method="bootstrap", 
                                                           n = 50, 
                                                           points_per_resample="sample_size", 
                                                           cores = 6,
                                                           verbose = TRUE, 
                                                           weight_func = NULL)

com_hypervolumes_randomized_bamboo_43_bw="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_bamboo_resample_bw"



#Montane
com_hypervolumes_randomized_montane_43_bw=hypervolume_resample("43_montane_resample_bw",obs_montane_hypervolume_gaussian_3,
                                                                method="bootstrap", 
                                                                n = 50, 
                                                                points_per_resample =43,  # make sure you keep the sample size identical to what you want to compare with 
                                                                cores = 6,
                                                                verbose = TRUE, 
                                                                weight_func = NULL)

com_hypervolumes_randomized_montane_43_bw="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_montane_resample_bw"



#High montane


com_hypervolumes_randomized_hmontane_43_bw=hypervolume_resample("43_hmontane_resample_bw", obs_h_montane_hypervolume_gaussian_3,
                                                             method="bootstrap", 
                                                             n = 50, 
                                                             points_per_resample =43,  # make sure you keep the sample size identical to what you want to compare with 
                                                             cores = 6,
                                                             verbose = TRUE, 
                                                             weight_func = NULL)
com_hypervolumes_randomized_hmontane_43_bw="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_hmontane_resample_bw"





#  Stats overlap [ observed ] Hypervolume ***flocks across the gradient------------------------------------

###_###_###_##_
#Hypervolume overlaps pairwise
###_###_###_##_

###_###_###_##_
#Hypervolume overlaps  between pairs of flock types or community vs flocks
###_###_###_##_

# set raw hypervolumes for comparison

set.seed(17)

understory_canopy_set<-hypervolume_set(obs_lowlands_hypervolume_gaussian_3,obs_canopy_hypervolume_gaussian_3, check.memory = FALSE)
understory_canopy_overlap<-hypervolume_overlap_statistics(understory_canopy_set)
summary(understory_canopy_overlap)
understory_canopy_set

lowlands_andes_set<-hypervolume_set(obs_lowlands_hypervolume_gaussian_3,obs_montane_hypervolume_gaussian_3, check.memory = FALSE)
hypervolume_overlap_statistics(lowlands_andes_set)
summary(lowlands_andes_set)

lowlands_h_andes_set<-hypervolume_set(obs_lowlands_hypervolume_gaussian_3,obs_h_montane_hypervolume_gaussian_3,check.memory = FALSE)
hypervolume_overlap_statistics(lowlands_h_andes_set)
summary(lowlands_h_andes_set)

understory_bamboo_set<-hypervolume_set(obs_lowlands_hypervolume_gaussian_3,obs_bamboo_hypervolume_gaussian_3, check.memory = FALSE)
hypervolume_overlap_statistics(understory_bamboo_set)
summary(understory_bamboo_set)
#
bamboo_canopy_set<-hypervolume_set(obs_bamboo_hypervolume_gaussian_3,obs_canopy_hypervolume_gaussian_3, check.memory = FALSE)
bamboo_canopy_overlap<-hypervolume_overlap_statistics(bamboo_canopy_set)
summary(bamboo_canopy_set)

canopy_andes_set<-hypervolume_set(obs_montane_hypervolume_gaussian_3,obs_canopy_hypervolume_gaussian_3,check.memory = FALSE)
hypervolume_overlap_statistics(canopy_andes_set)
summary(canopy_andes_set)

canopy_h_andes_set<-hypervolume_set(obs_h_montane_hypervolume_gaussian_3,obs_canopy_hypervolume_gaussian_3,check.memory = FALSE)
hypervolume_overlap_statistics(canopy_h_andes_set)
summary(canopy_andes_set)
#
bamboo_montane_set<-hypervolume_set(obs_bamboo_hypervolume_gaussian_3,obs_montane_hypervolume_gaussian_3, check.memory = FALSE)
hypervolume_overlap_statistics(bamboo_montane_set) # double check this values on the table
summary(bamboo_montane_set)

bamboo_h_andes_set<-hypervolume_set(obs_bamboo_hypervolume_gaussian_3,obs_h_montane_hypervolume_gaussian_3, check.memory = FALSE)
hypervolume_overlap_statistics(bamboo_h_andes_set)
summary(bamboo_h_andes_set)

montane_h_andes_set<-hypervolume_set(obs_montane_hypervolume_gaussian_3,obs_h_montane_hypervolume_gaussian_3, check.memory = FALSE)
montane_h_andes_overlap<-hypervolume_overlap_statistics(montane_h_andes_set)
summary(montane_h_andes_set)



###_###_###_##_
#Hypervolume overlaps between all flock types
###_###_###_##_
hv_list_flocks_gradient=hypervolume_join(obs_bamboo_hypervolume_gaussian_3,obs_canopy_hypervolume_gaussian_3,obs_lowlands_hypervolume_gaussian_3,obs_montane_hypervolume_gaussian_3,obs_h_montane_hypervolume_gaussian_3)
raw_all_flocks_intersection<-hypervolume_set_n_intersection(hv_list_flocks_gradient, num.points.max = NULL,verbose = TRUE, distance.factor = 1, check.hyperplane = FALSE)

saveRDS(raw_all_flocks_intersection, "data/data_analyses/stats_analyses/raw_all_flocks_intersection.rds") 


2249.776185/6263.44 # # percentage intersection between observed flock types



# [RANDOMIZED HYPERVOLUMENS]  STATS ** Community and ***flocks across the gradient ---------------------------------------------

com_hypervolumes_randomized_community_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_community_resample_bw"
hvs_hypervolumes_randomized_community_269=to_hv_list(com_hypervolumes_randomized_community_269)
volumes_community_269<-get_volume(hvs_hypervolumes_randomized_community_269)
mean(volumes_community_269)
sd(volumes_community_269)

com_hypervolumes_randomized_nonflocking_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_nonflocking_resample_bw"
hvs_hypervolumes_randomized_nonflocking_269=to_hv_list(com_hypervolumes_randomized_nonflocking_269)
volumes_nonflocking_269<-get_volume(hvs_hypervolumes_randomized_nonflocking_269)
mean(volumes_nonflocking_269)
sd(volumes_nonflocking_269)

com_hypervolumes_randomized_flocking_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_flocking_resample_bw"
hvs_hypervolumes_randomized_flocking_269=to_hv_list(com_hypervolumes_randomized_flocking_269)
volumes_flocking_269<-get_volume(hvs_hypervolumes_randomized_flocking_269)
mean(volumes_flocking_269)
sd(volumes_flocking_269)

# For flocks assemblages 

com_hypervolumes_randomized_bamboo_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_bamboo_resample_bw"
hvs_randomized_bamboo_43=to_hv_list(com_hypervolumes_randomized_bamboo_43)
volumes_bamboo_43<-get_volume(hvs_randomized_bamboo_43)
mean(volumes_bamboo_43)
sd(volumes_bamboo_43)
    
com_hypervolumes_randomized_terrafirme_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_terrafirme_resample_bw"
hvs_randomized_terrafirme_43=to_hv_list (com_hypervolumes_randomized_terrafirme_43)
volumes_terrafirme_43<-get_volume(hvs_randomized_terrafirme_43)
mean(volumes_terrafirme_43)
sd(volumes_terrafirme_43)  

com_hypervolumes_randomized_canopy_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_canopy_resample_bw"
hvs_randomized_canopy_43=to_hv_list(com_hypervolumes_randomized_canopy_43)
volumes_canopy_43<-get_volume(hvs_randomized_canopy_43)
mean(volumes_canopy_43)
sd(volumes_canopy_43)
    
com_hypervolumes_randomized_montane_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_montane_resample_bw"
hvs_randomized_montane_43=to_hv_list(com_hypervolumes_randomized_montane_43)
volumes_montane_43<-get_volume(hvs_randomized_montane_43)
mean(volumes_montane_43)
sd(volumes_montane_43)
    
com_hypervolumes_randomized_hmontane_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_hmontane_resample_bw"
hvs_randomized_h_montane_43=to_hv_list( com_hypervolumes_randomized_hmontane_43)
volumes_h_montane_43<-get_volume(hvs_randomized_h_montane_43)
mean(volumes_h_montane_43)
sd(volumes_h_montane_43)



####_###_###_## Hypervolume STATS OVERLAP ** Community and ***flocks across the gradient-----------------------------------------------

#[RANDOMIZED HYPERVOLUMENS]
#Creating the 50*50 hypervolumes is too much computationally, so I will create just groups of 10*10

# List of matched flocks reduced to 20 for the overlap computation 

com_hypervolumes_randomized_nonflocking_269_10="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_nonflocking_resample_bw_10"
#com_hypervolumes_randomized_community_269_10="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_community_resample_20"
com_hypervolumes_randomized_flocking_269_10="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_flocking_resample_bw_10"

4+5
# not used 
statistics_flocking_nonflocking269_bw=hypervolume_overlap_confidence(com_hypervolumes_randomized_flocking_269_10, com_hypervolumes_randomized_nonflocking_269_10,CI = .95, cores = 7) #43
saveRDS(statistics_flocking_nonflocking269_bw, "data/data_analyses/stats_analyses/statistics_flocking_nonflocking269_bw.rds") 
mean(statistics_flocking_nonflocking269_bw$distribution [,1])
sd(statistics_flocking_nonflocking269_bw$distribution [,1])
mean (statistics_flocking_nonflocking269_bw$distribution [,3]) # fraction unique flocking
sd (statistics_flocking_nonflocking269_bw$distribution [,3])
mean (statistics_flocking_nonflocking269_bw$distribution [,4])# zfraction unique no flocking
sd (statistics_flocking_nonflocking269_bw$distribution [,4])

#need to run these three estimation 8 hours each
statistics_community_nonflocking_269=hypervolume_overlap_confidence(com_hypervolumes_randomized_community_269_20, com_hypervolumes_randomized_nonflocking_269_20,CI = .95, cores = 7) #43
saveRDS(statistics_community_nonflocking_269, "data/data_analyses/statistics_community_nonflocking_269.rds") 
mean(statistics_community_nonflocking_269$distribution [,1])
sd(statistics_community_nonflocking_269$distribution [,1])
mean (statistics_community_nonflocking_269$distribution [,3]) # fraction unique terrafirme
sd (statistics_community_nonflocking_269$distribution [,3])
mean (statistics_community_nonflocking_269$distribution [,4]) # fraction unique canopy
sd (statistics_community_nonflocking_269$distribution [,4])

# only this one presented in te manuscriot to avoid redundance 

statistics_community_flocking_269=hypervolume_overlap_confidence(com_hypervolumes_randomized_community_269_20, com_hypervolumes_randomized_flocking_269_20,CI = .95, cores = 7) #43
saveRDS(statistics_community_flocking_269, "data/data_analyses/statistics_community_flocking_269.rds") 
mean(statistics_community_flocking_269$distribution [,1])
sd(statistics_community_flocking_269$distribution [,1])
mean (statistics_community_flocking_269$distribution [,3]) # fraction unique community
sd (statistics_community_flocking_269$distribution [,3])
mean (statistics_community_flocking_269$distribution [,4]) # fraction unique flocks
sd (statistics_community_flocking_269$distribution [,4])

###_###_###_###_###_###_###_###_###_###_###_###_###_###_
#Hypervolume STATS OVERLAP ** Community and ***flocks across the gradient
###_###_###_###_###_###_###_###_###_###_###_###_###_###_


com_hypervolumes_randomized_terrafirme_43_10="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_terrafirme_resample_bw_10"

com_hypervolumes_randomized_bamboo_43_10="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_bamboo_resample_bw_10"

com_hypervolumes_randomized_canopy_43_10="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_canopy_resample_bw_10"

com_hypervolumes_randomized_montane_43_10="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_montane_resample_bw_10"

com_hypervolumes_randomized_hmontane_43_10="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_hmontane_resample_bw_10"

statistics_terrafirme_canopy_43_bw=hypervolume_overlap_confidence(com_hypervolumes_randomized_terrafirme_43_10, com_hypervolumes_randomized_canopy_43_10,CI = .95, cores = 7) #43
saveRDS(statistics_terrafirme_canopy_43_bw, "data/data_analyses/stats_analyses/statistics_terrafirme_canopy_43_bw.rds") 
statistics_terrafirme_canopy_43_bw<-readRDS("data/data_analyses/stats_analyses/statistics_terrafirme_canopy_43_bw.rds")
mean(statistics_terrafirme_canopy_43_bw$distribution [,1])
sd(statistics_terrafirme_canopy_43_bw$distribution [,1])
mean (statistics_terrafirme_canopy_43_bw$distribution [,3]) # fraction unique terrafirme
sd (statistics_terrafirme_canopy_43_bw$distribution [,3])
mean (statistics_terrafirme_canopy_43_bw$distribution [,4]) # fraction unique canopy
sd (statistics_terrafirme_canopy_43_bw$distribution [,4])


statistics_terrafirme_bamboo_43_bw=hypervolume_overlap_confidence(com_hypervolumes_randomized_terrafirme_43_10,com_hypervolumes_randomized_bamboo_43_10,CI = .95, cores = 7) # 43
saveRDS(statistics_terrafirme_bamboo_43_bw, "data/data_analyses/stats_analyses/statistics_terrafirme_bamboo_43_bw.rds") 
statistics_terrafirme_bamboo_43_bw<-readRDS("data/data_analyses/stats_analyses/statistics_terrafirme_bamboo_43_bw.rds")
mean(statistics_terrafirme_bamboo_43_bw$distribution [,1])
sd(statistics_terrafirme_bamboo_43_bw$distribution [,1])
mean (statistics_terrafirme_bamboo_43_bw$distribution [,3]) # fraction unique terrafirme
sd (statistics_terrafirme_bamboo_43_bw$distribution [,3])
mean (statistics_terrafirme_bamboo_43_bw$distribution [,4]) # fraction unique bamboo
sd (statistics_terrafirme_bamboo_43_bw$distribution [,4])


statistics_terrafirme_montane_43_bw=hypervolume_overlap_confidence(com_hypervolumes_randomized_terrafirme_43_10,com_hypervolumes_randomized_montane_43_10,CI = .95, cores = 7) # 43
saveRDS(statistics_terrafirme_montane_43_bw,"data/data_analyses/stats_analyses/statistics_terrafirme_montane_43_bw.rds")
statistics_bamboo_montane_43_bw<-readRDS("data/data_analyses/stats_analyses/statistics_terrafirme_montane_43_bw.rds")
mean(statistics_terrafirme_montane_43_bw$distribution [,1])
sd(statistics_terrafirme_montane_43_bw$distribution [,1])
mean (statistics_terrafirme_montane_43_bw$distribution [,3]) # fraction unique terrafirme
sd (statistics_terrafirme_montane_43_bw$distribution [,3])
mean (statistics_terrafirme_montane_43_bw$distribution [,4]) # fraction unique montane
sd (statistics_terrafirme_montane_43_bw$distribution [,4])


statistics_terrafirme_h_montane_43_bw=hypervolume_overlap_confidence(com_hypervolumes_randomized_terrafirme_43_10,com_hypervolumes_randomized_hmontane_43_10,CI = .95, cores = 7) # 43
saveRDS(statistics_terrafirme_h_montane_43_bw,"data/data_analyses/stats_analyses/statistics_terrafirme_h_montane_43_bw.rds")
statistics_terrafirme_h_montane_43_bw<-readRDS("data/data_analyses/stats_analyses/statistics_terrafirme_h_montane_43_bw.rds")
mean(statistics_terrafirme_h_montane_43_bw$distribution [,1])
sd(statistics_terrafirme_h_montane_43_bw$distribution [,1])
mean (statistics_terrafirme_h_montane_43_bw$distribution [,3]) # fraction unique terrafirme
sd (statistics_terrafirme_h_montane_43_bw$distribution [,3])
mean (statistics_terrafirme_h_montane_43_bw$distribution [,4]) # fraction unique h_montane
sd (statistics_terrafirme_h_montane_43_bw$distribution [,4])


statistics_bamboo_canopy_43_bw=hypervolume_overlap_confidence(com_hypervolumes_randomized_bamboo_43_10, com_hypervolumes_randomized_canopy_43_10, CI = .95, cores = 7) # 43
saveRDS(statistics_bamboo_canopy_43_bw,"data/data_analyses/stats_analyses/statistics_bamboo_canopy_43_bw.rds")
statistics_bamboo_canopy_43_bw<-readRDS("data/data_analyses/stats_analyses/statistics_bamboo_canopy_43_bw.rds")
mean(statistics_bamboo_canopy_43_bw$distribution [,1])
sd(statistics_bamboo_canopy_43_bw$distribution [,1])
mean (statistics_bamboo_canopy_43_bw$distribution [,3]) # fraction unique bamboo
sd (statistics_bamboo_canopy_43_bw$distribution [,3])
mean (statistics_bamboo_canopy_43_bw$distribution [,4]) # fraction unique canopy
sd (statistics_bamboo_canopy_43_bw$distribution [,4])


statistics_bamboo_montane_43_bw=hypervolume_overlap_confidence(com_hypervolumes_randomized_bamboo_43_10,com_hypervolumes_randomized_montane_43_10,CI = .95, cores = 7) # 43
saveRDS(statistics_bamboo_montane_43_bw,"data/data_analyses/stats_analyses/statistics_bamboo_montane_43_bw.rds")
statistics_bamboo_montane_43_bw<-readRDS("data/data_analyses/stats_analyses/statistics_bamboo_montane_43_bw.rds")
mean(statistics_bamboo_montane_43_bw$distribution [,1])
sd(statistics_bamboo_montane_43_bw$distribution [,1])
mean (statistics_terrafirme_canopy_43_bw$distribution [,3]) # fraction unique bamboo
sd (statistics_bamboo_montane_43_bw$distribution [,3])
mean (statistics_bamboo_montane_43_bw$distribution [,4]) # fraction unique montane
sd (statistics_bamboo_montane_43_bw$distribution [,4])



statistics_bamboo_h_montane_43_bw=hypervolume_overlap_confidence(com_hypervolumes_randomized_bamboo_43_10,com_hypervolumes_randomized_hmontane_43_10,CI = .95, cores = 7) # 43
saveRDS(statistics_bamboo_h_montane_43_bw, "data/data_analyses/stats_analyses/statistics_bamboo_h_montane_43_bw.rds") 
statistics_terrafirme_bamboo_43<-readRDS("data/data_analyses/stats_analyses/statistics_bamboo_h_montane_43_bw.rds")
mean(statistics_bamboo_h_montane_43_bw$distribution [,1])
sd(statistics_bamboo_h_montane_43_bw$distribution [,1])
mean (statistics_bamboo_h_montane_43_bw$distribution [,3]) # fraction unique bamboo
sd (statistics_bamboo_h_montane_43_bw$distribution [,3])
mean (statistics_bamboo_h_montane_43_bw$distribution [,4]) # fraction unique high montane
sd (statistics_bamboo_h_montane_43_bw$distribution [,4])

statistics_canopy_montane_43_bw=hypervolume_overlap_confidence(com_hypervolumes_randomized_canopy_43_10,com_hypervolumes_randomized_montane_43_10,CI = .95, cores = 7) # 43
saveRDS(statistics_canopy_montane_43_bw,"data/data_analyses/stats_analyses/statistics_canopy_montane_43_bw.rds")
statistics_canopy_montane_43_bw<-readRDS("data/data_analyses/stats_analyses/statistics_canopy_montane_43_bw.rds")
mean(statistics_canopy_montane_43_bw$distribution [,1])
sd(statistics_canopy_montane_43_bw$distribution [,1])
mean (statistics_canopy_montane_43_bw$distribution [,3]) # fraction unique canopy
sd (statistics_canopy_montane_43_bw$distribution [,3])
mean (statistics_canopy_montane_43_bw$distribution [,4]) # fraction unique montane
sd (statistics_canopy_montane_43_bw$distribution [,4])

statistics_canopy_h_montane_43_bw=hypervolume_overlap_confidence(com_hypervolumes_randomized_canopy_43_10,com_hypervolumes_randomized_hmontane_43_10,CI = .95, cores = 7) # 43
saveRDS(statistics_canopy_h_montane_43_bw,"data/data_analyses/stats_analyses/statistics_canopy_h_montane_43_bw.rds")
statistics_canopy_h_montane_43_bw<-readRDS("data/data_analyses/stats_analyses/statistics_canopy_h_montane_43_bw.rds")
mean(statistics_canopy_h_montane_43_bw$distribution [,1])
sd(statistics_canopy_h_montane_43_bw$distribution [,1])
mean (statistics_canopy_h_montane_43_bw$distribution [,3]) # fraction unique canopy
sd (statistics_canopy_h_montane_43_bw$distribution [,3])
mean (statistics_canopy_h_montane_43_bw$distribution [,4]) # fraction unique h_montane
sd (statistics_canopy_h_montane_43_bw$distribution [,4])

statistcs_montane_hmontane_43_bw=hypervolume_overlap_confidence(com_hypervolumes_randomized_montane_43_10,com_hypervolumes_randomized_hmontane_43_10,CI = .95, cores = 7) # 43
saveRDS(statistcs_montane_hmontane_43_bw,"data/data_analyses/stats_analyses/statistcs_montane_hmontane_43_bw.rds")
statistcs_montane_hmontane_43_bw<-readRDS("data/data_analyses/stats_analyses/statistcs_montane_hmontane_43_bw.rds")
mean(statistcs_montane_hmontane_43_bw$distribution [,1])
sd(statistcs_montane_hmontane_43_bw$distribution [,1])
mean (statistcs_montane_hmontane_43_bw$distribution [,3]) # fraction unique montane
sd (statistcs_montane_hmontane_43_bw$distribution [,3])
mean (statistcs_montane_hmontane_43_bw$distribution [,4]) # fraction unique h_montane
sd (statistcs_montane_hmontane_43_bw$distribution [,4])


# Hypervolume plot  -------------------------------------------------------

install.packages("rgl")
library(rgl)

## S3 method for class 'HypervolumeList'
plot(x,
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=TRUE, show.density=TRUE,show.data=TRUE,
     names=NULL, show.legend=TRUE, limits=NULL,
     show.contour=TRUE, contour.lwd=1.5,
     contour.type='kde',
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1,
     contour.kde.level=1e-04,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     colors=rainbow(floor(length(x@HVList)*1.5),alpha=0.8),
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=0.75,cex.axis=0.75,cex.names=1.0,cex.legend=0.75,
     num.points.max.data = 1000, num.points.max.random = 2000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE,
     ...)


plot(hv,show.3d=TRUE)).

#### Plot 3D hypervolumes for mammals and birds

hypervolume_set(obs_lowlands_hypervolume_gaussian_3,obs_bamboo_hypervolume_gaussian_3,obs_montane_hypervolume_gaussian_3,obs_canopy_hypervolume_gaussian_3,obs_h_montane_hypervolume_gaussian_3, check.memory = FALSE)


set1<-hypervolume_set(obs_nonsocial_hypervolume_gaussian_3, obs_social_hypervolume_gaussian_3,check.memory = FALSE)

set2<-hypervolume_set(obs_nonsocial_hypervolume_gaussian_3, obs_social_hypervolume_gaussian_3,check.memory = FALSE)

set2<-hypervolume_set(obs_community_hypervolume_gaussian_3, obs_social_hypervolume_gaussian_3,check.memory = FALSE)




```{r}
set1@HVList$HV1 <- NULL
set1@HVList$HV2 <- NULL
set1@HVList$Union <- NULL

names <- c("Beak lenght", "Mass", "Tail lenght", "Wing Lenght","Hand-wind index","Foraging","Diet")

colnames(set1@HVList$Intersection@RandomPoints) <- names
colnames(set1@HVList$Unique_1@RandomPoints) <- names
colnames(set1@HVList$Unique_2@RandomPoints) <- names

colnames(set1@HVList$Intersection@RandomPoints)

# plot hypervolumes for mammals and birds
hypervolume::plot.HypervolumeList(set1,
                                  show.3d = TRUE, 
                                  plot.3d.axes.id = c(2,6,7), 
                                  show.random = FALSE, 
                                  show.data = TRUE, 
                                  show.centroid = FALSE,
                                  show.density = FALSE, 
                                  cex.random = 5, 
                                  colors = c("gray77", "cadetblue3","khaki3"))
```

cadetblue1




plot.HypervolumeList

#Data frame for Plotly Hyper Volume Plot
inter <- data.frame(bm_set_mi_1911vs2016@HVList$Intersection@RandomPoints)
inter <- inter %>% 
  mutate(section = "intersection")

uni1 <- data.frame(bm_set_mi_1911vs2016@HVList$Unique_1@RandomPoints) 
uni1 <- uni1 %>% 
  mutate(section = "Unique 1911")

uni2 <- data.frame(bm_set_mi_1911vs2016@HVList$Unique_2@RandomPoints)
uni2 <- uni2 %>% 
  mutate(section = "Unique 2016")

HV <- rbind(inter, uni1, uni2)
write.csv(HV,"HyperVolume.csv")

colnames(bm_set_mi_1911vs2016@HVList$Intersection@RandomPoints) <- names
colnames(bm_set_mi_1911vs2016@HVList$Unique_1@RandomPoints) <- names
colnames(bm_set_mi_1911vs2016@HVList$Unique_2@RandomPoints) <- names


library(plotly)

HV <- read.csv("HyperVolume.csv", h = T)

pal <- c("azure3","red", "cornflowerblue" )


fig <- plot_ly(HV, x = ~Comp.1, y = ~(-Comp.2), z = ~(-Comp.3), type = 'scatter3d', mode = 'markers', color = ~section, colors = pal, marker = list(size = 3,opacity = 0.6, line = list(color = 'grey66', width = 2)))

fig <- fig %>% layout(scene = list(xaxis = list(title = 'Body size'),
                                   yaxis = list(title = 'Dispersal ability'),
                                   zaxis = list(title = 'Habitat breadth')))

fig


######END OF USED FOR MANUSCRIPT







# Similar analyses for a less conservative camculation using band width 0.3
# PART 1 comparison community vs flcoks [ Observed hypervolumes]  (Controling for parameters (samples per point 2000))---------------------------------------------------
# Estimation bandwidth fixed 0.30 fixed/comparison community vs flcoks ------------------------------------------

# although later I realized that does not matter  really the samples per point it only changes teh decimals for example look at the following two 
obs_community_hypervolume_gaussian
obs_community_hypervolume_gaussian_1

### 1 #### Generate the observed hypervolumes for the community and Social
#Hypervolume community
#
set.seed(9) # Reproducible random number generator 
obs_community_hypervolume_gaussian_2<-hypervolume_gaussian(bird_traits_z_community, 
                                                           name="Community_2000_samples",
                                                           samples.per.point =2000,
                                                           weight = NULL,
                                                           sd.count= 3,
                                                           quantile.requested = 0.99,
                                                           quantile.requested.type = "probability")
#kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
kde.bandwidth= estimate_bandwidth(bird_traits_z_community, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)
saveRDS(obs_community_hypervolume_gaussian_1, "data/data_analyses/obs_community_hypervolume_gaussian_1_2000samples.rds") #SAVED
obs_community_hypervolume_gaussian_1<-readRDS("data/data_analyses/obs_community_hypervolume_gaussian_1_2000samples.rds")

#Hypervolume social
set.seed(9) #
obs_social_hypervolume_gaussian_2<-hypervolume_gaussian(bird_traits_z_social, 
                                                        name="All flocks_2000_samples",
                                                        samples.per.point =2000,
                                                        weight = NULL,
                                                        sd.count= 3,
                                                        quantile.requested = 0.99,
                                                        quantile.requested.type = "probability")
#kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
kde.bandwidth= estimate_bandwidth(bird_traits_z_social, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)
saveRDS(obs_social_hypervolume_gaussian_1, "data/data_analyses/obs_social_hypervolume_gaussian_1_2000samples.rds") #SAVED
obs_social_hypervolume_gaussian_1<-readRDS("data/data_analyses/obs_social_hypervolume_gaussian_1_2000samples.rds")


#Hypervolume nonsocial
set.seed(9)
obs_nonsocial_hypervolume_gaussian_2<-hypervolume_gaussian(bird_traits_z_nonsocial, 
                                                           name="non social_2000_samples",
                                                           samples.per.point =2000,
                                                           weight = NULL,
                                                           sd.count= 3,
                                                           quantile.requested = 0.99,
                                                           quantile.requested.type = "probability")
#kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5))
kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)
saveRDS(obs_nonsocial_hypervolume_gaussian, "data/data_analyses/obs_nonsocial_hypervolume_gaussian_1_2000samples.rds") #SAVED
obs_nonsocial_hypervolume_gaussian_1<-readRDS("data/data_analyses/obs_nonsocial_hypervolume_gaussian_1_2000samples.rds")

# flocks across the gradient fixed bandwidth 0.3 -----------------------------------------------------

# Obserevd hypervolumes controling for parameters (samples per point 2000)

# create the hypervolume  for lowlands flocks with the 2000 samples per point (see than here I a keeping the number of samples per point constant in al hypervolumes)
set.seed(9)
obs_lowlands_hypervolume_gaussian_1<-hypervolume_gaussian(bird_traits_z_und_lowlands, 
                                                          name="lowlands_2000_samples",
                                                          samples.per.point =2000,
                                                          weight = NULL,
                                                          sd.count= 3,
                                                          quantile.requested = 0.99,
                                                          quantile.requested.type = "probability",
                                                          kde.bandwidth= estimate_bandwidth(bird_traits_z_und_lowlands, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
                                                          #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)

saveRDS(obs_lowlands_hypervolume_gaussian_1, "data/data_analyses/obs_lowlands_hypervolume_gaussian_1_2000samples.rds") 
obs_lowlands_hypervolume_gaussian_1<-readRDS("data/data_analyses/obs_lowlands_hypervolume_gaussian_1_2000samples.rds")

# create the hypervolume for canopy flocks ( see than here I a keeping the number of samples per point constant )
set.seed(9)
obs_canopy_hypervolume_gaussian_1<-hypervolume_gaussian(bird_traits_z_canopy, 
                                                        name="canopy_2000_samples",
                                                        samples.per.point=2000,
                                                        weight = NULL,
                                                        sd.count= 3,
                                                        quantile.requested = 0.99,
                                                        quantile.requested.type = "probability",
                                                        kde.bandwidth= estimate_bandwidth(bird_traits_z_canopy, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
                                                        #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)

saveRDS(obs_canopy_hypervolume_gaussian_1, "data/data_analyses/obs_canopy_hypervolume_gaussian_1_2000samples.rds") 
obs_canopy_hypervolume_gaussian_1<-readRDS("data/data_analyses/obs_canopy_hypervolume_gaussian_1_2000samples.rds")

set.seed(9)
obs_bamboo_hypervolume_gaussian_1<-hypervolume_gaussian(bird_traits_z_bamboo, 
                                                        name="bamboo_2000_samples",
                                                        samples.per.point =2000,
                                                        weight = NULL,
                                                        sd.count= 3,
                                                        quantile.requested = 0.99,
                                                        quantile.requested.type = "probability",
                                                        #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5))
                                                        kde.bandwidth= estimate_bandwidth(bird_traits_z_bamboo, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)

saveRDS(obs_bamboo_hypervolume_gaussian_1, "data/data_analyses/obs_bamboo_hypervolume_gaussian_1_2000samples.rds") 
obs_bamboo_hypervolume_gaussian_1<-readRDS("data/data_analyses/obs_bamboo_hypervolume_gaussian_1_2000samples.rds")

set.seed(9)
obs_montane_hypervolume_gaussian_1<-hypervolume_gaussian(bird_traits_z_andes, 
                                                         name="montane_2000_samples",
                                                         samples.per.point =2000,
                                                         weight = NULL,
                                                         sd.count= 3,
                                                         quantile.requested = 0.99,
                                                         quantile.requested.type = "probability",
                                                         #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5))
                                                         kde.bandwidth= estimate_bandwidth(bird_traits_z_andes, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)

saveRDS(obs_montane_hypervolume_gaussian_1, "data/data_analyses/obs_montane_hypervolume_gaussian_1_2000samples.rds") 
obs_montane_hypervolume_gaussian_1<-readRDS("data/data_analyses/obs_montane_hypervolume_gaussian_1_2000samples.rds")

set.seed(9)
obs_h_montane_hypervolume_gaussian_1<-hypervolume_gaussian(bird_traits_z_h_andes, 
                                                           name="h_montane_2000_samples",
                                                           samples.per.point =2000,
                                                           weight = NULL,
                                                           sd.count= 3,
                                                           quantile.requested = 0.99,
                                                           quantile.requested.type = "probability",
                                                           #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5))
                                                           kde.bandwidth= estimate_bandwidth(bird_traits_z_h_andes, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)
saveRDS(obs_h_montane_hypervolume_gaussian_1, "data/data_analyses/obs_h_montane_hypervolume_gaussian_1_2000samples.rds") 
obs_h_montane_hypervolume_gaussian_1<-readRDS("data/data_analyses/obs_h_montane_hypervolume_gaussian_1_2000samples.rds") 

# subPart 1_ Community_hypervolume matching sample size bd 0.3----------------------

#Hypervolume with matching sample size (minimum diversity in a flocking assemblage 269)

#Commmunity
com_hypervolumes_randomized_community_269=hypervolume_resample("269_community_resample", obs_community_hypervolume_gaussian_1,
                                                               method="bootstrap",
                                                               n = 100, 
                                                               points_per_resample = "269", 
                                                               cores = 6,
                                                               verbose = TRUE, 
                                                               weight_func = NULL)

com_hypervolumes_randomized_community_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_community_resample"


#Flocking

com_hypervolumes_randomized_flocking_269=hypervolume_resample("269_flocking_resample", obs_social_hypervolume_gaussian_1,
                                                              method="bootstrap",
                                                              n = 100, 
                                                              points_per_resample = "269", 
                                                              cores = 6,
                                                              verbose = TRUE,
                                                              weight_func = NULL)
com_hypervolumes_randomized_flocking_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_flocking_resample"

#No-flocking 

com_hypervolumes_randomized_nonflocking_269=hypervolume_resample("269_nonflocking_resample", obs_nonsocial_hypervolume_gaussian_1,
                                                                 method="bootstrap",
                                                                 n = 100, 
                                                                 points_per_resample = "269", 
                                                                 cores = 6,
                                                                 verbose = TRUE, 
                                                                 weight_func = NULL)

com_hypervolumes_randomized_nonflocking_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_nonflocking_resample"


# # Flock types_Hypervolume with matching sample size bw 0.3 --------------

# Hypervolume with matching sample size (minimum diversity in a flock assemblage 43)

#Terrafirme
# make sure you keep the sample size identical to what you want to compare with 
com_hypervolumes_randomized_terrafirme_43=hypervolume_resample("43_terrafirme_resample", obs_lowlands_hypervolume_gaussian_1,
                                                               method="bootstrap",
                                                               n = 50, 
                                                               points_per_resample = "43", 
                                                               cores = 6,
                                                               verbose = TRUE, 
                                                               weight_func = NULL)

com_hypervolumes_randomized_terrafirme_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_terrafirme_resample"

#Canopy

com_hypervolumes_randomized_canopy_43=hypervolume_resample("43_canopy_resample", 
                                                           obs_canopy_hypervolume_gaussian_1, 
                                                           method="bootstrap", 
                                                           n = 50, 
                                                           points_per_resample =43,  
                                                           cores = 6,
                                                           verbose = TRUE, 
                                                           weight_func = NULL)



com_hypervolumes_randomized_canopy_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_canopy_resample"

#Bamboo

com_hypervolumes_randomized_bamboo_43=hypervolume_resample("43_bamboo_resample", 
                                                           obs_bamboo_hypervolume_gaussian_1, 
                                                           method="bootstrap", 
                                                           n = 50, 
                                                           points_per_resample ="sample_size", 
                                                           cores = 6,
                                                           verbose = TRUE, 
                                                           weight_func = NULL)

com_hypervolumes_randomized_bamboo_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_bamboo_resample"



#Montane
com_hypervolumes_randomized_montane_43=hypervolume_resample("43_montane_resample", 
                                                            obs_montane_hypervolume_gaussian_1, 
                                                            method="bootstrap", 
                                                            n = 50, 
                                                            points_per_resample=43,  # make sure you keep the sample size identical to what you want to compare with 
                                                            cores = 6,
                                                            verbose = TRUE, 
                                                            weight_func = NULL)

com_hypervolumes_randomized_montane_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_montane_resample"



#High montane


com_hypervolumes_randomized_hmontane_43=hypervolume_resample("43_hmontane_resample", 
                                                             obs_h_montane_hypervolume_gaussian_1, 
                                                             method="bootstrap", 
                                                             n = 50, 
                                                             points_per_resample =43,  # make sure you keep the sample size identical to what you want to compare with 
                                                             cores = 6,
                                                             verbose = TRUE, 
                                                             weight_func = NULL)
com_hypervolumes_randomized_hmontane_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_hmontane_resample"


# List of matched assemblages 

com_hypervolumes_randomized_nonflocking_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_nonflocking_resample"
com_hypervolumes_randomized_community_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_community_resample"
com_hypervolumes_randomized_flocking_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_flocking_resample"


# List of matched flocks  

com_hypervolumes_randomized_terrafirme_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_terrafirme_resample"

com_hypervolumes_randomized_bamboo_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_bamboo_resample"

com_hypervolumes_randomized_canopy_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_canopy_resample"

com_hypervolumes_randomized_montane_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_montane_resample"

com_hypervolumes_randomized_hmontane_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_hmontane_resample"



# # Hypervolume ( acroos the gradient) with matching sample size bw 03  --------


#Terrafirme
# make sure you keep the sample size identical to what you want to compare with 
com_hypervolumes_randomized_terrafirme_43=hypervolume_resample("43_terrafirme_resample", obs_lowlands_hypervolume_gaussian_1,
                                                               method="bootstrap",
                                                               n = 50, 
                                                               points_per_resample = "43", 
                                                               cores = 6,
                                                               verbose = TRUE, 
                                                               weight_func = NULL)

com_hypervolumes_randomized_terrafirme_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_terrafirme_resample"

#Canopy

com_hypervolumes_randomized_canopy_43=hypervolume_resample("43_canopy_resample", 
                                                           obs_canopy_hypervolume_gaussian_1, 
                                                           method="bootstrap", 
                                                           n = 50, 
                                                           points_per_resample =43,  
                                                           cores = 6,
                                                           verbose = TRUE, 
                                                           weight_func = NULL)



com_hypervolumes_randomized_canopy_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_canopy_resample"

#Bamboo

com_hypervolumes_randomized_bamboo_43=hypervolume_resample("43_bamboo_resample", 
                                                           obs_bamboo_hypervolume_gaussian_1, 
                                                           method="bootstrap", 
                                                           n = 50, 
                                                           points_per_resample ="sample_size", 
                                                           cores = 6,
                                                           verbose = TRUE, 
                                                           weight_func = NULL)

com_hypervolumes_randomized_bamboo_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_bamboo_resample"



#Montane
com_hypervolumes_randomized_montane_43=hypervolume_resample("43_montane_resample", 
                                                            obs_montane_hypervolume_gaussian_1, 
                                                            method="bootstrap", 
                                                            n = 50, 
                                                            points_per_resample=43,  # make sure you keep the sample size identical to what you want to compare with 
                                                            cores = 6,
                                                            verbose = TRUE, 
                                                            weight_func = NULL)

com_hypervolumes_randomized_montane_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_montane_resample"



#High montane


com_hypervolumes_randomized_hmontane_43=hypervolume_resample("43_hmontane_resample", 
                                                             obs_h_montane_hypervolume_gaussian_1, 
                                                             method="bootstrap", 
                                                             n = 50, 
                                                             points_per_resample =43,  # make sure you keep the sample size identical to what you want to compare with 
                                                             cores = 6,
                                                             verbose = TRUE, 
                                                             weight_func = NULL)
com_hypervolumes_randomized_hmontane_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_hmontane_resample"


# List of matched assemblages 

com_hypervolumes_randomized_nonflocking_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_nonflocking_resample"
com_hypervolumes_randomized_community_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_community_resample"
com_hypervolumes_randomized_flocking_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_flocking_resample"


# List of matched flocks  

com_hypervolumes_randomized_terrafirme_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_terrafirme_resample"

com_hypervolumes_randomized_bamboo_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_bamboo_resample"

com_hypervolumes_randomized_canopy_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_canopy_resample"

com_hypervolumes_randomized_montane_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_montane_resample"

com_hypervolumes_randomized_hmontane_43="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_hmontane_resample"


### Stats for 0.6 bw 

#Creating the 50*50 hypervolumes is too much computationally, so I will create just groups of 20*20

# List of matched flocks reduced to 20 for the overlap computation 

com_hypervolumes_randomized_nonflocking_269_20="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_nonflocking_resample_20"
com_hypervolumes_randomized_community_269_20="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_community_resample_20"
com_hypervolumes_randomized_flocking_269_20="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/269_flocking_resample_20"

#need to run these three estimation 8 hours each
statistics_community_nonflocking_269=hypervolume_overlap_confidence(com_hypervolumes_randomized_community_269_20, com_hypervolumes_randomized_nonflocking_269_20,CI = .95, cores = 7) #43
saveRDS(statistics_community_nonflocking_269, "data/data_analyses/statistics_community_nonflocking_269.rds") 
mean(statistics_community_nonflocking_269$distribution [,1])
sd(statistics_community_nonflocking_269$distribution [,1])
mean (statistics_community_nonflocking_269$distribution [,3]) # fraction unique terrafirme
sd (statistics_community_nonflocking_269$distribution [,3])
mean (statistics_community_nonflocking_269$distribution [,4]) # fraction unique canopy
sd (statistics_community_nonflocking_269$distribution [,4])


statistics_community_flocking_269=hypervolume_overlap_confidence(com_hypervolumes_randomized_community_269_20, com_hypervolumes_randomized_flocking_269_20,CI = .95, cores = 7) #43
saveRDS(statistics_community_flocking_269, "data/data_analyses/statistics_community_flocking_269.rds") 
mean(statistics_community_flocking_269$distribution [,1])
sd(statistics_community_flocking_269$distribution [,1])
mean (statistics_community_flocking_269$distribution [,3]) # fraction unique community
sd (statistics_community_flocking_269$distribution [,3])
mean (statistics_community_flocking_269$distribution [,4]) # fraction unique flocks
sd (statistics_community_flocking_269$distribution [,4])


statistics_flocking_nonflocking269=hypervolume_overlap_confidence(com_hypervolumes_randomized_flocking_269_20, com_hypervolumes_randomized_nonflocking_269_20,CI = .95, cores = 7) #43
saveRDS(statistics_flocking_nonflocking269, "data/data_analyses/statistics_flocking_nonflocking269.rds") 
mean(statistics_flocking_nonflocking269$distribution [,1])
sd(statistics_flocking_nonflocking269$distribution [,1])
mean (statistics_flocking_nonflocking269$distribution [,3]) # fraction unique flocking
sd (statistics_flocking_nonflocking269$distribution [,3])
mean (statistics_flocking_nonflocking269$distribution [,4]) # fraction unique no flocking
sd (statistics_flocking_nonflocking269$distribution [,4])

com_hypervolumes_randomized_terrafirme_43_10="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_terrafirme_resample_10"

com_hypervolumes_randomized_bamboo_43_10="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_bamboo_resample_10"

com_hypervolumes_randomized_canopy_43_10="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_canopy_resample_10"

com_hypervolumes_randomized_montane_43_10="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_montane_resample_10"

com_hypervolumes_randomized_hmontane_43_10="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/43_hmontane_resample_10"


statistics_terrafirme_canopy_43=hypervolume_overlap_confidence(com_hypervolumes_randomized_terrafirme_43_10, com_hypervolumes_randomized_canopy_43_10,CI = .95, cores = 7) #43
saveRDS(statistics_terrafirme_canopy_43, "data/data_analyses/statistics_terrafirme_canopy_43.rds") 
statistics_terrafirme_canopy_43<-readRDS("data/data_analyses/statistics_terrafirme_canopy_43.rds")
mean(statistics_terrafirme_canopy_43$distribution [,1])
sd(statistics_terrafirme_canopy_43$distribution [,1])
mean (statistics_terrafirme_canopy_43$distribution [,3]) # fraction unique terrafirme
sd (statistics_terrafirme_canopy_43$distribution [,3])
mean (statistics_terrafirme_canopy_43$distribution [,4]) # fraction unique canopy
sd (statistics_terrafirme_canopy_43$distribution [,4])



statistics_terrafirme_bamboo_43=hypervolume_overlap_confidence(com_hypervolumes_randomized_terrafirme_43_10,com_hypervolumes_randomized_bamboo_43_10,CI = .95, cores = 7) # 43
saveRDS(statistics_terrafirme_bamboo_43, "data/data_analyses/statistics_terrafirme_bamboo_43.rds") 
statistics_terrafirme_bamboo_43<-readRDS("data/data_analyses/statistics_terrafirme_bamboo_43.rds")
mean(statistics_terrafirme_bamboo_43$distribution [,1])
sd(statistics_terrafirme_bamboo_43$distribution [,1])
mean (statistics_terrafirme_bamboo_43$distribution [,3]) # fraction unique terrafirme
sd (statistics_terrafirme_bamboo_43$distribution [,3])
mean (statistics_terrafirme_bamboo_43$distribution [,4]) # fraction unique bamboo
sd (statistics_terrafirme_bamboo_43$distribution [,4])


statistics_terrafirme_montane_43=hypervolume_overlap_confidence(com_hypervolumes_randomized_terrafirme_43_10,com_hypervolumes_randomized_montane_43_10,CI = .95, cores = 7) # 43
saveRDS(statistics_terrafirme_montane_43,"data/data_analyses/statistics_terrafirme_montane_43.rds")
statistics_bamboo_montane_43<-readRDS("data/data_analyses/statistics_terrafirme_montane_43.rds")
mean(statistics_terrafirme_montane_43$distribution [,1])
sd(statistics_terrafirme_montane_43$distribution [,1])
mean (statistics_terrafirme_montane_43$distribution [,3]) # fraction unique terrafirme
sd (statistics_terrafirme_montane_43$distribution [,3])
mean (statistics_terrafirme_montane_43$distribution [,4]) # fraction unique montane
sd (statistics_terrafirme_montane_43$distribution [,4])


statistics_terrafirme_h_montane_43=hypervolume_overlap_confidence(com_hypervolumes_randomized_terrafirme_43_10,com_hypervolumes_randomized_hmontane_43_10,CI = .95, cores = 7) # 43
saveRDS(statistics_terrafirme_h_montane_43,"data/data_analyses/statistics_terrafirme_h_montane_43.rds")
statistics_terrafirme_h_montane_43<-readRDS("data/data_analyses/statistics_terrafirme_h_montane_43.rds")
mean(statistics_terrafirme_h_montane_43$distribution [,1])
sd(statistics_terrafirme_h_montane_43$distribution [,1])
mean (statistics_terrafirme_h_montane_43$distribution [,3]) # fraction unique terrafirme
sd (statistics_terrafirme_h_montane_43$distribution [,3])
mean (statistics_terrafirme_h_montane_43$distribution [,4]) # fraction unique h_montane
sd (statistics_terrafirme_h_montane_43$distribution [,4])


statistics_bamboo_canopy_43=hypervolume_overlap_confidence(com_hypervolumes_randomized_bamboo_43_10, com_hypervolumes_randomized_canopy_43_10, CI = .95, cores = 7) # 43
saveRDS(statistics_bamboo_canopy_43,"data/data_analyses/statistics_bamboo_canopy_43.rds")
statistics_bamboo_canopy_43<-readRDS("data/data_analyses/statistics_bamboo_canopy_43.rds")
mean(statistics_bamboo_canopy_43$distribution [,1])
sd(statistics_bamboo_canopy_43$distribution [,1])
mean (statistics_bamboo_canopy_43$distribution [,3]) # fraction unique bamboo
sd (statistics_bamboo_canopy_43$distribution [,3])
mean (statistics_bamboo_canopy_43$distribution [,4]) # fraction unique canopy
sd (statistics_bamboo_canopy_43$distribution [,4])


statistics_bamboo_montane_43=hypervolume_overlap_confidence(com_hypervolumes_randomized_bamboo_43_10,com_hypervolumes_randomized_montane_43_10,CI = .95, cores = 7) # 43
saveRDS(statistics_bamboo_montane_43,"data/data_analyses/statistics_bamboo_montane_43.rds")
statistics_bamboo_montane_43<-readRDS("data/data_analyses/statistics_bamboo_montane_43.rds")
mean(statistics_bamboo_montane_43$distribution [,1])
sd(statistics_bamboo_montane_43$distribution [,1])
mean (statistics_terrafirme_canopy_43$distribution [,3]) # fraction unique bamboo
sd (statistics_bamboo_montane_43$distribution [,3])
mean (statistics_bamboo_montane_43$distribution [,4]) # fraction unique montane
sd (statistics_bamboo_montane_43$distribution [,4])



statistics_bamboo_h_montane_43=hypervolume_overlap_confidence(com_hypervolumes_randomized_bamboo_43_10,com_hypervolumes_randomized_hmontane_43_10,CI = .95, cores = 7) # 43
saveRDS(statistics_terrafirme_bamboo_43, "data/data_analyses/statistics_bamboo_h_montane_43.rds") 
statistics_terrafirme_bamboo_43<-readRDS("data/data_analyses/statistics_bamboo_h_montane_43.rds")
mean(statistics_bamboo_h_montane_43$distribution [,1])
sd(statistics_bamboo_h_montane_43$distribution [,1])
mean (statistics_bamboo_h_montane_43$distribution [,3]) # fraction unique bamboo
sd (statistics_bamboo_h_montane_43$distribution [,3])
mean (statistics_bamboo_h_montane_43$distribution [,4]) # fraction unique high montane
sd (statistics_bamboo_h_montane_43$distribution [,4])

statistics_canopy_montane_43=hypervolume_overlap_confidence(com_hypervolumes_randomized_canopy_43_10,com_hypervolumes_randomized_montane_43_10,CI = .95, cores = 7) # 43
saveRDS(statistics_canopy_montane,"data/data_analyses/statistics_canopy_montane_43.rds")
statistics_canopy_montane_43<-readRDS("data/data_analyses/statistics_canopy_montane_43.rds")
mean(statistics_canopy_montane_43$distribution [,1])
sd(statistics_canopy_montane_43$distribution [,1])
mean (statistics_canopy_montane_43$distribution [,3]) # fraction unique canopy
sd (statistics_canopy_montane_43$distribution [,3])
mean (statistics_canopy_montane_43$distribution [,4]) # fraction unique montane
sd (statistics_canopy_montane_43$distribution [,4])

statistics_canopy_h_montane_43=hypervolume_overlap_confidence(com_hypervolumes_randomized_canopy_43_10,com_hypervolumes_randomized_hmontane_43_10,CI = .95, cores = 7) # 43
saveRDS(statistics_canopy_h_montane,"data/data_analyses/statistics_canopy_h_montane_43.rds")
statistics_canopy_h_montane_43<-readRDS("data/data_analyses/statistics_canopy_h_montane_43.rds")
mean(statistics_canopy_h_montane_43$distribution [,1])
sd(statistics_canopy_h_montane_43$distribution [,1])
mean (statistics_canopy_h_montane_43$distribution [,3]) # fraction unique canopy
sd (statistics_canopy_h_montane_43$distribution [,3])
mean (statistics_canopy_h_montane_43$distribution [,4]) # fraction unique h_montane
sd (statistics_canopy_h_montane_43$distribution [,4])

statistcs_montane_hmontane_43=hypervolume_overlap_confidence(com_hypervolumes_randomized_montane_43_10,com_hypervolumes_randomized_hmontane_43_10,CI = .95, cores = 7) # 43
saveRDS(statistcs_montane_hmontane_43,"data/data_analyses/statistcs_montane_hmontane_43.rds")
statistcs_montane_hmontane_43<-readRDS("data/data_analyses/statistcs_montane_hmontane_43.rds")
mean(statistcs_montane_hmontane_43$distribution [,1])
sd(statistcs_montane_hmontane_43$distribution [,1])
mean (statistcs_montane_hmontane_43$distribution [,3]) # fraction unique montane
sd (statistcs_montane_hmontane_43$distribution [,3])
mean (statistcs_montane_hmontane_43$distribution [,4]) # fraction unique h_montane
sd (statistcs_montane_hmontane_43$distribution [,4])

