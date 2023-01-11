#######################################################################################
### Manuscript-Flocks in a community context                                         ###
###                                                                                 ###
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

###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
## Section 2 Data loading------------------------------------------------------------------
###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_

# Traits for Manu bird species raw; File incorporating traits form Manu, Avonet, diet and foraging.

#### Data formatted for PCA analyses

#traits_birds_binomial<-read.csv("data/data_analyses/4.TIDY_all_traits_compiled_birds_manu_binomial.csv") 



## Section 3 Data preparation --------------------------------------------------------
## Section 3.1 Needed to generate the PCA data #### Diet reduction PCoA  ####
###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_

#### Calculate synthetic diet trait (principal component from a PCoA of 10 diet categories)

diet_all <- traits_birds_continous %>% 
  # select diet data
  dplyr::select(species_checked, contains("diet")) %>% 
  # add species names to rownames (needed for gowdis function)
  tibble::column_to_rownames(var = "species_checked") %>% 
  as.data.frame()

# calculate species x species gower distance matrix based on traits
diet_gd <- FD::gowdis(diet_all)
str(diet_gd)
# to visualize the matrix as.matrix(diet_gd)

# perform principal coordinates analysis (PCoA)
diet_pco <- ade4::dudi.pco(diet_gd, scannf = FALSE)
pc_diet <- diet_pco$tab
summary(diet_pco)

#alternatively do the principal coordinates analysis (PCoA) using 'ape'package pcoa function
#matrix<-ape::pcoa(diet_gd)
#diet_pco_alt<-(matrix$vectors)
#summary(diet_pco_alt)
# R/ Projected inertia Axis1 = 42.234
# R/ Projected inertia Axis2 =  23.154
# principle component axes
pcomp_diet <- as.data.frame(diet_pco$tab[,1:4]) %>% 
  tibble::rownames_to_column(var = "species_checked")

# diet category projection
n <- nrow(diet_all)
points_stand <- scale(diet_pco$tab[,1:2])
S <- cov(diet_all, points_stand) # loadings of the variables in eacg component
U <- S %*% diag((diet_pco$eig[1:2]/(n-1))^(-0.5))
colnames(U) <- colnames(diet_pco$tab[,1:2])

# diet categoires (see Wilman et al., 2014)
U <- as.data.frame(U) %>% 
  mutate(trait = c("Inv", "Vend", "Vect", "Vfish", "Vunk", "Scav", "Fruit", "Nect", "Seed", "Planto"))
# Inv - Invertebrates
# Vend - Vertebrate endotherms 
# Vect - Vertebrate ectotherms
# Vfish - Fish # Vunk - Vertebrate unknown or general
# Scav - Scavenge 
# Fruit - Fruit
# Nect - Nectar 
# Seed - Seeds
# Planto - Other plant material
summary(pcomp_diet)
# scale diet category arrows
mult <- min(
  (max(pcomp_diet$A2) - min(pcomp_diet$A2)/(max(U$A2)-min(U$A2))),
  (max(pcomp_diet$A1) - min(pcomp_diet$A1)/(max(U$A1)-min(U$A1)))
)

U <- U %>% 
  mutate(v1 = 0.0003 * mult * A1) %>% 
  mutate(v2 = 0.0003 * mult * A2)

# plot diet PCoA
pcoa_diet<- ggplot(pcomp_diet, aes(A1, A2)) +
  # set up plot
  geom_hline(yintercept = 0, size = 0.2, lty = 2, colour = "grey") + 
  geom_vline(xintercept = 0, size = 0.2, lty = 2, colour = "grey") +
  # add origin lines
  #geom_text(alpha = 0.4, size = 1, aes(label = binomial))
  geom_point( size = 4) +
  # add species
  coord_equal() +
  geom_segment(data = U, aes(x = 0, y = 0, xend = v1, yend = v2), colour = "black",size = 3) +
  # add arrows
  geom_text(data = U, aes(x = v1, y = v2, label = trait), size = 8, colour = "black",
            nudge_y = c(rep(0, 6), 0.005, 0.005, 0.0005, -0.004), 
            nudge_x = c(0.005, rep(0, 7), -0.009, 0)) +
  # add arrow labels
  labs(x = "PC1 (42.7%)", y = "PC2 (23.4%)") +
  ggtitle("Diet")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "grey", size=4),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 30))

pcoa_diet

ggsave("figures/fig_manuscript/Sumpl.PCoA_Diet_plot.pdf",   width=15, height=15, units="in")




###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
## Section 3.2 Needed to generate the PCA data #### Foraging data reduction PCoA  ####
###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_

foraging_all <- traits_birds_continous %>% 
  # select diet data
  dplyr::select(species_checked, contains("ForStrat")) %>% 
  # add species names to rownames (needed for gowdis function)
  tibble::column_to_rownames(var = "species_checked") %>% 
  as.data.frame()

# calculate species x species gower distance matrix based on traits
foraging_gd <- FD::gowdis(foraging_all)

str(foraging_gd)

# perform principal coordinates analysis (PCoA)
foraging_pco <- ade4::dudi.pco(foraging_gd, scannf = FALSE)
pc_foraging<- foraging_pco$tab

summary(foraging_pco)

# Projected inertia Axis1 = 40.9
# Projected inertia Axis2 =  25.44

# principle component axes
pcomp_foraging <- as.data.frame(foraging_pco$tab[,1:4]) %>% 
  tibble::rownames_to_column(var = "species_checked")

#% add this line to rename variables, but will conflict with plot output so demoved to create the plot only 
#rename(A1_foraging=A1,A2_foraging=A2,A3_foraging=A3,A4_foraging=A4)

str(pcomp_foraging)

# diet category projection
n <- nrow(foraging_all)
points_stand <- scale(foraging_pco$tab[,1:2])
S <- cov(foraging_all, points_stand)
U <- S %*% diag((foraging_pco$eig[1:2]/(n-1))^(-0.5))
colnames(U) <- colnames(foraging_pco$tab[,1:2])

# diet categoires (see Wilman et al., 2014)
U <- as.data.frame(U) %>% mutate(trait = c("below_surf", "around_surf", "ground", "understory", "mid_high", "canopy", "aerial"))


# scale diet category arrows
mult <- min(
  (max(pcomp_foraging$A2) - min(pcomp_foraging$A2)/(max(U$A2)-min(U$A2))),
  (max(pcomp_foraging$A1) - min(pcomp_foraging$A1)/(max(U$A1)-min(U$A1)))
)

U <- U %>% 
  mutate(v1 = 0.0003 * mult * A1) %>% 
  mutate(v2 = 0.0003 * mult * A2)

# plot diet PCoA

pcoa_foraging <- ggplot(pcomp_foraging, aes(A1, A2)) +
  # set up plot
  geom_hline(yintercept = 0, size = 0.2, lty = 2, colour = "grey") + 
  geom_vline(xintercept = 0, size = 0.2, lty = 2, colour = "grey") +
  # add origin lines
  #geom_text(alpha = 0.4, size = 1, aes(label = binomial))
  geom_point(size=4) +
  # add species
  coord_equal() +
  geom_segment(data = U, aes(x = 0, y = 0, xend = v1, yend = v2), size=3) +
  # add arrows
  geom_text(data = U, aes(x = v1, y = v2, label = trait), size = 8, colour = "black",
            nudge_y = c(rep(0, 6), 0.005, 0.005, 0.0005, -0.004), 
            nudge_x = c(0.005, rep(0, 7), -0.009, 0)) +
  # add arrow labels
  labs(x = "PC1 (41.6%)", y = "PC2 (26.2%)") +
  ggtitle("Foraging")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "grey", size=4),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 30))


pcoa_foraging

ggsave("figures/fig_manuscript/Supl.PCoA_Foraging_plot.pdf",
       width=15, height=15, units="in")


###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_
## Section 3.3 Adding synthetic data to traits -----------------------------
###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_###_

# join diet "synthetic" data with the traits database 
# we used teh column species checked because it has the "old taxonomy" 
# using the checked column allow us to join our data set with the foraging and diet data whic both have different taxonomies
str(traits_birds_continous)
str(pcomp_diet)
# adding diet 
traits_birds_continous_diet<-full_join(traits_birds_continous, pcomp_diet, by="species_checked")

#adding foraging
traits_birds_continous_diet_foraging<-full_join(traits_birds_continous_diet,pcomp_foraging, by="species_checked")

View(traits_birds_continous_diet_foraging)



#### ----Testing for the difference in individual traits----------------------------C-----------------------------

###_###_###_###_###_
#df_traits<-df_traits%>% select-c(X.1,X)

##  Wilcoxon rank sum test with continuity correction (non-parametric)

View(df_traits)
df_traits<-na.omit(df_traits[,1:35]) # Make sure there are not NAs

##  Wilcoxon rank sum test with continuity correction
#if p>0.05 we fail to reject (accept) the null H of equal means
#if p<0.05 we conclude that there are difference between social and son social species


# a) Mass * sig difference


df_traits  %>%                                        # Specify data frame
  group_by(sociality) %>%                         # Specify group indicator
  summarise_at(vars(tail_tidy),              # Specify column
               list(name = mean, st.dev=sd))  

test_mass <- wilcox.test(df_traits$mass_tidy ~ df_traits$sociality)
test_mass

## data:  df_traits$mass_tidy by df_traits$sociality
## W = 64600, p-value = 9.433e-05
## alternative hypothesis: true location shift is not equal to 0



# b)Wing   * sig differenc
test_wing <- wilcox.test(df_traits$wing_tidy ~ df_traits$sociality)
test_wing
##
## W = 65710, p-value = 1.358e-05
##

#c)Hand-wing index
test_index <- wilcox.test(df_traits$wing_index_tidy ~ df_traits$sociality)
test_index

## W = 65785, p-value = 1.183e-05

# d)Beak * sig difference

test_beak <- wilcox.test(df_traits$bill_tidy ~ df_traits$sociality)
test_beak

## W = 70614, p-value = 2.615e-10
## alternative hypothesis: true location shift is not equal to 0
##

#e)Tail
test_tail <- wilcox.test(df_traits$tail_tidy ~ df_traits$sociality)
test_tail
##W = 57302, p-value = 0.3299

#d)Tarsus
test_tarsus <- wilcox.test(df_traits$tarsus_tidy ~ df_traits$sociality)
test_tarsus

#W = 59924, p-value = 0.04268

#The p-value is 0.021. Therefore, at the 5% significance level, we reject the null hypothesis 
#and we conclude that grades are significantly different between girls and boys.


## 
##  Wilcoxon rank sum test with continuity correction
## 
## data:  dat$Grade by dat$Sex
## W = 31.5, p-value = 0.02056
## alternative hypothesis: true location shift is not equal to 0



#### Section 4 PCA--------------------------------------------------------------
# Note: Manu species list 681; Traits list 677 [ We exclude 4 species for which we have not diet or foraging data available]
## The following species do not have diet info and have NA, we excluded the NAs
# no diet Lepicolaptes fatimaliae
# no diet data for Psittacara leucophthalmus##
# no diet data for  Pyrrhura roseifrons
# we ended not using teh foraging data from Pigot
# EXCLUDE Eucrepornis callinota because seems to be an id error

# Validate that the data is align with taxonomy
#Make sure taxonomy is the same that in the Manu data base
#df_traits_list<-df_traits %>% distinct(species)
#manu_detections<-read.csv("data/data_analyses/1.Manu_bird_species_detections_tidy_taxonomy_18052022.csv")
#manu_detections_list<-manu_detections %>% distinct( species_taxonomy_SACC_2021)
#df_traits_taxonomy_check<-right_join(df_traits, manu_detections_list, by=c("species"="species_taxonomy_SACC_2021"))
#View(df_traits_taxonomy_check)
###

# Data traits for PCA
# the final file with traits of Manu including diet and foraging from PCoA
###_###_###_###_###_
###_###_###_###_###_
#df_traits<-df_traits%>% select-c(X.1,X)

#View(df_traits)
df_traits<-na.omit(df_traits[,1:35]) # Make sure there are not NAs

#Warning!!!!! make sure you use the column species and not "species checked". Species is the 2021 taxonomy

# ## Section 4.1: PCA community## --------------------------------------------
###_###_###_###_###_###_###_###_###_###_###
# Select data
df_traits_selected<-df_traits %>% 
  select(sociality,species,mass_tidy,bill_tidy,wing_tidy,tail_tidy,wing_index_tidy, A1.x,A2.x,A1.y,A2.y) %>% 
rename(A1=A1.x,A2=A2.x,A1_foraging=A1.y,A2_foraging=A2.y)

#hist(df_traits_selected$mass_tidy, breaks = 15) # normality 
#shapiro.test(df_traits_selected$bill_tidy)# iF P<0.05 data deviate from normality, all deviate from normality

#Reformat the structure
bird_traits1<-df_traits_selected %>% mutate(log_weight=log(mass_tidy),
                                            log_tail=log(tail_tidy),
                                            log_wing=log(wing_tidy),
                                            log_beak=log(bill_tidy),
                                            log_w_index=log(wing_index_tidy)) 

#View(bird_traits1)

bird_traits2<-bird_traits1 %>% select(!c(sociality,mass_tidy,A2_foraging,A2,tail_tidy, wing_tidy,bill_tidy,wing_index_tidy ))
names(bird_traits2)
bird_traits3 <- bird_traits2[order(bird_traits2$species),]
bird_traits4 <- data.frame(bird_traits3[,-1], row.names=bird_traits3[,1])
bird_traits_df<-bird_traits4 %>% rename(Diet=A1,
                                      Foraging=A1_foraging,
                                      BeakLength=log_beak,
                                      WingLength=log_wing,
                                      TailLength=log_tail,
                                      HandWingIndex=log_w_index,
                                      BodyMass=log_weight) 


str(bird_traits_df) #677 species
## Single PCA with all traits and species for the community
## Run PCA on the non-scaled traits and let the package scale the traits, cor=TRUE.
## Alternatively you can manually scaled.

pca_birds <- princomp(bird_traits_df, cor = TRUE, scores = TRUE) # cor-TRUE this z-scale the traits 
summary(pca_birds)

#Scores
pca_scores_total<- as.data.frame(pca_birds$scores) %>% 
  tibble::rownames_to_column("species") %>% 
  left_join(df_traits_selected, by= "species")
summary(pca_scores_total)

#Loadings
pca_loadings<-pca_birds$loadings # loading tells you the contribution of each variable to the component
#or
pca_loadings_total <- as.data.frame(unclass(pca_birds$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()
###
#fviz_eig(pca_birds) # Visualize eigenvalues (scree plot). #Show the percentage of variances explained by each principal component.

# scalar to adjust arrow size
scale <-5

pca_loadings_total_scalar<- pca_loadings_total %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")),list(scalar=~.*scale))  
  # posh names
  #mutate(trait = c("Beak length", "Wing length", "Tail length", "Hand-wing index","Diet","Foraging","Body mass"))
# posh names 
#mutate(trait = c("Beak length", "Wing length", "Tail length", "Wing index", "Diet","Diet2","Body mass"))

# Eigenvalues
eig.val <- get_eigenvalue(pca_birds)
eig.val
res.var <- get_pca_var(pca_birds)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(pca_birds)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 
# Results for Variables using factoExtra
#Read some theory for interpretation here
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# ## Section 4.2: PCA for community vs social species (in a community context##) --------------------------------------------
###_###_###_###_###_###_###_###_###_###_###
# This means we will be using the PCA for the community and just plotting on top all the flocking species, using the same P.Components than the community
# Social species in a community context
#Create new database for the whole community and  to compare social and non social species in teh same plot

pca_birds_community<-pca_scores_total  %>%
  select(species,Comp.1, Comp.2, Comp.3,Comp.4,Comp.5,Comp.6,Comp.7,sociality) 

pca_birds_social<-pca_birds_community %>% filter(sociality=="yes")
pca_birds_nonsocial<-pca_birds_community %>% filter(sociality=="no")


#View(pca_birds_community)
#View(pca_birds_social)

# ## Section 4.3: PCA for community vs each flock type (in a community context## --------------------------------------------
# Creates multidimensional plots for each flock type using the whole community PCA multidimensional space
# on the flock types type= 1 underestory 2 bambo 3 canopy 4 andean 5 high andean
unique(pca_birds_community$species)
##Lowlands understory
lowlands_understory<-flocks %>% filter(flock_type_ref=="1")
lowlands_understory<-unique(lowlands_understory$species_taxonomy_SACC_2021)
lowlands_understory<-as.data.frame(lowlands_understory)
species_lowlands<-rename(lowlands_understory,c(species="lowlands_understory"))
PCA_und_lowlands<-inner_join(pca_birds_community,species_lowlands, by="species")

##Lowlands canopy
l_canopy<-flocks %>% filter(flock_type_ref=="3")
l_canopy<-unique(l_canopy$species_taxonomy_SACC_2021)
l_canopy<-as.data.frame(l_canopy)
species_l_canopy<-rename(l_canopy,c(species="l_canopy"))
PCA_l_canopy<-inner_join(pca_birds_community,species_l_canopy, by="species")

##Bamboo
bamboo<-flocks %>% filter(flock_type_ref=="2")
bamboo<-unique(bamboo$species_taxonomy_SACC_2021)
bamboo<-as.data.frame(bamboo)
species_bamboo<-rename(bamboo,c(species="bamboo"))
PCA_bamboo<-inner_join(pca_birds_community,species_bamboo, by="species")

##Andes
andes<-flocks %>% filter(flock_type_ref=="4")
andes<-unique(andes$species_taxonomy_SACC_2021)
andes<-as.data.frame(andes)
species_andes<-rename(andes,c(species="andes"))
PCA_andes<-inner_join(pca_birds_community,species_andes, by="species")

##High Andes
high_andes<-flocks %>% filter(flock_type_ref=="5")
list_flocks_high_andes<-unique(high_andes$species_taxonomy_SACC_2021)
list_flocks_high_andes<-as.data.frame(list_flocks_high_andes)
species_h_andes<-rename(list_flocks_high_andes,c(species="list_flocks_high_andes"))
PCA_h_andes<-inner_join(pca_birds_community,species_h_andes, by="species")

# ## Section 4.4: Kernel density estimation for the community and each flock (n a community context)## --------------------------------------------
###__###__######__###__######__###__###
#### kernel density estimation for the community and each flock type
###_ ###_######__###__######__###__######__###__###

## all community
# kernel density estimation for each period
k_community<- pca_scores_total%>% 
  # extract first two principal components
  dplyr::select(.,species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "species")

hpi_community<- Hpi(x = k_community)

k_est_community <- kde(x = k_community, H = hpi_community, compute.cont = TRUE)  

den_community<- list( k_est_community$eval.points[[1]],  k_est_community$eval.points[[2]],
                      k_est_community$estimate)
names(den_community) <- c("x", "y", "z")
dimnames(den_community$z) <- list(den_community$x, den_community$y)
dcc_community <- reshape2::melt(den_community$z)

#Kernel function 
k_fun <- function(df, prob) {
  dx <- diff(df$x[1:2])
  dy <- diff(df$y[1:2])
  sz <- sort(df$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}

# 0.5 probability kernel
k_50_community <- k_fun(df = den_community, prob = 0.50) # 0.5 probability kernel
k_95_community<- k_fun(df = den_community, prob = 0.90)# 0.95 probability kernel
k_99_community <- k_fun(df = den_community, prob = 0.99)#0.99 probability kernel

#Social birds
###_###_###
k_social<- pca_birds_social%>% 
  # extract first two principal components
  dplyr::select(.,species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "species")

hpi_social<- Hpi(x = k_social)
k_est_social <- kde(x = k_social, H = hpi_social, compute.cont = TRUE)  
den_social<- list( k_est_social$eval.points[[1]],  k_est_social$eval.points[[2]],
                   k_est_social$estimate)
names(den_social) <- c("x", "y", "z")
dimnames(den_social$z) <- list(den_social$x, den_social$y)
dcc_social <- reshape2::melt(den_social$z)

# 0.5 probability kernel
k_50_social <- k_fun(df = den_social, prob = 0.50) # 0.5 probability kernel
k_95_social<- k_fun(df = den_social, prob = 0.90)# 0.95 probability kernel
k_99_social<- k_fun(df = den_social, prob = 0.99)#0.99 probability kernel


#non_Social birds
###_###_###
k_nonsocial<- pca_birds_nonsocial%>% 
  # extract first two principal components
  dplyr::select(.,species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "species")

hpi_nonsocial<- Hpi(x = k_nonsocial)
k_est_nonsocial <- kde(x = k_nonsocial, H = hpi_nonsocial, compute.cont = TRUE)  
den_nonsocial<- list( k_est_nonsocial$eval.points[[1]],  k_est_nonsocial$eval.points[[2]],
                   k_est_nonsocial$estimate)
names(den_nonsocial) <- c("x", "y", "z")
dimnames(den_nonsocial$z) <- list(den_nonsocial$x, den_nonsocial$y)
dcc_nonsocial <- reshape2::melt(den_nonsocial$z)

# 0.5 probability kernel
k_50_nonsocial <- k_fun(df = den_nonsocial, prob = 0.50) # 0.5 probability kernel
k_95_nonsocial<- k_fun(df = den_nonsocial, prob = 0.90)# 0.95 probability kernel
k_99_nonsocial<- k_fun(df = den_nonsocial, prob = 0.99)#0.99 probability kernel

####_###
#lowlands
###_###_###
k_lowlands<-PCA_und_lowlands%>% 
  # extract first two principal components
  dplyr::select(.,species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "species")

hpi_lowlands<- Hpi(x = k_lowlands)
k_est_lowlands<- kde(x = k_lowlands, H = hpi_lowlands, compute.cont = TRUE)  
den_lowlands<- list( k_est_lowlands$eval.points[[1]],  k_est_lowlands$eval.points[[2]],
                     k_est_lowlands$estimate)
names(den_lowlands) <- c("x", "y", "z")
dimnames(den_lowlands$z) <- list(den_lowlands$x, den_lowlands$y)
dcc_lowlands<- reshape2::melt(den_lowlands$z)

# 0.5 probability kernel
k_50_lowlands<- k_fun(df = den_lowlands, prob = 0.50) # 0.5 probability kernel
k_95_lowlands<- k_fun(df = den_lowlands, prob = 0.90)# 0.95 probability kernel
k_99_lowlands<- k_fun(df = den_lowlands, prob = 0.99)#0.99 probability kernel

####_###
#canopy
###_###_###
k_canopy<-PCA_l_canopy%>% 
  # extract first two principal components
  dplyr::select(.,species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "species")

hpi_canopy<- Hpi(x = k_canopy)
k_est_canopy<- kde(x = k_canopy, H = hpi_canopy, compute.cont = TRUE)  
den_canopy<- list( k_est_canopy$eval.points[[1]],  k_est_canopy$eval.points[[2]],
                   k_est_canopy$estimate)
names(den_canopy) <- c("x", "y", "z")
dimnames(den_canopy$z) <- list(den_canopy$x, den_canopy$y)
dcc_canopy<- reshape2::melt(den_canopy$z)

# 0.5 probability kernel
k_50_canopy<- k_fun(df = den_canopy, prob = 0.50) # 0.5 probability kernel
k_95_canopy<- k_fun(df = den_canopy, prob = 0.90)# 0.95 probability kernel
k_99_canopy<- k_fun(df = den_canopy, prob = 0.99)#0.99 probability kernel

####_###
#bamboo
###_###_###
k_bamboo<-PCA_bamboo%>% 
  # extract first two principal components
  dplyr::select(.,species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "species")

hpi_bamboo<- Hpi(x = k_bamboo)
k_est_bamboo<- kde(x = k_bamboo, H = hpi_bamboo, compute.cont = TRUE)  
den_bamboo<- list( k_est_bamboo$eval.points[[1]],  k_est_bamboo$eval.points[[2]],
                   k_est_bamboo$estimate)
names(den_bamboo) <- c("x", "y", "z")
dimnames(den_bamboo$z) <- list(den_bamboo$x, den_bamboo$y)
dcc_bamboo<- reshape2::melt(den_bamboo$z)

# 0.5 probability kernel
k_50_bamboo<- k_fun(df = den_bamboo, prob = 0.50) # 0.5 probability kernel
k_95_bamboo<- k_fun(df = den_bamboo, prob = 0.90)# 0.95 probability kernel
k_99_bamboo<- k_fun(df = den_bamboo, prob = 0.99)#0.99 probability kernel

####_###
#ANDES
###_###_###
k_andes<-PCA_andes%>% 
  # extract first two principal components
  dplyr::select(.,species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "species")

hpi_andes<- Hpi(x = k_andes)
k_est_andes<- kde(x = k_andes, H = hpi_andes, compute.cont = TRUE)  
den_andes<- list( k_est_andes$eval.points[[1]],  k_est_andes$eval.points[[2]],
                  k_est_andes$estimate)
names(den_andes) <- c("x", "y", "z")
dimnames(den_andes$z) <- list(den_andes$x, den_andes$y)
dcc_andes<- reshape2::melt(den_andes$z)

# 0.5 probability kernel
k_50_andes<- k_fun(df = den_andes, prob = 0.50) # 0.5 probability kernel
k_95_andes<- k_fun(df = den_andes, prob = 0.90)# 0.95 probability kernel
k_99_andes<- k_fun(df = den_andes, prob = 0.99)#0.99 probability kernel

####_###
#H_ANDES
###_###_###
k_H_andes<-PCA_h_andes%>% 
  # extract first two principal components
  dplyr::select(.,species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "species")

hpi_H_andes<- Hpi(x = k_H_andes)
k_est_H_andes<- kde(x = k_H_andes, H = hpi_H_andes, compute.cont = TRUE)  
den_H_andes<- list( k_est_H_andes$eval.points[[1]],  k_est_H_andes$eval.points[[2]],
                    k_est_H_andes$estimate)
names(den_H_andes) <- c("x", "y", "z")
dimnames(den_H_andes$z) <- list(den_H_andes$x, den_H_andes$y)
dcc_H_andes<- reshape2::melt(den_H_andes$z)

# 0.5 probability kernel
k_50_H_andes<- k_fun(df = den_H_andes, prob = 0.50) # 0.5 probability kernel
k_95_H_andes<- k_fun(df = den_H_andes, prob = 0.90)# 0.95 probability kernel
k_99_H_andes<- k_fun(df = den_H_andes, prob = 0.99)#0.99 probability kernel

###_###_###_###_###_###_###_###
# ## Section 5 PCA for flocks independent of community ----------------------------------
###_###_###_###_###_###_###_###_
#Generating INDEPENDENT PCA, NOTE that this do not include all the community species and dataset are restricted to the subject of analyses
# I repited the steps that we did for the PCA, but for social species only
# this is the final file with traits of Manu including diet and foraging from PCoA

df_traits<-read.csv("data/data_analyses/4.df_traits_manu_birds.csv")
#%>% select(-c(X.1, X))

df_traits_selected<-df_traits %>% 
  select(sociality,species,mass_tidy,bill_tidy,wing_tidy,tail_tidy,wing_index_tidy, A1.x,A2.x,A1.y,A2.y) %>% 
  rename(A1=A1.x,A2=A2.x,A1_foraging=A1.y,A2_foraging=A2.y ) 

#Reformat the structure
bird_traits1<-df_traits_selected %>% mutate(log_weight=log(mass_tidy),
                                            log_tail=log(tail_tidy),
                                            log_wing=log(wing_tidy),
                                            log_beak=log(bill_tidy),
                                            log_w_index=log(wing_index_tidy)) 
bird_traits2<-bird_traits1 %>% select(!c(sociality,mass_tidy,A2_foraging,A2,tail_tidy, wing_tidy,bill_tidy,wing_index_tidy ))
names(bird_traits2)
bird_traits3 <- bird_traits2[order(bird_traits2$species),]
bird_traits4 <- data.frame(bird_traits3[,-1], row.names=bird_traits3[,1])
bird_traits_df<-bird_traits4 %>% rename(Diet=A1,
                                        Foraging=A1_foraging,
                                        BeakLength=log_beak,
                                        WingLength=log_wing,
                                        TailLength=log_tail,
                                        HandWingIndex=log_w_index,
                                        BodyMass=log_weight) 

##### PCA for social species only (independent PCA)
# this exclude species that are not social and creates a pcA with new coordinates 

df_traits_selected1<-df_traits_selected %>% select(species, sociality)

bird_traits_social<- as.data.frame(bird_traits_df) %>% 
  tibble::rownames_to_column("species") %>% 
  left_join(df_traits_selected1, by="species" )%>% 
  filter(sociality=="yes")

bird_traits_social1 <- data.frame(bird_traits_social[,-1], row.names=bird_traits_social[,1])
bird_traits_social2<-bird_traits_social1 %>% select(-c(sociality))

#PCA social species
pca_birds_social_raw <- princomp(bird_traits_social2, cor = TRUE, scores = TRUE)
summary(pca_birds_social_raw) 

#Scores
pca_scores_total_social_raw<- as.data.frame(pca_birds_social_raw$scores) %>% 
  tibble::rownames_to_column("species") %>% 
  left_join(bird_traits3, by= "species")

# loadings
pca_loadings_social_raw<-pca_birds_social_raw$loadings # loading tells you the contribution of each variable to the component
#or
pca_loadings_total_social_raw <- as.data.frame(unclass(pca_birds_social_raw$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

str(pca_loadings_social_raw)

# Eigenvalues
eig.val <- get_eigenvalue(pca_birds_social_raw)
eig.val
res.var <- get_pca_var(pca_birds_social_raw)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(pca_birds_social_raw)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

# scalar to adjust arrow size
scale <-5

pca_loadings_total_scalar_raw<- pca_loadings_total_social_raw %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")),list(scalar=~.*scale))

  # posh names
#mutate(trait = c("Beak length", "Wing length", "Tail length", "Hand-wing Index","Diet","Foraging","Body mass"))
# posh names when including tewo diets or foraging variables
#mutate(trait = c("Beak length", "Wing length", "Tail length", "Wing index", "Diet","Diet2","Body mass"))

# Kernel for social species in the independent pCA 
k_raw_social<-pca_scores_total_social_raw%>% 
  # extract first two principal components
  dplyr::select(.,species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "species")

hpi_raw_social<- Hpi(x = k_raw_social)
k_est_raw_social<- kde(x = k_raw_social, H = hpi_raw_social, compute.cont = TRUE)  
den_raw_social<- list( k_est_raw_social$eval.points[[1]],  k_est_raw_social$eval.points[[2]],
                       k_est_raw_social$estimate)
names(den_raw_social) <- c("x", "y", "z")
dimnames(den_raw_social$z) <- list(den_raw_social$x, den_raw_social$y)
dcc_raw_social<- reshape2::melt(den_raw_social$z)

# 0.5 probability kernel
k_50_raw_social<- k_fun(df = den_raw_social, prob = 0.50) # 0.5 probability kernel
k_95_raw_social<- k_fun(df = den_raw_social, prob = 0.90)# 0.95 probability kernel
k_99_raw_social<- k_fun(df = den_raw_social, prob = 0.99)#0.99 probability kernel

summary(pca_birds_social_raw)
# colour palette
col_pal <- colorRampPalette(c("red", "darkgoldenrod1", "white"))(200)
#col_pal <- colorRampPalette(c("grey30","grey60", "grey80", "white"))(200)

# plot
ggplot(dcc_raw_social, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  #geom_raster(aes(fill = value)) +
  #scale_fill_gradientn(colours = rev(col_pal)) +
  #lines
  geom_hline(yintercept = 0, colour = "gray90",linetype="solid") + 
  geom_vline(xintercept = 0, colour = "gray90",linetype="solid") + 
  # points for species
  geom_point(data =  pca_scores_total_social_raw, aes(x = Comp.1, y = Comp.2), size = 2, alpha = 1, colour ="turquoise4") +
  # probability kernels
  #geom_contour(aes(z = value), breaks = k_50_raw_social, colour = "grey40", size = 0.5) +
  #geom_contour(aes(z = value), breaks = k_95_raw_social, colour = "grey50", size = 0.5) +
  #geom_contour(aes(z = value), breaks = k_99_raw_social, colour = "grey70", size = 0.5) +
  #coord_equal() +
  # add arrows
  geom_segment(data = pca_loadings_total_scalar_raw, aes(x = 0, y = 0, xend = Comp.1_scalar, yend = Comp.2_scalar), size=1.1, lineend = "square") +
  # add dashed arrows ends
  geom_segment(data = pca_loadings_total_scalar_raw, aes(x = 0, y = 0, xend = -Comp.1_scalar, yend = -Comp.2_scalar), lty = 5, colour = "grey50") +
  # add arrow labels
  geom_text(data = pca_loadings_total_scalar_raw, aes(x = Comp.1_scalar, y = Comp.2_scalar, label = trait), nudge_x = c(1, 1, 1, 0, 0.5,0,-0.5,1), nudge_y = c(0.2, 0, 0, -0.5,0.5,0.5,0.2,0),size = 7) +
  # axis labels - see comp_var
  labs(x = "PC1 (%)", y = "PC2 (%)") +
  ylim(-5,5)+
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 30))


# PCA for nonsocial species only ( Independent PCA)------------------------------------------

# this exclude species that social and creates a PCA with new coordinates 

df_traits<-read.csv("data/data_analyses/4.df_traits_manu_birds.csv")
#%>% select(-c(X.1, X))

df_traits_selected<-df_traits %>% 
  select(sociality,species,mass_tidy,bill_tidy,wing_tidy,tail_tidy,wing_index_tidy, A1.x,A2.x,A1.y,A2.y) %>% 
  rename(A1=A1.x,A2=A2.x,A1_foraging=A1.y,A2_foraging=A2.y ) 

#Reformat the structure
bird_traits1<-df_traits_selected %>% mutate(log_weight=log(mass_tidy),
                                            log_tail=log(tail_tidy),
                                            log_wing=log(wing_tidy),
                                            log_beak=log(bill_tidy),
                                            log_w_index=log(wing_index_tidy)) 
bird_traits2<-bird_traits1 %>% select(!c(sociality,mass_tidy,A2_foraging,A2,tail_tidy, wing_tidy,bill_tidy,wing_index_tidy ))
names(bird_traits2)
bird_traits3 <- bird_traits2[order(bird_traits2$species),]
bird_traits4 <- data.frame(bird_traits3[,-1], row.names=bird_traits3[,1])
bird_traits_df<-bird_traits4 %>% rename(Diet=A1,
                                        Foraging=A1_foraging,
                                        BeakLength=log_beak,
                                        WingLength=log_wing,
                                        TailLength=log_tail,
                                        HandWingIndex=log_w_index,
                                        BodyMass=log_weight) 



# Extract the species list 
species_list_sociality<-df_traits_selected %>% select(species, sociality)

bird_traits_nonsocial<- as.data.frame(bird_traits_df) %>% 
  tibble::rownames_to_column("species") %>% 
  left_join(species_list_sociality, by="species" )%>% 
  filter(sociality=="no")

bird_traits_nonsocial1 <- data.frame(bird_traits_nonsocial[,-1], row.names=bird_traits_nonsocial[,1])
bird_traits_nonsocial2 <-bird_traits_nonsocial1 %>% select(-c(sociality))

#PCA social species
pca_birds_nonsocial_raw <- princomp(bird_traits_nonsocial2, cor = TRUE, scores = TRUE)
summary(pca_birds_nonsocial_raw) 

#Scores
pca_scores_total_nonsocial_raw<- as.data.frame(pca_birds_nonsocial_raw$scores) %>% 
  tibble::rownames_to_column("species") %>% 
  left_join(bird_traits3, by= "species")

# loadings
pca_loadings_nonsocial_raw<-pca_birds_nonsocial_raw$loadings # loading tells you the contribution of each variable to the component
#or
pca_loadings_total_nonsocial_raw <- as.data.frame(unclass(pca_birds_nonsocial_raw$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

str(pca_loadings_social_raw)

# Eigenvalues
eig.val <- get_eigenvalue(pca_birds_nonsocial_raw)
eig.val
res.var <- get_pca_var(pca_birds_nonsocial_raw)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(pca_birds_nonsocial_raw)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

# scalar to adjust arrow size
scale <-5

pca_loadings_total_scalar_raw_nonsocial<- pca_loadings_total_nonsocial_raw %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")),list(scalar=~.*scale))

# posh names
#mutate(trait = c("Beak length", "Wing length", "Tail length", "Hand-wing Index","Diet","Foraging","Body mass"))
# posh names when including tewo diets or foraging variables
#mutate(trait = c("Beak length", "Wing length", "Tail length", "Wing index", "Diet","Diet2","Body mass"))

# Kernel for social species in the independent pCA 
k_raw_nonsocial<-pca_scores_total_nonsocial_raw%>% 
  # extract first two principal components
  dplyr::select(.,species, Comp.1, Comp.2) %>% 
  tibble::column_to_rownames(var = "species")

hpi_raw_nonsocial<- Hpi(x = k_raw_nonsocial)
k_est_raw_nonsocial<- kde(x = k_raw_nonsocial, H = hpi_raw_nonsocial, compute.cont = TRUE)  
den_raw_nonsocial<- list( k_est_raw_nonsocial$eval.points[[1]],  k_est_raw_nonsocial$eval.points[[2]],
                       k_est_raw_nonsocial$estimate)
names(den_raw_nonsocial) <- c("x", "y", "z")
dimnames(den_raw_nonsocial$z) <- list(den_raw_nonsocial$x, den_raw_nonsocial$y)
dcc_raw_nonsocial<- reshape2::melt(den_raw_nonsocial$z)

# 0.5 probability kernel
k_50_raw_nonsocial<- k_fun(df = den_raw_nonsocial, prob = 0.50) # 0.5 probability kernel
k_95_raw_nonsocial<- k_fun(df = den_raw_nonsocial, prob = 0.90)# 0.95 probability kernel
k_99_raw_nonsocial<- k_fun(df = den_raw_nonsocial, prob = 0.99)#0.99 probability kernel

summary(pca_birds_nonsocial_raw)
# colour palette
col_pal <- colorRampPalette(c("red", "darkgoldenrod1", "white"))(200)
#col_pal <- colorRampPalette(c("grey30","grey60", "grey80", "white"))(200)

# plot
non_social<-ggplot(dcc_raw_nonsocial, aes(x = Var1, y = Var2)) +
  # coloured probabilty background
  #geom_raster(aes(fill = value)) +
  #scale_fill_gradientn(colours = rev(col_pal)) +
  #lines
  geom_hline(yintercept = 0, colour = "gray90",linetype="solid") + 
  geom_vline(xintercept = 0, colour = "gray90",linetype="solid") + 
  # points for species
  geom_point(data =  pca_scores_total_nonsocial_raw, aes(x = Comp.1, y = Comp.2), size = 2, alpha = 1, colour ="turquoise4") +
  # probability kernels
  #geom_contour(aes(z = value), breaks = k_50_raw_social, colour = "grey40", size = 0.5) +
  #geom_contour(aes(z = value), breaks = k_95_raw_social, colour = "grey50", size = 0.5) +
  #geom_contour(aes(z = value), breaks = k_99_raw_social, colour = "grey70", size = 0.5) +
  #coord_equal() +
  # add arrows
  geom_segment(data = pca_loadings_total_scalar_raw_nonsocial, aes(x = 0, y = 0, xend = Comp.1_scalar, yend = Comp.2_scalar), size=1.1, lineend = "square") +
  # add dashed arrows ends
  geom_segment(data = pca_loadings_total_scalar_raw_nonsocial, aes(x = 0, y = 0, xend = -Comp.1_scalar, yend = -Comp.2_scalar), lty = 5, colour = "grey50") +
  # add arrow labels
  geom_text(data = pca_loadings_total_scalar_raw_nonsocial, aes(x = Comp.1_scalar, y = Comp.2_scalar, label = trait), nudge_x = c(0,0,0,0,0,0,0), nudge_y = c(0,0,0,0,0,0,0),size = 7) +
  # axis labels - see comp_var
  labs(x = "PC1 (47.9%)", y = "PC2 (22.9%)") +
  ylim(-5,5)+
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white', colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 30))

pdf(file="figures/non_social.pdf", width = 8.5, height = 6)
grid.arrange(non_social,nrow=1)
dev.off()

# ## Section 5.1 PCA for flocks ALONG THE GRADIENT  independent of community ----------------------------------

str(flocks)
##Lowlands understory
lowlands_understory<-flocks %>% filter(flock_type_ref=="1")
lowlands_understory<-unique(lowlands_understory$species_taxonomy_SACC_2021)
lowlands_understory<-as.data.frame(lowlands_understory)
species_lowlands<-rename(lowlands_understory,c(species="lowlands_understory"))

# right_join (bird_traits3,species_lowlands, by="species") %>% View() # Double check that we have traits data for all the species 

##_##_##_##_##_##_##_##_##_##_##_##_##
##### lowland PCA Independent
##_##_##_##_##_##_##_##_##_##_##_##_##

# this exclude species that are not social and creates a pcA with new coordinates 
bird_traits_lowlands<- as.data.frame(bird_traits_df) %>% 
  tibble::rownames_to_column("species") %>% 
  right_join(species_lowlands, by="species" )

bird_traits_lowlands1 <- data.frame(bird_traits_lowlands[,-1], row.names=bird_traits_lowlands[,1])
pca_und_lowlands_raw<- princomp(bird_traits_lowlands1, cor = TRUE, scores = TRUE)  
summary(pca_und_lowlands_raw)

#Scores
pca_scores_total_lowlands_raw<- as.data.frame(pca_und_lowlands_raw$scores) %>% 
  tibble::rownames_to_column("species") %>% 
  left_join(df_traits_selected, by= "species")

# loadings
pca_loadings_lowlands_raw<-pca_und_lowlands_raw$loadings # loading tells you the contribution of each variable to the component
#or
pca_loadings_total_lowlands_raw <- as.data.frame(unclass(pca_und_lowlands_raw$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

# scalar to adjust arrow size
scale <-5

pca_loadings_total_scalar_lowlands_raw<- pca_loadings_total_lowlands_raw %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")),list(scalar=~.*scale))

##_##_##_##_##_##_##_##_##_##_##_##_##
##Lowlands canopy
##_##_##_##_##_##_##_##_##_##_##_##_##
l_canopy<-flocks %>% filter(flock_type_ref=="3")
l_canopy<-unique(l_canopy$species_taxonomy_SACC_2021)
l_canopy<-as.data.frame(l_canopy)
species_l_canopy<-rename(l_canopy,c(species="l_canopy"))

##### Independent canopy PCA
# this exclude species that are not social and creates a pcA with new coordinates 
bird_traits_canopy<- as.data.frame(bird_traits_df) %>% 
  tibble::rownames_to_column("species") %>% 
  right_join(species_l_canopy, by="species" )

bird_traits_canopy1 <- data.frame(bird_traits_canopy[,-1], row.names=bird_traits_canopy[,1])
pca_canopy_raw<- princomp(bird_traits_canopy1, cor = TRUE, scores = TRUE)  
summary(pca_canopy_raw)

#Scores
pca_scores_total_canopy_raw<- as.data.frame(pca_canopy_raw$scores) %>% 
  tibble::rownames_to_column("species") %>% 
  left_join(df_traits_selected, by= "species")

# loadings
pca_loadings_canopy_raw<-pca_canopy_raw$loadings # loading tells you the contribution of each variable to the component
#or
pca_loadings_total_canopy_raw <- as.data.frame(unclass(pca_canopy_raw$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

#scalar
scale <-5

pca_loadings_total_scalar_canopy_raw<- pca_loadings_total_canopy_raw %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")),list(scalar=~.*scale))

##_##_##_##_##_##_##_##_##_##_##_##_##
##Bamboo
##_##_##_##_##_##_##_##_##_##_##_##_##

bamboo<-flocks %>% filter(flock_type_ref=="2")
bamboo<-unique(bamboo$species_taxonomy_SACC_2021)
bamboo<-as.data.frame(bamboo)
species_bamboo<-rename(bamboo,c(species="bamboo"))

##### Independednt bamboo PCA
# this exclude species that are not social and creates a pcA with new coordinates 
bird_traits_bamboo<-as.data.frame(bird_traits_df)%>% 
  tibble::rownames_to_column("species") %>% 
  right_join(species_bamboo, by="species" )

bird_traits_bamboo<-na.omit(bird_traits_bamboo[,1:8]) # Make sure there are not NAs


bird_traits_bamboo1 <- data.frame(bird_traits_bamboo[,-1], row.names=bird_traits_bamboo[,1])
pca_bamboo_raw<- princomp(bird_traits_bamboo1, cor = TRUE, scores = TRUE)  
summary(pca_bamboo_raw)

#Scores
pca_scores_total_bamboo_raw<- as.data.frame(pca_bamboo_raw$scores) %>% 
  tibble::rownames_to_column("species") %>% 
  left_join(df_traits_selected, by= "species")

#Loadings
pca_loadings_bamboo_raw<-pca_bamboo_raw$loadings # loading tells you the contribution of each variable to the component
#or
pca_loadings_total_bamboo_raw <- as.data.frame(unclass(pca_bamboo_raw$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

#scalar
scale <-5

pca_loadings_total_scalar_bamboo_raw<- pca_loadings_total_bamboo_raw %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")),list(scalar=~.*scale))

##_##_##_##_##_##_##_##_##_##_##_##_##
##Andes
##_##_##_##_##_##_##_##_##_##_##_##_##

andes<-flocks %>% filter(flock_type_ref=="4")
andes<-unique(andes$species_taxonomy_SACC_2021)
andes<-as.data.frame(andes)
species_andes<-rename(andes,c(species="andes"))

##### Independednt andes PCA
# this exclude species that are not social and creates a pcA with new coordinates 
bird_traits_andes<- as.data.frame(bird_traits_df) %>% 
  tibble::rownames_to_column("species") %>% 
  right_join(species_andes, by="species" )

bird_traits_andes1 <- data.frame(bird_traits_andes[,-1], row.names=bird_traits_andes[,1])
pca_andes_raw<- princomp(bird_traits_andes1, cor = TRUE, scores = TRUE)  
summary(pca_andes_raw)

#Scores
pca_scores_total_andes_raw<- as.data.frame(pca_andes_raw$scores) %>% 
  tibble::rownames_to_column("species") %>% 
  left_join(df_traits_selected, by= "species")

# loadings
pca_loadings_andes_raw<-pca_andes_raw$loadings # loading tells you the contribution of each variable to the component
#or
pca_loadings_total_andes_raw <- as.data.frame(unclass(pca_andes_raw$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

#scalar
scale <-5

pca_loadings_total_scalar_andes_raw<- pca_loadings_total_andes_raw %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")),list(scalar=~.*scale))



##_##_##_##_##_##_##_##_##_##_##_##_##
##High Andes
##_##_##_##_##_##_##_##_##_##_##_##_##
high_andes<-flocks %>% filter(flock_type_ref=="5")
list_flocks_high_andes<-unique(high_andes$species_taxonomy_SACC_2021)
list_flocks_high_andes<-as.data.frame(list_flocks_high_andes)
species_h_andes<-rename(list_flocks_high_andes,c(species="list_flocks_high_andes"))

##### Independednt andes PCA
# this exclude species that are not social and creates a pcA with new coordinates 
bird_traits_h_andes<- as.data.frame(bird_traits_df) %>% 
  tibble::rownames_to_column("species") %>% 
  right_join(species_h_andes, by="species" )

bird_traits_h_andes1 <- data.frame(bird_traits_h_andes[,-1], row.names=bird_traits_h_andes[,1])
pca_h_andes_raw<- princomp(bird_traits_h_andes1, cor = TRUE, scores = TRUE)  
summary(pca_h_andes_raw)

#Scores
pca_scores_total_h_andes_raw<- as.data.frame(pca_h_andes_raw$scores) %>% 
  tibble::rownames_to_column("species") %>% 
  left_join(df_traits_selected, by= "species")

# loadings
pca_loadings_h_andes_raw<-pca_h_andes_raw$loadings # loading tells you the contribution of each variable to the component
#or
pca_loadings_total_h_andes_raw <- as.data.frame(unclass(pca_h_andes_raw$loadings)) %>% 
  tibble::rownames_to_column("trait") %>% 
  dplyr::bind_rows()

#scalar
scale <-5

pca_loadings_total_scalar_h_andes_raw<-pca_loadings_total_h_andes_raw  %>% 
  # rescale for arrow sizes
  dplyr::mutate_at(vars(contains("Comp")),list(scalar=~.*scale))

### Section 6 Hypervolume -----------------------------------------------------
###_###_###_###_###_###
# Computationally heavy 
# Following Mammola et al advice for best practices 
#"Based on the results of this study, it is possible to draw some best practices
#for constructing and comparing n‐dimensional hypervolumes,
#with a special reference to the method by Blonder et al. (2014, 2018)"

#• As long as the correlation between overlap and distance metrics can be low for both datatypes considered (Figure 3), the most diligent
#approach would be to routinely calculate and report at least one overlap and one distance metric
#To facilitate replicability of the analyses, the publication of hypervolumeobjects themselves is also recommended—for example, as an
#RDS file, using the basic R function ‘saveRDS’ (R Core Team, 2017).
#Note : # [Shouldl use scaled traits] Scale the traits first
###_###_###_###_###_###

# Use the data set taht we are using for the PCA and scaled and then filter per elevation
df_traits<-read.csv("data/data_analyses/4.df_traits_manu_birds.csv")

df_traits_selected<-df_traits %>% 
  select(sociality,species,mass_tidy,bill_tidy,wing_tidy,tail_tidy,wing_index_tidy, A1.x,A2.x,A1.y,A2.y) %>% 
  rename(A1=A1.x,A2=A2.x,A1_foraging=A1.y,A2_foraging=A2.y ) 

#Reformat the structure
bird_traits1<-df_traits_selected %>% mutate(log_weight= log(mass_tidy),
                                            log_tail=log(tail_tidy),
                                            log_wing=log(wing_tidy),
                                            log_beak=log(bill_tidy),
                                            log_w_index=log(wing_index_tidy))  

bird_traits2<-bird_traits1 %>% select(!c(sociality,mass_tidy,A2_foraging,A2,bill_tidy,wing_tidy,tail_tidy,wing_index_tidy ))
bird_traits3 <- bird_traits2[order(bird_traits2$species),]
bird_traits_formatted<-bird_traits3 %>% rename(Diet=A1,
                                        Foraging=A1_foraging,
                                        BeakLength=log_beak,
                                        WingLength=log_wing,
                                        TailLength=log_tail,
                                        HandWingIndex=log_w_index,
                                        BodyMass=log_weight) 



#View(bird_traits1)

 
#z- Transform all traits ( scale the traits manually)
# function to z-transform data
scale_z <- function(x){
  (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
}

bird_traits_z <- bird_traits_formatted %>% 
  mutate(mass_z = scale_z(BodyMass)) %>%  # body mass was previously Log10-transformed see "bird_traits3"
  mutate(beak_z = scale_z(BeakLength)) %>% 
  mutate(wing_z = scale_z(WingLength)) %>% 
  mutate(tail_z = scale_z(TailLength )) %>% 
  mutate(wing_ind_z = scale_z(HandWingIndex)) %>% 
  mutate(Foraging_z = scale_z(Foraging))%>% 
  mutate(Diet_z = scale_z(Diet))%>% 
  as.data.frame

bird_traits_z_selected<- bird_traits_z %>% 
  select(species,beak_z,mass_z,tail_z,wing_z,wing_ind_z, Foraging_z, Diet_z) # double check that the species is included as the raw names


#_#_#_#_#_#_#Create traits dataframes of interest_#_#_#_#_#_#
# Select z(scaled traits) for each hypervolume of interest, and 
#remember that species has to be the row names
# Community
bird_traits_z_community<- data.frame(bird_traits_z_selected[,-1], row.names=bird_traits_z_selected[,1])
##_#_#_#_

# Social
flocks_list_species<-flocks_list %>% rename(species=species_taxonomy_SACC_2021) %>% select(-flocking)
bird_traits_z_s<-inner_join(bird_traits_z_selected,flocks_list_species, by="species")
#
bird_traits_z_social<- data.frame(bird_traits_z_s[,-1], row.names=bird_traits_z_s[,1])
#_#_#_#_

# Non-Social
flocks_list_species<-flocks_list %>% rename(species=species_taxonomy_SACC_2021) %>% select(-flocking)
bird_traits_z_ns<-anti_join(bird_traits_z_selected,flocks_list_species, by="species")
#
bird_traits_z_nonsocial<- data.frame(bird_traits_z_ns[,-1], row.names=bird_traits_z_ns[,1])
#_#_#_#_

##Lowlands understory
species_lowlands
bird_traits_z_l<-inner_join(bird_traits_z_selected,species_lowlands, by="species")
bird_traits_z_und_lowlands<-data.frame(bird_traits_z_l[,-1], row.names=bird_traits_z_l[,1])

##Lowlands canopy
species_l_canopy
bird_traits_z_c<-inner_join(bird_traits_z_selected,species_l_canopy, by="species")
bird_traits_z_canopy<-data.frame(bird_traits_z_c[,-1], row.names=bird_traits_z_c[,1])

##Bamboo
species_bamboo
bird_traits_z_b<-inner_join(bird_traits_z_selected,species_bamboo, by="species")
bird_traits_z_bamboo<-data.frame(bird_traits_z_b[,-1], row.names=bird_traits_z_b[,1])
##Andes
species_andes
bird_traits_z_a<-inner_join(bird_traits_z_selected,species_andes, by="species")
bird_traits_z_andes<-data.frame(bird_traits_z_a[,-1], row.names=bird_traits_z_a[,1])

##High Andes
species_h_andes
bird_traits_z_ha<-inner_join(bird_traits_z_selected,species_h_andes, by="species")
bird_traits_z_h_andes<-data.frame(bird_traits_z_ha[,-1], row.names=bird_traits_z_ha[,1])

####Warning double check similarity withusing the Components instead of the traits example from gomez

###_###_###_##_
# Observed hypervolumes ---------------------------------------------------
###_###_###_##_
#package(hypervolume)

# Hypervolumes
# The value of sd.count should generally not be modified, and should never be decreased from its default value of 3
#samples.per.point: This parameter determines the number of uniform random points used to represent values of h(x). 
#The default value of this parameter is chosen via a heuristic approach to increase with the square root of the dimensionality of the analysis, in order to
#provide more robust estimates of higher-dimensional shapes that may have complex forms. the value is chosen on a per-point
#basis, so that hypervolumes with different number of data points have the same overall sampling effort. 
#The investigator can further increase this parameter if desired. Larger values of this parameter produce more robust results at higher computational cost.

# Comparing hypervolumnes 
#When and how can I compare two hypervolumes? 
#Hypervolumes can only be compared if they are constructed using the same axes (both number of axes and identify of axes). 
#The volume of a hypervolume is in units with dimensionality equal to the dimensionality of the axes; while it appears to be just a scalar number its units will change. 
#Thus a 3-dimensional volume of '11.2' is not comparable to a 4-dimensional volume of '65.8': it is neither smaller nor larger, but simply incomparable.
#Some care should also be taken when comparing the volumes of hypervolumes of the same dimensionality 
#if a fixed kernel bandwidth was used to construct them. 
#In this case, the volume of the hypervolume is approximately linearly proportional to the number of observations in the dataset. 
#This is because each new data point contributes approximately the same amount of volume, unless it overlaps with a previous data point. 
#The issue is then that the largest hypervolumes will be those constructed from the largest number of data points. 
#This may reflect the true structure of your data, but if it does not, you should instead proceed with a null-modeling procedure where you compare the observed hypervolume to that of a distribution of null hypervolumes constructed by resampling an identical number of observations. 
#Instead of reporting a raw hypervolume you can report a deviation hypervolume, e.g. a z-score.

#Hypervolume community
set.seed(9) # Reproducible random number generator 
obs_community_hypervolume<-hypervolume_svm(bird_traits_z_community, 
                                           name="Community")
                                           

saveRDS(obs_community_hypervolume, "data/data_analyses/obs_community_hypervolume_7_traits_svm.rds") #SAVED
obs_community_hypervolume<-readRDS("data/data_analyses/obs_community_hypervolume_7_traits_svm.rds")

obs_community_hypervolume_gaussian<-hypervolume_gaussian(bird_traits_z_community, 
                                           name="Community",
                                           weight = NULL,
                                           sd.count= 3,
                                           quantile.requested = 0.99,
                                           quantile.requested.type = "probability",
                                           #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
                                           kde.bandwidth= estimate_bandwidth(bird_traits_z_community, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
                                           )
 saveRDS(obs_community_hypervolume_gaussian, "data/data_analyses/obs_community_hypervolume_7_traits_gaussian.rds") #SAVED
 
 summary( obs_community_hypervolume_gaussian)
 
 #Hypervolume social
 
set.seed(9) #
obs_social_hypervolume<-hypervolume_svm(bird_traits_z_social,
                                        name="All flocks")
saveRDS(obs_social_hypervolume, "data/data_analyses/obs_social_hypervolume_7_traits_svm.rds") #SAVED
obs_social_hypervolume<-readRDS("data/data_analyses/obs_social_hypervolume_7_traits_svm.rds")

obs_social_hypervolume_gaussian<-hypervolume_gaussian(bird_traits_z_social, 
                                                      name="All flocks",
                                                      weight = NULL,
                                                      sd.count= 3,
                                                      quantile.requested = 0.99,
                                                      quantile.requested.type = "probability",
                                                      #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
                                                      kde.bandwidth= estimate_bandwidth(bird_traits_z_social, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)
saveRDS(obs_social_hypervolume_gaussian, "data/data_analyses/obs_social_hypervolume_7_traits_gaussian.rds") #SAVED

summary(obs_social_hypervolume_gaussian)
hypervolume_variable_importance(obs_social_hypervolume_gaussian)

#Hypervolume nonsocial
set.seed(9)
obs_nonsocial_hypervolume<-hypervolume_svm(bird_traits_z_nonsocial, 
                                                name="non social")

saveRDS(obs_nonsocial_hypervolume, "data/data_analyses/obs_nonsocial_hypervolume_7_traits_svm.rds") #SAVED
obs_nonsocial_hypervolume<-readRDS("data/data_analyses/obs_nonsocial_hypervolume_7_traits_z.rds")

obs_nonsocial_hypervolume_gaussian<-hypervolume_gaussian(bird_traits_z_nonsocial, 
                                                         name="non social",
                                                         weight = NULL,
                                                         sd.count= 3,
                                                         quantile.requested = 0.99,
                                                         quantile.requested.type = "probability",
                                                         #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5))
                                                         kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
                                                         )
saveRDS(obs_nonsocial_hypervolume_gaussian, "data/data_analyses/obs_nonsocial_hypervolume_7_traits_gaussian.rds") #SAVED
obs_nonsocial_hypervolume_gaussian<-readRDS("data/data_analyses/obs_nonsocial_hypervolume_7_traits_gaussian.rds")
summary(obs_nonsocial_hypervolume_gaussian)

#hypervolume lowlands 61 1.855701/93.911370/88

set.seed(9)
obs_lowlands_hypervolume<-hypervolume_svm(bird_traits_z_und_lowlands, name="Lowlands understory flocks")#SAVED
saveRDS(obs_lowlands_hypervolume, "data/data_analyses/obs_lowlands_hypervolume_7_traits_svm.rds")

obs_lowlands_hypervolume_gaussian<-hypervolume_gaussian(bird_traits_z_und_lowlands, 
                                                         name="lowlands",
                                                         weight = NULL,
                                                         sd.count= 3,
                                                         quantile.requested = 0.99,
                                                         quantile.requested.type = "probability",
                                                         kde.bandwidth= estimate_bandwidth(bird_traits_z_und_lowlands, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
                                                         #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)

saveRDS(obs_lowlands_hypervolume_gaussian, "data/data_analyses/obs_lowlands_hypervolume_7_traits_gaussian.rds") # SAVED





#hypervolume canopy 110 7.543544/150/147

set.seed(9)
obs_canopy_hypervolume<-hypervolume_svm(bird_traits_z_canopy, name="Lowlands canopy flocks") #SAVED
saveRDS(obs_canopy_hypervolume, "data/data_analyses/obs_canopy_hypervolume_7_traits_svm.rds")


obs_canopy_hypervolume_gaussian<-hypervolume_gaussian(bird_traits_z_canopy, 
                                                        name="canopy",
                                                        weight = NULL,
                                                        sd.count= 3,
                                                        quantile.requested = 0.99,
                                                        quantile.requested.type = "probability",
                                                        #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5))
                                                        kde.bandwidth= estimate_bandwidth(bird_traits_z_canopy, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)

saveRDS(obs_canopy_hypervolume_gaussian, "data/data_analyses/obs_canopy_hypervolume_7_traits_gaussian.rds") # SAVED


#hypervolume bamboo 43  0.603039/66.611306/66.93

set.seed(9)
obs_bamboo_hypervolume<-hypervolume_svm(bird_traits_z_bamboo, name="Lowlands Bamboo flocks")
saveRDS(obs_bamboo_hypervolume, "data/data_analyses/obs_bamboo_hypervolume_7_traits_svm.rds")#SAVED

obs_bamboo_hypervolume_gaussian<-hypervolume_gaussian(bird_traits_z_bamboo, 
                                                      name="bamboo",
                                                      weight = NULL,
                                                      sd.count= 3,
                                                      quantile.requested = 0.99,
                                                      quantile.requested.type = "probability",
                                                      #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5))
                                                      kde.bandwidth= estimate_bandwidth(bird_traits_z_bamboo, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)

 saveRDS(obs_bamboo_hypervolume_gaussian, "data/data_analyses/obs_bamboo_hypervolume_7_traits_gaussian.rds") # SAVED


#hypervolume montane 119 4.514920/ 149.6  /151.04
set.seed(9)


obs_montane_hypervolume<-hypervolume_svm(bird_traits_z_andes, name="Montane flocks")
saveRDS(obs_montane_hypervolume, "data/data_analyses/obs_montane_hypervolume_7_traits_svm.rds") # SAVED

obs_montane_hypervolume_gaussian<-hypervolume_gaussian(bird_traits_z_andes, 
                                                      name="montane",
                                                      weight = NULL,
                                                      sd.count= 3,
                                                      quantile.requested = 0.99,
                                                      quantile.requested.type = "probability",
                                                      #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5))
                                                      kde.bandwidth= estimate_bandwidth(bird_traits_z_andes, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)
saveRDS(obs_montane_hypervolume_gaussian, "data/data_analyses/obs_montane_hypervolume_7_traits_gaussian.rds") # SAVED


#hypervolume high montane  75  1.805658/ /104.26   svm/gau/gau 2000
set.seed(9) 
obs_h_montane_hypervolume<-hypervolume_svm(bird_traits_z_h_andes, name="High montane flocks")
saveRDS(obs_h_montane_hypervolume, "data/data_analyses/obs_h_montane_hypervolume_7_traits_svm.rds") 


obs_h_montane_hypervolume_gaussian<-hypervolume_gaussian(bird_traits_z_h_andes, 
                                                       name="h_montane",
                                                       weight = NULL,
                                                       sd.count= 3,
                                                       quantile.requested = 0.99,
                                                       quantile.requested.type = "probability",
                                                       #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5))
                                                       kde.bandwidth= estimate_bandwidth(bird_traits_z_h_andes, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)
saveRDS(obs_h_montane_hypervolume_gaussian, "data/data_analyses/obs_h_montane_hypervolume_7_traits_gaussian.rds") # SAVED

# Estimating variable importance for each hypervolume
#Larger values indicate that a variable makes a proportionally higher contribution to the overall volume.

varimp_community<-hypervolume_variable_importance(obs_community_hypervolume_gaussian,verbose=TRUE)
varimp_social<-hypervolume_variable_importance(obs_social_hypervolume_gaussian,verbose=TRUE)  
varimp_nonsocial<-hypervolume_variable_importance(obs_nonsocial_hypervolume_gaussian,verbose=TRUE)

varimp_h_montane<-hypervolume_variable_importance(obs_h_montane_hypervolume_gaussian,verbose=TRUE)
varimp_montane<-hypervolume_variable_importance(obs_montane_hypervolume_gaussian,verbose=TRUE)

varimp_canopy<-hypervolume_variable_importance(obs_canopy_hypervolume_gaussian,verbose=TRUE)
varimp_lowlands<-hypervolume_variable_importance(obs_lowlands_hypervolume_gaussian,verbose=TRUE)
varimp_bamboo<-hypervolume_variable_importance(obs_bamboo_hypervolume_gaussian,verbose=TRUE)

varimp_montane=varimp_h_montane
varimp_lowlands=varimp1
barplot(varimp,ylab='Importance',xlab='Variable')
barplot(varimp1,ylab='Importance',xlab='Variable')

###_###_###_###_###_###_###_
# Raw intersections
#Calculating overlap between the hypervolumes total=all social shared=intersection
###_###_###_###_###_###_###_


hv_list_general=hypervolume_join(obs_nonsocial_hypervolume_gaussian,obs_social_hypervolume_gaussian)

all_species_intersection<-hypervolume_set_n_intersection(hv_list_general, num.points.max = NULL,
                                                         verbose = TRUE, distance.factor = 1, check.hyperplane = FALSE)


# Finds the intersection of multiple hypervolumes.
#Using this function is likely faster and more accurate than iteratively applying hypervolume_set to hypervolume pairs, as this function does
#not iteratively perform downsampling.Stores all the points from the input hypervolumes in a single set. Then uses the inclusion test
#approach to identify and store points from this set that are within each individual resampled hypervolume, successively. All the points that are common to all the tests are grouped, resampled and
#used to generate the hypervolume corresponding to the intersection.


###_###_###_##_
#Hypervolume overlaps between all flock types
###_###_###_##_
hv_list_flocks_gradient=hypervolume_join(obs_bamboo_hypervolume_gaussian,obs_canopy_hypervolume_gaussian,obs_lowlands_hypervolume_gaussian,obs_montane_hypervolume_gaussian,obs_h_montane_hypervolume_gaussian)
all_flocks_intersection<-hypervolume_set_n_intersection(hv_list_flocks_gradient, num.points.max = NULL,verbose = TRUE, distance.factor = 1, check.hyperplane = FALSE)
raw_all_flocks_intersection=all_flocks_intersection

saveRDS(raw_all_flocks_intersection, "data/data_analyses/raw_all_flocks_intersection.rds") 

# summarise volumes example
hypervolume::get_volume(hv_set)
plot(hv1, countour.type='kde')
plot(hypervolume_join(hv1,hv2, hv3),contour.lwd=2,num.points.max.data = 442336 , num.points.max.random = 50000,cex.axis=1,colors = c("red","blue","green"),show.random = FALSE,show.centroid = TRUE,show.legend =TRUE)

hv1 = hypervolume_box(penguins_adelie,name='Adelie')
hv2 = hypervolume_box(penguins_chinstrap,name='Chinstrap')
hv_set <- hypervolume_set(hv1, hv2, check.memory=FALSE)
hypervolume_overlap_statistics(hv_set)

#ow do I animate or save a 3D plot? 
#Two-dimensional plots can be saved using standard R commands. 
#However three-dimensional plots use the RGL library and must be saved differently.
#To save a snapshot of a plot, run your normal plotting commands, then: rgl.bringtotop(); rgl.snapshot('test.png'). 
#If you would instead like to save an animated GIF of a rotating hypervolume, you can run movie3d(spin3d(),duration=5,movie='mymovie',dir='./').


# read: hypervolumes
obs_community_hypervolume <- readRDS("data/data_analyses/obs_community_hypervolume_7_traits_gaussian.rds")
obs_nonsocial_hypervolume <- readRDS("data/data_analyses/obs_nonsocial_hypervolume_7_traits_gaussian.rds")
obs_social_hypervolume <- readRDS("data/data_analyses/obs_social_hypervolume_7_traits_gaussian.rds")
#obs_lowlands_hypervolume<- readRDS("data/data_analyses/obs_lowlands_hypervolume_7_traits_gaussian.rds.rds")
#obs_canopy_hypervolume<- readRDS("data/data_analyses/obs_canopy_hypervolume_7_traits_gaussian.rds")
#obs_bamboo_hypervolume<- readRDS("data/data_analyses/obs_bamboo_hypervolume_7_traits_gaussian.rds")
#obs_montane_hypervolume<- readRDS("data/data_analyses/obs_montane_hypervolume_7_traits_gaussian.rds.rds")
#obs_h_montane_hypervolume<- readRDS("data/data_analyses/obs_h_montane_hypervolume_7_traits_gaussian.rds")

###_###_###_##_
#Hypervolume overlaps  between pairs of flock types or community vs flocks
###_###_###_##_

# set raw hypervolumes for comparison

set.seed(9)

social_nonsocial_set<-hypervolume_set(obs_nonsocial_hypervolume_gaussian,obs_social_hypervolume_gaussian,check.memory = FALSE)
social_nonsocial_overlap<-hypervolume_overlap_statistics(social_nonsocial_set)
summary(social_nonsocial_overlap)
social_nonsocial_set

social_community_set<-hypervolume_set(obs_community_hypervolume_gaussian,obs_social_hypervolume_gaussian,check.memory = FALSE)
social_community_overlap<-hypervolume_overlap_statistics(social_community_set)
summary(social_community_set)
social_community_set

nonsocial_community_set<-hypervolume_set(obs_community_hypervolume_gaussian,obs_nonsocial_hypervolume_gaussian,check.memory = FALSE)
nonsocial_community_overlap<-hypervolume_overlap_statistics(nonsocial_community_set)
summary(social_community_set)
social_community_set


understory_canopy_set<-hypervolume_set(obs_lowlands_hypervolume_gaussian,obs_canopy_hypervolume_gaussian, check.memory = FALSE)
understory_canopy_overlap<-hypervolume_overlap_statistics(understory_canopy_set)
summary(understory_canopy_overlap)
understory_canopy_set

lowlands_andes_set<-hypervolume_set(obs_lowlands_hypervolume_gaussian,obs_montane_hypervolume_gaussian, check.memory = FALSE)
hypervolume_overlap_statistics(lowlands_andes_set)
summary(lowlands_andes_set)

lowlands_h_andes_set<-hypervolume_set(obs_lowlands_hypervolume_gaussian,obs_h_montane_hypervolume_gaussian,check.memory = FALSE)
hypervolume_overlap_statistics(lowlands_h_andes_set)
summary(lowlands_h_andes_set)

understory_bamboo_set<-hypervolume_set(obs_lowlands_hypervolume_gaussian,obs_bamboo_hypervolume_gaussian, check.memory = FALSE)
hypervolume_overlap_statistics(understory_bamboo_set)
summary(understory_bamboo_set)

bamboo_canopy_set<-hypervolume_set(obs_bamboo_hypervolume_gaussian,obs_canopy_hypervolume_gaussian, check.memory = FALSE)
bamboo_canopy_overlap<-hypervolume_overlap_statistics(bamboo_canopy_set)
summary(bamboo_canopy_set)

canopy_andes_set<-hypervolume_set(obs_montane_hypervolume_gaussian,obs_canopy_hypervolume_gaussian,check.memory = FALSE)
hypervolume_overlap_statistics(canopy_andes_set)
summary(canopy_andes_set)

canopy_h_andes_set<-hypervolume_set(obs_h_montane_hypervolume_gaussian,obs_canopy_hypervolume_gaussian,check.memory = FALSE)
hypervolume_overlap_statistics(canopy_h_andes_set)
summary(canopy_andes_set)

bamboo_montane_set<-hypervolume_set(obs_bamboo_hypervolume_gaussian,obs_montane_hypervolume_gaussian, check.memory = FALSE)
bamboo_montane_overlap<-hypervolume_overlap_statistics(bamboo_montane_set) # double check this values on the table
summary(bamboo_montane_set)

bamboo_h_andes_set<-hypervolume_set(obs_bamboo_hypervolume_gaussian,obs_h_montane_hypervolume_gaussian, check.memory = FALSE)
bamboo_h_andes_overlap<-hypervolume_overlap_statistics(bamboo_h_andes_set)
summary(bamboo_h_andes_set)

montane_h_andes_set<-hypervolume_set(obs_montane_hypervolume_gaussian,obs_h_montane_hypervolume_gaussian, check.memory = FALSE)
montane_h_andes_overlap<-hypervolume_overlap_statistics(montane_h_andes_set)
summary(montane_h_andes_set)

hypervolume_overlap_test(obs_montane_hypervolume_gaussian,obs_h_montane_hypervolume_gaussian)

#Notes: #jaccard Jaccard similarity (volume of intersection of 1 and 2 divided by volume of union of 1 and 2) # used for the analyses
#sorensen Sorensen similarity (twice the volume of intersection of 1 and 2 divided by volume of 1 plus volume of 2)
#frac_unique_1 Unique fraction 1 (volume of unique component of 1 divided by volume of 1))
#frac_unique_2 Unique fraction 2 (volume of unique component of 2 divided by volume of 2))

# summarise volumes
hypervolume::get_volume(montane_h_andes_set)
hypervolume::get_volume(bamboo_h_andes_set)
hypervolume::get_volume(bamboo_montane_set)
hypervolume::get_volume(canopy_h_andes_set)
hypervolume::get_volume(lowlands_andes_set)

hypervolume::get_volume(understory_canopy_set)

## Generate a cluster analyses of the Functional differences between flocks 

#method Jaccard aglomeration method used "average" 

hypervolume()


sim<-read.csv("data/data_analyses/5.flock_dissimilarity_jaccard.csv", header=TRUE, row.names=1)

# USING DENDEXTEND PACKAGE TO REDO THE CLUSTER ANALYSES

library(dendextend)
hypervolume_overlap_test()

sim<-read.csv("data/data_analyses/1.similarity_matrix_flocks.csv", header=TRUE, row.names=1)
sim1<-as.matrix(sim)
true.dist.mat <- as.dist(sim1) # coerce a object to a distance matrix. The smaller the distance teh smaller the dissimilarity
hc<-hclust(true.dist.mat,method = "average")
dend <- as.dendrogram(hc)


layout(t(c(1,1,1,2,2)))

layout(t(c(1,1,1,1,1)))

dend %>% set("branches_k_color",value=c(3,5,8,4,2,6),k=6) %>%
  set("labels_cex",1) %>% 
  set("labels_cex",1) %>% 
  plot(horiz = TRUE,main = "Cluster of mixed-species flocks along Manu elevational gradient", size=2) 
dend %>% rect.dendrogram(k=6, horiz = TRUE,
                         border = 8, lty = 5, lwd = 2)


par(mfrow = c(1,1))




par(mfrow = c(1,2))
dend %>% set("branches_k_color", k = 6) %>% plot(main = "Nice defaults")
dend %>% set("branches_k_color", value = 3:9, k = 6) %>% 
  plot(main = "Controlling branches' colors\n(via clustering)")


dend15 %>% rect.dendrogram(k=3, 
                           border = 8, lty = 5, lwd = 2)







####more clusters


dend %>% highlight_branches_col %>% plot(main = "Coloring branches")



dend %>% set("labels_cex", 2) %>% set("labels_col", value = c(1,2,3,4,5,6), k=6) %>% 
  plot(main = "Color labels \nper cluster")
abline(h = 2, lty = 2)


sim<-read.csv("data/data_analyses/5.flock_dissimilarity_jaccard.csv", header=TRUE, row.names=1)
true.dist.mat <- as.dist(sim) # coerce a object to a distance matrix. The smaller the distance teh smaller the dissimilarity
hc<-hclust(true.dist.mat,method = "complete")
dend <- as.dendrogram(hc)
dend2 <- as.dendrogram(hc)

par(mfrow = c(1,1))

dend %>% highlight_branches_col %>% plot(main = "Coloring branches")
dend %>% highlight_branches %>% plot(main = "Emphasizing color\n and line-width")
dend %>% highlight_branches_lwd %>% plot(main = "Emphasizing line-width")

dl <- dendlist(dend, dend2)
tanglegram(dl, sort = TRUE, common_subtrees_color_lines = FALSE, highlight_distinct_edges  = FALSE, highlight_branches_lwd = FALSE)


# Hypervolume resampling --------------------------------------------------

# Correcting for sample size 
# Since hypervolumes are influenced by sample size, and the number of species in each flock type varies strongly accross flcoks types. 
#For the comparison of individual hypervolumens, we will select random samples of species from teh larger of each paried hypervolume to match
# the number os species in the smaller hypervolume. And run each comparative anlyses in 50 random subsets of teh larger hypervolume  usingthe function hypervolume_resample.
# Then we calculated  the statistics including jaccard overalp with the 95% Confidence interval

# Create the hypervolume for social species with 2000 samples per point


#Hypervolume community
set.seed(9) # Reproducible random number generator 
obs_community_hypervolume_gaussian_1<-hypervolume_gaussian(bird_traits_z_community, 
                                                        name="Community_2000_samples",
                                                        samples.per.point =2000,
                                                         weight = NULL,
                                                         sd.count= 3,
                                                         quantile.requested = 0.99,
                                                         quantile.requested.type = "probability",
                                                         #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
                                                         kde.bandwidth= estimate_bandwidth(bird_traits_z_community, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)
saveRDS(obs_community_hypervolume_gaussian_1, "data/data_analyses/obs_community_hypervolume_gaussian_1_2000samples.rds") #SAVED
obs_community_hypervolume_gaussian_1<-readRDS("data/data_analyses/obs_community_hypervolume_gaussian_1_2000samples.rds")
#Hypervolume social

set.seed(9) #
obs_social_hypervolume_gaussian_1<-hypervolume_gaussian(bird_traits_z_social, 
                                                      name="All flocks_2000_samples",
                                                      samples.per.point =2000,
                                                      weight = NULL,
                                                      sd.count= 3,
                                                      quantile.requested = 0.99,
                                                      quantile.requested.type = "probability",
                                                      #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
                                                      kde.bandwidth= estimate_bandwidth(bird_traits_z_social, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)
saveRDS(obs_social_hypervolume_gaussian, "data/data_analyses/obs_social_hypervolume_gaussian_1_2000samples.rds") #SAVED
obs_social_hypervolume_gaussian_1<-readRDS("data/data_analyses/obs_social_hypervolume_gaussian_1_2000samples.rds")

#Hypervolume nonsocial
set.seed(9)
obs_nonsocial_hypervolume_gaussian1<-hypervolume_gaussian(bird_traits_z_nonsocial, 
                                                         name="non social_2000_samples",
                                                         samples.per.point =2000,
                                                         weight = NULL,
                                                         sd.count= 3,
                                                         quantile.requested = 0.99,
                                                         quantile.requested.type = "probability",
                                                         #kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5))
                                                         kde.bandwidth= estimate_bandwidth(bird_traits_z_nonsocial, method="fixed", value = c(0.3,0.3,0.3,0.3,0.3,0.3,0.3))
)
saveRDS(obs_nonsocial_hypervolume_gaussian, "data/data_analyses/obs_nonsocial_hypervolume_gaussian_1_2000samples.rds") #SAVED
obs_nonsocial_hypervolume_gaussian1<-readRDS("data/data_analyses/obs_nonsocial_hypervolume_gaussian_1_2000samples.rds")

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

# Hypervolume with MATCHING sample sizes ----------------------------------
# Commmunity ves social and non social 
# This creates objects in my directory that I can use from the paths 

# We will us the observed for the pairwise comparison 
#hypervolumes_randomized_social=hypervolume_resample("social_resample", 
                                                        obs_social_hypervolume_gaussian_1, 
                                                        method="bootstrap", 
                                                        n = 1, 
                                                        points_per_resample = "sample_size", 
                                                        cores = 6,
                                                        verbose = TRUE, 
                                                        weight_func = NULL)

hypervolumes_randomized_social="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/social_resample"

hypervolumes_randomized_community_match_social=hypervolume_resample("community_resample_match_social", 
                                                    obs_community_hypervolume_gaussian_1, 
                                                    method="bootstrap", 
                                                    n = 50, 
                                                    points_per_resample=269, 
                                                    cores = 6,
                                                    verbose = TRUE, 
                                                    weight_func = NULL)

hypervolumes_randomized_community_match_social="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/community_resample_match_social"

hypervolumes_randomized_nonsocial_match_social=hypervolume_resample("nonsocial_resample_match_social", 
                                                                    obs_nonsocial_hypervolume_gaussian1, 
                                                                    method="bootstrap", 
                                                                    n = 50, 
                                                                    points_per_resample=269, 
                                                                    cores = 6,
                                                                    verbose = TRUE, 
                                                                    weight_func = NULL)

hypervolumes_randomized_nonsocial_match_social= "/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/nonsocial_resample_match_social"

#hypervolumes_randomized_nonsocial=hypervolume_resample("nonsocial_resample", 
                                                       obs_nonsocial_hypervolume_gaussian1, 
                                                       method="bootstrap", 
                                                       n =1, 
                                                       points_per_resample="sample_size", 
                                                       cores = 6,
                                                       verbose = TRUE, 
                                                       weight_func = NULL)

hypervolumes_randomized_nonsocial="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/nonsocial_resample"

hypervolumes_randomized_community_match_nonsocial=hypervolume_resample("community_resample_match_non_social", 
                                                                    obs_community_hypervolume_gaussian_1, 
                                                                    method="bootstrap", 
                                                                    n = 100, 
                                                                    points_per_resample=408, 
                                                                    cores = 6,
                                                                    verbose = TRUE, 
                                                                    weight_func = NULL)


hypervolumes_randomized_community_match_nonsocial="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/community_resample_match_non_social"

#Understory 
# create the distributions of hypervolumes for flock type understory (matching the understory sample size 61)

# make sure you keep the sample size identical to what you want to compare with 
#hypervolumes_randomized_understory=hypervolume_resample("understory_resample", 
                                                         obs_lowlands_hypervolume_gaussian_1, 
                                                         method="bootstrap", 
                                                         n = 1, 
                                                         points_per_resample = "sample_size", 
                                                         cores = 6,
                                                         verbose = TRUE, 
                                                         weight_func = NULL)
#hypervolumes_randomized_understory="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/understory_resample"

hypervolumes_randomized_canopy_match_understory=hypervolume_resample("canopy_resample_match_understory", 
                                                    obs_canopy_hypervolume_gaussian_1, 
                                                    method="bootstrap", 
                                                    n = 50, 
                                                    points_per_resample =61,  
                                                    cores = 6,
                                                    verbose = TRUE, 
                                                    weight_func = NULL)
hypervolumes_randomized_canopy_match_understory="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/canopy_resample_match_understory"

hypervolumes_randomized_montane_match_understory=hypervolume_resample("montane_resample_match_understory", 
                                                    obs_montane_hypervolume_gaussian_1, 
                                                    method="bootstrap", 
                                                    n = 50, 
                                                    points_per_resample =61,  # make sure you keep the sample size identical to what you want to compare with 
                                                    cores = 6,
                                                    verbose = TRUE, 
                                                    weight_func = NULL)
hypervolumes_randomized_montane_match_understory="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/montane_resample_match_understory"

hypervolumes_randomized_h_montane_match_understory=hypervolume_resample("h_montane_resample_match_understory", 
                                                     obs_h_montane_hypervolume_gaussian_1, 
                                                     method="bootstrap", 
                                                     n = 50, 
                                                     points_per_resample =61,  # make sure you keep the sample size identical to what you want to compare with 
                                                     cores = 6,
                                                     verbose = TRUE, 
                                                     weight_func = NULL)
hypervolumes_randomized_h_montane_match_understory="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/h_montane_resample_match_understory"
#Understory stats
# create the distributions of hypervolumes for flock type understory (matching the bambo sample size 43)

#### ###Bamboo
hypervolumes_randomized_bamboo=hypervolume_resample("bamboo_resample", 
                                                    obs_bamboo_hypervolume_gaussian_1, 
                                                    method="bootstrap", 
                                                    n = 1, 
                                                    points_per_resample ="sample_size", 
                                                    cores = 6,
                                                    verbose = TRUE, 
                                                    weight_func = NULL)
hypervolumes_randomized_bamboo="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/bamboo_resample"

#### ###Bamboo to compare will all flocks in the sommunity. 
# This is to be able to compare eauql numbers from all the flocks 
hypervolumes_randomized_bamboo_overall=hypervolume_resample("bamboo_resample_overall", 
                                                                     obs_bamboo_hypervolume_gaussian_1, 
                                                                     method="bootstrap", 
                                                                     n = 50, 
                                                                     points_per_resample ="sample_size", 
                                                                     cores = 6,
                                                                     verbose = TRUE, 
                                                                     weight_func = NULL)
hypervolumes_randomized_bamboo_overall="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/bamboo_resample_overall"

hypervolumes_randomized_understory_match_bamboo=hypervolume_resample("understory_resample_match_bamboo", 
                                                        obs_lowlands_hypervolume_gaussian_1, 
                                                        method="bootstrap", 
                                                        n = 50, 
                                                        points_per_resample =43, 
                                                        cores = 6,
                                                        verbose = TRUE, 
                                                        weight_func = NULL)

hypervolumes_randomized_understory_match_bamboo="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/understory_resample_match_bamboo"

hypervolumes_randomized_canopy_match_bamboo=hypervolume_resample("canopy_resample_match_bamboo", 
                                                                     obs_canopy_hypervolume_gaussian_1, 
                                                                     method="bootstrap", 
                                                                     n = 50, 
                                                                     points_per_resample =43,  
                                                                     cores = 6,
                                                                     verbose = TRUE, 
                                                                     weight_func = NULL)

hypervolumes_randomized_canopy_match_bambo="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/canopy_resample_match_bamboo"

hypervolumes_randomized_montane_match_bamboo=hypervolume_resample("montane_resample_match_bamboo", 
                                                                      obs_montane_hypervolume_gaussian_1, 
                                                                      method="bootstrap", 
                                                                      n = 50, 
                                                                      points_per_resample =43,  # make sure you keep the sample size identical to what you want to compare with 
                                                                      cores = 6,
                                                                      verbose = TRUE, 
                                                                      weight_func = NULL)

hypervolumes_randomized_montane_match_bamboo="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/montane_resample_match_bamboo"

hypervolumes_randomized_h_montane_match_bamboo=hypervolume_resample("h_montane_resample_match_bamboo", 
                                                                         obs_h_montane_hypervolume_gaussian_1, 
                                                                         method="bootstrap", 
                                                                         n = 50, 
                                                                         points_per_resample =43,  # make sure you keep the sample size identical to what you want to compare with 
                                                                         cores = 6,
                                                                         verbose = TRUE, 
                                                                         weight_func = NULL)

hypervolumes_randomized_h_montane_match_bamboo="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/h_montane_resample_match_bamboo"
                                                        
# Canopy

#hypervolumes_randomized_canopy= hypervolume_resample("canopy_resample", 
                                                                     obs_canopy_hypervolume_gaussian_1, 
                                                                     method="bootstrap", 
                                                                     n = 1, 
                                                                     points_per_resample ="sample_size",  
                                                                     cores = 6,
                                                                     verbose = TRUE, 
                                                                     weight_func = NULL)



hypervolumes_randomized_canopy= "~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/canopy_resample"


hypervolumes_randomized_canopy_overall=hypervolume_resample("canopy_resample_overall", 
                                                            obs_canopy_hypervolume_gaussian_1, 
                                                            method="bootstrap", 
                                                            n = 50, 
                                                            points_per_resample ="sample_size", 
                                                            cores = 6,
                                                            verbose = TRUE, 
                                                            weight_func = NULL)
hypervolumes_randomized_canopy_overall="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/canopy_resample_overall"


hypervolumes_randomized_montane_match_canopy=hypervolume_resample("montane_resample_match_canopy", 
                                                                  obs_montane_hypervolume_gaussian_1, 
                                                                  method="bootstrap", 
                                                                  n = 50, 
                                                                  points_per_resample =110,  # make sure you keep the sample size identical to what you want to compare with 
                                                                  cores = 6,
                                                                  verbose = TRUE, 
                                                                  weight_func = NULL)
hypervolumes_randomized_montane_match_canopy="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/montane_resample_match_canopy"

# Montane

#hypervolumes_randomized_montane= hypervolume_resample("montane_resample", 
                                                      obs_montane_hypervolume_gaussian_1, 
                                                     method="bootstrap", 
                                                     n = 1, 
                                                     points_per_resample ="sample_size",  
                                                     cores = 6,
                                                     verbose = TRUE, 
                                                     weight_func = NULL)

hypervolumes_randomized_montane="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/montane_resample"

#High montane #75

#hypervolumes_randomized_h_montane=hypervolume_resample("h_montane_resample", 
                                                                    obs_h_montane_hypervolume_gaussian_1, 
                                                                    method="bootstrap", 
                                                                    n = 1, 
                                                                    points_per_resample ="sample_size",  # make sure you keep the sample size identical to what you want to compare with 
                                                                    cores = 6,
                                                                    verbose = TRUE, 
                                                                    weight_func = NULL)

hypervolumes_randomized_h_montane="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/h_montane_resample"

hypervolumes_randomized_canopy_match_h_montane=hypervolume_resample("canopy_resample_match_h_montane", 
                                                       obs_canopy_hypervolume_gaussian_1, 
                                                       method="bootstrap", 
                                                       n = 50, 
                                                       points_per_resample =75,  # make sure you keep the sample size identical to what you want to compare with 
                                                       cores = 6,
                                                       verbose = TRUE, 
                                                       weight_func = NULL)

hypervolumes_randomized_canopy_match_h_montane="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/canopy_resample_match_h_montane"

hypervolumes_randomized_montane_match_h_montane=hypervolume_resample("montane_resample_match_h_montane", 
                                                                    obs_montane_hypervolume_gaussian_1, 
                                                                    method="bootstrap", 
                                                                    n = 50, 
                                                                    points_per_resample =75,  # make sure you keep the sample size identical to what you want to compare with 
                                                                    cores = 6,
                                                                    verbose = TRUE, 
                                                                    weight_func = NULL)


# Hypervolume STATS -------------------------------------------------------
## Summary stats

# Randomly estimated and observed
#For the hypervolumes that are going to be the limiting on sample size, we just need the observed one, ( example undestory_resmple, bamboo_resample)
#however we need to create a folder and a path for the calculation of the confidence intervals, so we will have both files as teh path ans as the file

#Obserevd social 
hypervolumes_observed_social_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/social_resample"
obs_social_hypervolume<-readRDS("/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/social_resample/resample 1.rds")
get_volume(obs_social_hypervolume)

hypervolumes_randomized_community_match_social_269="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/community_resample_match_social"
hvs_community_match_social_269 = to_hv_list (hypervolumes_randomized_community_match_social_269)
volumes_community_match_social_269<-get_volume(hvs_community_match_social_269 )
mean(volumes_community_match_social_269)
sd(volumes_community_match_social_269)

hypervolumes_randomized_nonsocial_match_social_269= "/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/nonsocial_resample_match_social"
hvs_nonsocial_match_social_269=to_hv_list(hypervolumes_randomized_nonsocial_match_social_269)
volumes_nonsocial_match_social_269<-get_volume(hvs_nonsocial_match_social_269)
mean(volumes_nonsocial_match_social_269)
sd(volumes_nonsocial_match_social_269)

#Runnning the code below now !!!!
statistics_social_community=hypervolume_overlap_confidence(hypervolumes_observed_social_269, hypervolumes_randomized_community_match_social_269,CI = .95, cores = 8) #269
saveRDS(statistics_social_community, "data/data_analyses/statistics_social_community.rds") 
statistics_social_community<-readRDS("data/data_analyses/statistics_social_community.rds")
social_community_mean<-mean (statistics_social_community$distribution [,1])
social_community_sd<-sd (statistics_social_community$distribution [,1])
statistics_social_community$jaccard
statistics_social_community$frac_unique_2

social_community_frac_unique_community<-mean (statistics_social_community$distribution [,4])
social_community_frac_unique_community_sd<-sd(statistics_social_community$distribution [,4])

#social_community_frac_unique_social<-mean (statistics_social_community$distribution [,3]) not needed because it is nested so shoudl not  have any unique fraction
#social_community_frac_unique_social_sd<-sd (statistics_social_community$distribution [,3]) not needed because it is nested so shoudl not  have any unique fraction

statistics_social_nonsocial=hypervolume_overlap_confidence(hypervolumes_observed_social_269, hypervolumes_randomized_nonsocial_match_social_269,CI = .95, cores = 8) #269
saveRDS(statistics_social_nonsocial, "data/data_analyses/statistics_social_nonsocial.rds") 
statistics_social_nonsocial<-readRDS("data/data_analyses/statistics_social_nonsocial.rds")
social_nonsocial_mean<-mean (statistics_social_nonsocial$distribution [,1]) # Jaccard similarity
social_nonsocial_sd<-sd (statistics_social_nonsocial$distribution [,1])
statistics_social_nonsocial$jaccard
statistics_social_nonsocial$frac_unique_1

social_nonsocial_frac_unique_social<-mean (statistics_social_nonsocial$distribution [,3])
social_nonsocial_frac_unique_social_sd<-sd (statistics_social_nonsocial$distribution [,3])

social_nonsocial_frac_unique_nonsocial<-mean (statistics_social_nonsocial$distribution [,4])
social_nonsocial_frac_unique_nonsocial_sd<-sd (statistics_social_nonsocial$distribution [,4])


##_###_##

hypervolumes_randomized_nonsocial="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/nonsocial_resample"
obs_nonsocial_hypervolume<-readRDS("/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/nonsocial_resample/resample 1.rds")
get_volume(obs_nonsocial_hypervolume)

# do not use this one use he one taht is commparable between the 3 groups 269 species 

#hypervolumes_randomized_community_match_nonsocial_408="/Users/jennymunoz/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/community_resample_match_non_social"
#hvs_community_match_nonsocial_408=to_hv_list(hypervolumes_randomized_community_match_nonsocial_408)
#volummes_community_match_nonsocial_408<-get_volume(hvs_community_match_nonsocial_408)
#mean(volummes_community_match_nonsocial_408)
#sd(volummes_community_match_nonsocial_408)

##Hypervolume by sample size for flock types
######Flocks:observed is for the flock that have few species so they are the limiting factor adn we match samples from the othe flocks to take into account sample size

#Observed understory
hypervolumes_observed_understory_61="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/understory_resample"
obs_understory_hypervolume<-readRDS("~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/understory_resample/resample 1.rds")
get_volume(obs_understory_hypervolume)

hypervolumes_randomized_canopy_match_understory_61="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/canopy_resample_match_understory"
hvs_canopy_match_understory_61 = to_hv_list (hypervolumes_randomized_canopy_match_understory_61)
volumes_canopy_match_understory_61<-get_volume(hvs_canopy_match_understory_61)
mean(volumes_canopy_match_understory_61)
sd(volumes_canopy_match_understory_61)

hypervolumes_randomized_montane_match_understory_61="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/montane_resample_match_understory"
hvs_montane_match_understory_61 = to_hv_list (hypervolumes_randomized_montane_match_understory_61)
volumes_montane_match_understory_61=get_volume(hvs_montane_match_understory_61)
mean(volumes_montane_match_understory_61)
sd(volumes_montane_match_understory_61)

hypervolumes_randomized_h_montane_match_understory_61="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/h_montane_resample_match_understory"
hvs_h_montane_match_understory_61= to_hv_list (hypervolumes_randomized_h_montane_match_understory_61)
volumnes_h_montane_match_understory_61=get_volume(hvs_h_montane_match_understory_61)
mean(volumnes_h_montane_match_understory_61)
sd(volumnes_h_montane_match_understory_61)

#Randomized bamboo
hypervolumes_observed_bamboo_43="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/bamboo_resample_overall"
hvs_bamboo_43=to_hv_list(hypervolumes_observed_bamboo_43)
volumes_bamboo_43=get_volume(hvs_bamboo_43)
mean(volumes_bamboo_43)
sd(volumes_bamboo_43)

#Observed bamboo
# using only one hypervolume
hypervolumes_observed_bamboo_43="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/bamboo_resample"
obs_bamboo_hypervolume<-readRDS("~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/bamboo_resample/resample 1.rds")
get_volume(obs_bamboo_hypervolume)

hypervolumes_randomized_understory_match_bamboo_43="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/understory_resample_match_bamboo"
hvs_understory_match_bamboo_43=to_hv_list(hypervolumes_randomized_understory_match_bamboo_43)
volumes_understory_match_bamboo_43=get_volume(hvs_understory_match_bamboo_43)
mean(volumes_understory_match_bamboo_43)
sd(volumes_understory_match_bamboo_43)


hypervolumes_randomized_canopy_match_bambo_43="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/canopy_resample_match_bamboo"
hvs_canopy_match_bambo_43=to_hv_list(hypervolumes_randomized_canopy_match_bambo_43)
volumes_canopy_match_bambo_43=get_volume(hvs_canopy_match_bambo_43)
mean(volumes_canopy_match_bambo_43)
sd(volumes_canopy_match_bambo_43)

hypervolumes_randomized_montane_match_bamboo_43="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/montane_resample_match_bamboo"
hvs_montane_match_bamboo_43=to_hv_list(hypervolumes_randomized_montane_match_bamboo_43)
volumes_montane_match_bamboo_43=get_volume(hvs_montane_match_bamboo_43)
mean(volumes_montane_match_bamboo_43)
sd(volumes_montane_match_bamboo_43)

hypervolumes_randomized_h_montane_match_bamboo_43="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/h_montane_resample_match_bamboo"
hvs_h_montane_match_bamboo_43=to_hv_list(hypervolumes_randomized_h_montane_match_bamboo_43)
volumes_h_montane_match_bamboo_43=get_volume(hvs_h_montane_match_bamboo_43)
mean(volumes_h_montane_match_bamboo_43)
sd(volumes_h_montane_match_bamboo_43)

#Observed canopy
hypervolumes_observed_canopy_110="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/canopy_resample"
obs_canopy_hypervolume<-readRDS("~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/canopy_resample/resample 1.rds")
get_volume(obs_canopy_hypervolume)

#Randomized canopy

hypervolumes_randomized_canopy_overall="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/canopy_resample_overall"
hvs_canopy_110=to_hv_list(hypervolumes_randomized_canopy_overall)
volumes_canopy_110=get_volume(hvs_canopy_110)
mean(volumes_canopy_110)
sd(volumes_canopy_110)


hypervolumes_randomized_montane_match_canopy_110="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/montane_resample_match_canopy"
hvs_montane_match_canopy_110=to_hv_list(hypervolumes_randomized_montane_match_canopy_110)
volumes_montane_match_canopy_110=get_volume(hvs_montane_match_canopy_110)
mean(volumes_montane_match_canopy_110)
sd(volumes_montane_match_canopy_110)

#Observed montane
hypervolumes_observed_montane_119="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/montane_resample"
obs_montane_hypervolume<-readRDS("~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/montane_resample/resample 1.rds")
get_volume(obs_montane_hypervolume)

#Observed h_montane
hypervolumes_observed_h_montane_75="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/h_montane_resample"
obs_h_montane_hypervolume<-readRDS("~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/h_montane_resample/resample 1.rds")
get_volume(obs_h_montane_hypervolume)


hypervolumes_randomized_canopy_match_h_montane_75="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/canopy_resample_match_h_montane"
hvs_canopy_match_h_montane_75=to_hv_list(hypervolumes_randomized_canopy_match_h_montane_75)
volumes_canopy_match_h_montane_75=get_volume(hvs_canopy_match_h_montane_75)
mean(volumes_canopy_match_h_montane_75)
sd(volumes_canopy_match_h_montane_75)

hypervolumes_randomized_montane_match_h_montane_75="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/montane_resample_match_h_montane"
hvs_montane_match_h_montane_75=to_hv_list(hypervolumes_randomized_montane_match_h_montane_75)
volumes_montane_match_h_montane_75=get_volume(hvs_montane_match_h_montane_75)
mean(volumes_montane_match_h_montane_75)
sd(volumes_montane_match_h_montane_75)

####
#Some key lines of code
hvs = to_hv_list (hypervolumes_randomized_community_match_nonsocial)
hvs@HVList # the list in one object
hvs@HVList[[30]]
volumes<-get_volume(hvs)
mean(volumes)
#get_centroid(hvs)

#Paths to files for analyses
#Observed understory
hypervolumes_observed_understory_61="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/understory_resample"
hypervolumes_randomized_canopy_match_understory_61="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/canopy_resample_match_understory"
hypervolumes_randomized_montane_match_understory_61="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/montane_resample_match_understory"
hypervolumes_randomized_h_montane_match_understory_61="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/h_montane_resample_match_understory"
#Observed bamboo
hypervolumes_observed_bamboo_43="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/bamboo_resample"
hypervolumes_randomized_understory_match_bamboo_43="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/understory_resample_match_bamboo"
hypervolumes_randomized_canopy_match_bambo_43="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/canopy_resample_match_bamboo"
hypervolumes_randomized_montane_match_bamboo_43="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/montane_resample_match_bamboo"
hypervolumes_randomized_h_montane_match_bamboo_43="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/h_montane_resample_match_bamboo"
#Observed canopy
hypervolumes_observed_canopy_110="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/canopy_resample"
hypervolumes_randomized_montane_match_canopy_110="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/montane_resample_match_canopy"
#Observed montane
hypervolumes_observed_montane_113="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/montane_resample"
#Observed h_montane
hypervolumes_observed_h_montane_75="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/h_montane_resample"
hypervolumes_randomized_canopy_match_h_montane_75="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/canopy_resample_match_h_montane"
hypervolumes_randomized_montane_match_h_montane_75="~/Desktop/PhD/Phd by chapters/Chapter 1-Community context/Chapter1_flocks_in_community_context/Objects/montane_resample_match_h_montane"


# Hypervolume stats overlap -----------------------------------------------

##STATISTICS###_##_##_##_##_##_##_##_##_##_##_##_##_##_##_##_##_##_##_##_##_##_##_##
# Statistics on hypevolume overlpa using the randomized samples (controling for sample size)
statistics_understory_canopy=hypervolume_overlap_confidence(hypervolumes_observed_understory_61, hypervolumes_randomized_canopy_match_understory_61,CI = .95, cores = 10) #43
saveRDS(statistics_understory_canopy, "data/data_analyses/statistics_understory_canopy.rds") 
statistics_understory_canopy<-readRDS("data/data_analyses/statistics_understory_canopy.rds")
understory_canopy_mean<-mean (statistics_understory_canopy$distribution [,1])
understory_canopy_sd<-sd (statistics_understory_canopy$distribution [,1])
statistics_understory_canopy$jaccard

statistics_understory_montane=hypervolume_overlap_confidence(hypervolumes_observed_understory_61,hypervolumes_randomized_montane_match_understory_61,CI = .95, cores = 10) # 43
saveRDS(statistics_understory_montane, "data/data_analyses/statistics_understory_montane.rds") 
statistics_understory_montane<-readRDS("data/data_analyses/statistics_understory_montane.rds")
understory_montane_mean<-mean (statistics_understory_montane$distribution [,1])
understory_montane_sd<-sd (statistics_understory_montane$distribution [,1])

statistics_understory_h_montane=hypervolume_overlap_confidence(hypervolumes_observed_understory_61, hypervolumes_randomized_h_montane_match_understory_61,CI = .95, cores = 8) #43
saveRDS(statistics_understory_h_montane, "data/data_analyses/statistics_understory_h_montane.rds") 
statistics_understory_h_montane= readRDS("data/data_analyses/statistics_understory_h_montane.rds")
mean (statistics_understory_h_montane$distribution [,1])
sd (statistics_understory_h_montane$distribution [,1])
statistics_understory_h_montane$jaccard

statistics_bamboo_understory=hypervolume_overlap_confidence(hypervolumes_observed_bamboo_43, hypervolumes_randomized_understory_match_bamboo_43,CI = .95, cores = 20) # sample size 43
saveRDS(statistics_bamboo_understory,"data/data_analyses/statistics_bamboo_understory.rds")
statistics_bamboo_understory=readRDS("data/data_analyses/statistics_bamboo_understory.rds")
mean(statistics_bamboo_understory$distribution [,1])
sd(statistics_bamboo_understory$distribution [,1])

statistics_bamboo_canopy=hypervolume_overlap_confidence(hypervolumes_observed_bamboo_43,hypervolumes_randomized_canopy_match_bambo_43,CI = .95, cores = 6) # sample size 43
saveRDS(statistics_bamboo_canopy,"data/data_analyses/statistics_bamboo_canopy.rds")
statistics_bamboo_canopy<-readRDS("data/data_analyses/statistics_bamboo_canopy.rds")
mean(statistics_bamboo_canopy$distribution [,1])
sd(statistics_bamboo_canopy$distribution [,1])


statistics_bamboo_montane=hypervolume_overlap_confidence(hypervolumes_observed_bamboo_43, hypervolumes_randomized_montane_match_bamboo_43,CI = .95, cores = 6) # sample size 43
saveRDS(statistics_bamboo_montane,"data/data_analyses/statistics_bamboo_montane.rds")
statistics_bamboo_montane<-readRDS("data/data_analyses/statistics_bamboo_montane.rds")
mean(statistics_bamboo_montane$distribution [,1])
sd(statistics_bamboo_montane$distribution [,1])

statistics_bamboo_h_montane=hypervolume_overlap_confidence(hypervolumes_observed_bamboo_43, hypervolumes_randomized_h_montane_match_bamboo_43,CI = .95, cores = 6) # sample size 43
saveRDS(statistics_bamboo_h_montane,"data/data_analyses/statistics_bamboo_h_montane.rds")
statistics_bamboo_h_montane<-readRDS("data/data_analyses/statistics_bamboo_h_montane.rds")
mean(statistics_bamboo_h_montane$distribution [,1])
sd(statistics_bamboo_h_montane$distribution [,1])

statistics_canopy_montane=hypervolume_overlap_confidence(hypervolumes_observed_canopy_110, hypervolumes_randomized_montane_match_canopy_110,CI = .95, cores = 6) # sample size 110
saveRDS(statistics_canopy_montane,"data/data_analyses/statistics_canopy_montane.rds")
statistics_canopy_montane<-readRDS("data/data_analyses/statistics_canopy_montane.rds")
mean(statistics_canopy_montane$distribution [,1])
sd(statistics_canopy_montane$distribution [,1])

statistics_h_montane_canopy =hypervolume_overlap_confidence(hypervolumes_observed_h_montane_75,hypervolumes_randomized_canopy_match_h_montane_75,CI = .95, cores = 6) #sample size 75
saveRDS(statistics_h_montane_canopy,"data/data_analyses/statistics_h_montane_canopy.rds")
statistics_h_montane_canopy<-readRDS("data/data_analyses/statistics_h_montane_canopy.rds")
mean(statistics_h_montane_canopy$distribution [,1])
sd(statistics_h_montane_canopy$distribution [,1])

statistics_h_montane_montane =hypervolume_overlap_confidence(hypervolumes_observed_h_montane_75,hypervolumes_randomized_montane_match_h_montane_75,CI = .95, cores = 6) # sample size 75
saveRDS(statistics_h_montane_montane,"data/data_analyses/statistics_h_montane_montane.rds")
statistics_h_montane_montane<-readRDS("data/data_analyses/statistics_h_montane_montane.rds")
mean(statistics_h_montane_montane$distribution [,1])
sd(statistics_h_montane_montane$distribution [,1])



## End(Not run)
