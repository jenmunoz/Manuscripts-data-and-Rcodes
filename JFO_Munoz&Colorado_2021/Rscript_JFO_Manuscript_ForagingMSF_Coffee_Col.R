###############################JENNY MUNOZ#####################################
###############################################################################
## Proyect Foraging ecology of residents and migrants
## R-code for manuscript Generalized linear models
## Jenny Munoz
#### last update: February 20 2017
################################################################################

#Assumptions": 
#check the R version
sessionInfo()
# Install: 
install.packages("car")
install.packages("vegan")
install.packages("AICcmodavg")
install.packages("gridExtra ")
install.packages("tidyverse")
# Loading packages --------------------------------------------------------
##general packages  and visualization
library(lattice)
library(ggplot2)
library(plotly)
library(car)
library(visreg)
library(dplyr)
#library(warbleR)
library(vegan)
library(permute)
library(sjPlot) #Make tables of the glm stimates https://rdrr.io/cran/sjPlot/man/sjt.glm.html
#Model selection
install.packages("AICcmodavg")
library(AICcmodavg)

library(tidyverse)
library(gridExtra)

#########################################################################################################

# Reading and cleaning the data -------------------------------------------

# Read the data from the file, View the first few lines
#foraging<-read.csv(file.choose("foraging_all_ species_coffee_2011"), stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
#head(foraging)
#summary(foraging)
#str(foraging)
#which(is.na(foraging)==T)
#foraging

#################### REANALIZING THE DATA FOR MANUSCRIPT ----------------------------------------------------
##############################################################################################
###Need to include duration of the observation as an offset because glm-poisoon only work with discrete data 

###################################
#Reading the data
###################################

foraging<-read.csv(file.choose("foraging_all_species_offset"), stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))  

head(foraging)
summary(foraging)
str(foraging)
which(is.na(foraging)==T)
foraging

#Converting variables as needed
foraging$sociality<- as.factor(foraging$sociality) 
foraging$movement<-as.numeric(foraging$movement)
foraging$capture<-as.numeric(foraging$capture)
foraging$capture<-as.integer(foraging$capture)
foraging$flocksizespecies<- as.numeric(foraging$flocksizespecies)
foraging$flocksizeind<- as.numeric(foraging$flocksizeind)
foraging$dayofseason<- as.integer(foraging$dayofseason)

###################################
# Filtering the data
###################################

library(dplyr)

fusca<-filter(foraging, species=="Sethophaga fusca")
canadensis<-filter(foraging, species=="Cardelina canadensis")
cerulea<-filter(foraging,species=="Setophaga cerulea")
peregrina<-filter(foraging,species=="Oreothlypis peregrina")
guira<-filter(foraging,species=="Hemithraupis guira")
chrysops<-filter(foraging,species=="Zimmerius chrysops")
pitiayumi<-filter(foraging,species=="Parula pitiayumi")
canadensis<-filter(foraging, species=="Cardelina canadensis")


###################################
#####################Checking model assumptions and data distributions
###################################

##plot a frequency data distribution 
###R##it seems as a poisson 
hist(fusca$foragingrate)
hist(cerulea$foragingrate)
hist(canadensis$foragingrate)
hist(peregrina$foragingrate)
hist(chrysops$foragingrate)
hist(guira$foragingrate)
hist(pitiayumi$foragingrate)


#checking the normality of the data ### seems like it is non-normal
qqnorm(fusca$foragingrate)
qqline(fusca$foragingrate, lty=2)
qqnorm(cerulea$foragingrate)
qqnorm(canadensis$foragingrate)
qqnorm(peregrina$foragingrate)
qqnorm(chrysops$foragingrate)
qqnorm(guira$foragingrate)
qqnorm(pitiayumi$foragingrate)

shapiro.test(flocksd$Number_of_species) # iF P<0.05 daa deviate from normality, all deviate from normality

shapiro.test(fusca$foragingrate)
shapiro.test(cerulea$foragingrate)
shapiro.test(canadensis$foragingrate)
shapiro.test(peregrina$foragingrate)
shapiro.test(chrysops$foragingrate)
shapiro.test(guira$foragingrate)
shapiro.test(pitiayumi$foragingrate)

# test for overdisspersion 
# is the variance is higher that the mean the data is overdispersed
tapply(fusca$foragingrate, fusca$sociality, mean)
tapply(fusca$foragingrate, fusca$sociality, var)

tapply(canadensis$foragingrate, canadensis$sociality, mean)
tapply(canadensis$foragingrate, canadensis$sociality, var)

tapply(cerulea$foragingrate, cerulea$sociality, mean)
tapply(cerulea$foragingrate, cerulea$sociality, var)

tapply(peregrina$foragingrate, peregrina$sociality, mean)
tapply(peregrina$foragingrate, peregrina$sociality, var)

# in chrysops is not overdisspersed, then poisson distribution
tapply(chrysops$foragingrate, chrysops$sociality, mean)
tapply(chrysops$foragingrate, chrysops$sociality, var)

tapply(guira$foragingrate, guira$sociality, mean)
tapply(guira$foragingrate, guira$sociality, var)

tapply(pitiayumi$foragingrate, pitiayumi$sociality, mean)
tapply(pitiayumi$foragingrate, pitiayumi$sociality, var)
# in chrysops is the only speces which data is not overdisspersed, then poisson distribution, for all the other I will apply quasipoisson

##because they deviate from normality, I will run GLM instead of lm
# Data is overdispersed then I will use quasi poisson distribution

###################################################################################################
#1. Mixed species flocks features
###################################################################################################

#Filter only the flock data
flocks<-filter(foraging,sociality=="Flock")
flockfeatures<-na.omit(flocks)                             #omit NA values
summarise_each(flockfeatures,funs(mean))                   #Summarise all the columns
summarise(flockfeatures,mean(flocksizespecies))            #Summarise flock size 
summarise(flockfeatures,mean(flocksizeind))                #mean
summarise(flockfeatures,sd(flocksizespecies))              #standard deviation
summarise(flockfeatures,sd(flocksizeind)) 

summarise(flockfeatures,max(flocksizespecies))             #max and mins
summarise(flockfeatures,min(flocksizespecies))
summarise(flockfeatures,max(flocksizeind))
summarise(flockfeatures,min(flocksizeind))

flocksfusca<-filter(fusca,sociality=="Flock")
a<-summarise(flocksfusca,mean(foragingrate))


#Normality of the variables
hist(flockfeatures$flocksizespecies)
hist(flockfeatures$flocksizeind)
qqnorm(flockfeatures$flocksizeind) #normality
qqline(flockfeatures$flocksizeind, lty=2)
qqnorm(flockfeatures$flocksizespecies) #normality
qqline(flockfeatures$flocksizespecies, lty=2)

#Correlation test
#The alternative hypothesis of interest is that the flocksize is positively associated with the flockind.
cor.test(flockfeatures$flocksizespecies, flockfeatures$flocksize)

###################################
#2## Foraging rate vs sociality -----------------------------------------------------------
###################################

##### Model selection 
glm1<-glm(capture~sociality+dayofseason, data=fusca,family=poisson(link="log"), offset=log(minutes))

glm1<-glm(capture~sociality+dayofseason, data=canadensis,family=poisson(link="log"), offset=log(minutes))
glm1<-glm(capture~sociality+dayofseason, data=cerulea,family=poisson(link="log"), offset=log(minutes)) #Overdisspersion to high model innacurat, small sample size
glm1<-glm(capture~sociality+dayofseason, data=peregrina,family=poisson(link="log"), offset=log(minutes))
glm1<-glm(capture~sociality+dayofseason, data=pitiayumi,family=poisson(link="log"), offset=log(minutes))
glm1<-glm(capture~sociality+dayofseason, data=guira,family=poisson(link="log"), offset=log(minutes))
glm1<-glm(capture~sociality+dayofseason, data=chrysops1,family=poisson(link="log"), offset=log(minutes))

glm2<-update(glm1, . ~ . - dayofseason) #social context
glm3<-update(glm1, . ~ . - sociality) #day of season
glm4<-update(glm3, . ~ . - dayofseason)
glmQ1<-update(glm1, family=quasipoisson(link="log"))
glmQ2<-update(glm2, family=quasipoisson(link="log"))
glmQ3<-update(glm3, family=quasipoisson(link="log"))
glmQ4<-update(glm4, family=quasipoisson(link="log"))

glm5<-glm(foragingrate~1, data =fusca, family=poisson(link="log")) #Null model equivalent to glm2 long version

# Summary of the model allow me to interpret the estimates of the parameters (e.g effect sizes) in the model and the difference from cero and between them.
# The summary give me the effect size. 
summary(glm1)
summary(glmQ1)

#####overdisperssion parameters 
dfun<-function(glmQ1q){with(glmQ1,sum((weights * residuals^2)[weights > 0])/df.residual)}
dfun(glmQ1)
summary(glmQ1)

#fusca=2.20
#cerulea=5.34
#peregrina=4.41~4
#canadensis=1.76
#pitiayumi=2.29
#guira=3.37
#chrysops=1.06

#### Model selection using  AICcmodavg
library(AICcmodavg)

###Model used  (glm1Q,glm2Q,glm3Q,glm4Q), for all the species which have overdispersed datta, remember to use the glm models instead of glmQ
aictab(list(glm1,glm2,glm3,glm4),
       modnames=c("Socialcontext+Dayofseason",
                  "Socialcontext",
                  "Dayofseason",
                  "Intercept"),
       c.hat=3.37)

summary (glmQ2)

#### To make the table for AIC cfor species that are not overdisperssed

aictab(list(glm1,glm2,glm3,glm4),
       modnames=c("Socialcontext+Dayofseason",
                  "Socialcontext",
                  "Dayofseason",
                  "Intercept"))

######SUMMARY TABLE
#A beautiful table in html for the parameters of glm models
#For the parameter stimates we need to consider the models with the Quasipoisson distribution, the points estimates are identical
# to the model with poisson distribution but the standard error and confidence intervals are wider!
sjt.glm(glmQ1,glmQ2,glmQ3,
        depvar.labels = c("Model1: Socialcontext+Dayofseason","Model2:Socialcontext","Model3:Day of season"),
        pred.labels = NULL,
        show.aicc=TRUE, 
        show.family=TRUE, 
        group.pred = FALSE,
        exp.coef = FALSE,  # if true Regression and cof.intervals are exponentiaded st.error are not in the unstransformed scale
        p.numeric=TRUE,
        robust = TRUE,
        show.se = TRUE,
        show.r2=TRUE,
        show.dev=TRUE,
        show.chi2=TRUE,
        cell.spacing = 0.001,
        sep.column = FALSE)
###################################################################################
#Anova test for differences between groups
##################################################################################

## Anova allow me to do hyphothesis testing
# The ANOVA atable compares the fit of two models, the null model foraging~1 vs foraging~sociality+....
## f test appropiate for quasipoisson distributions, and type3 anova order orthe term does not matter
# test the null H that there is no differences in the foraging rate  among flocking and non-flocking individuals
#Note: Use anova(test = "F"), ie an F test, if testing hypotheses that use gaussian, quasibinomial or quasipoisson link functions.
#This is because the quasi-likelihood is not a real likelihood and the generalized log-likelihood ratio test is not accurate FROM dolph class
#Anova type=3 rom the car package, in that case order of the terms in the model does not matter and neither does hirarchy

#anova(model, test="Chi") using the "best model" I guess i can use any of the two qasi or poisson since the stimates are the same?
### For the migrants 
glm1<-glm(capture~sociality+dayofseason, data=fusca1,family=poisson(link="log"), offset=log(minutes))
glmQ1<-update(glm1, family=quasipoisson(link="log"))
Anova(glmQ1, type = 3, test="F")

glm1<-glm(capture~sociality+dayofseason, data=peregrina,family=poisson(link="log"), offset=log(minutes))
summary(glm1)
summary(glmQ1)
glmQ1<-update(glm1, family=quasipoisson(link="log"))
Anova(glmQ1, type = 3, test="F")

glm1<-glm(capture~sociality+dayofseason, data=canadensis,family=poisson(link="log"), offset=log(minutes))
glmQ1<-update(glm1, family=quasipoisson(link="log"))
glm2<-update(glm1, . ~ . - dayofseason) #social context
glmQ2<-update(glm2, family=quasipoisson(link="log"))

Anova(glmQ2, type = 3, test="F")

glm1<-glm(capture~sociality+dayofseason, data=cerulea,family=poisson(link="log"), offset=log(minutes))
glmQ1<-update(glm1, family=quasipoisson(link="log"))
Anova(glmQ1, type = 3, test="F")

#For the residents
glm1<-glm(capture~sociality+dayofseason, data=chrysops,family=poisson(link="log"), offset=log(minutes))
glm2<-update(glm1, . ~ . - dayofseason) #social context
Anova(glm2, type = 3, test="F")
summary(glm2)

glm1<-glm(capture~sociality+dayofseason, data=pitiayumi,family=poisson(link="log"), offset=log(minutes))
glmQ1<-update(glm1, family=quasipoisson(link="log"))
glm2<-update(glm1, . ~ . - dayofseason) #social context
glmQ2<-update(glm2, family=quasipoisson(link="log"))
Anova(glmQ2, type = 3, test="F")


glm1<-glm(capture~sociality+dayofseason, data=guira,family=poisson(link="log"), offset=log(minutes))
glmQ1<-update(glm1, family=quasipoisson(link="log"))
glm2<-update(glm1, . ~ . - dayofseason) #social context
glmQ2<-update(glm2, family=quasipoisson(link="log"))
Anova(glmQ2, type = 3, test="F")

#####Movement rate 

##############################################################################
# Part 2 Movement rate -----------------------------------------------------------

glm1m<-glm(movement~sociality+dayofseason, data=fusca,family=poisson(link="log"), offset=log(minutes))
glm1m<-glm(movement~sociality+dayofseason, data=canadensis,family=poisson(link="log"), offset=log(minutes))
glm1m<-glm(movement~sociality+dayofseason, data=cerulea,family=poisson(link="log"), offset=log(minutes))
glm1m<-glm(movement~sociality+dayofseason, data=peregrina,family=poisson(link="log"), offset=log(minutes))
glm1m<-glm(movement~sociality+dayofseason, data=pitiayumi,family=poisson(link="log"), offset=log(minutes))
glm1m<-glm(movement~sociality+dayofseason, data=guira,family=poisson(link="log"), offset=log(minutes))
glm1m<-glm(movement~sociality+dayofseason, data=chrysops,family=poisson(link="log"), offset=log(minutes))

glm2m<-update(glm1m, . ~ . - dayofseason) #social context
glm3m<-update(glm1m, . ~ . - sociality) #day of season
glm4m<-update(glm3m, . ~ . - dayofseason)
glmQ1m<-update(glm1m, family=quasipoisson(link="log"))
glmQ2m<-update(glm2m, family=quasipoisson(link="log"))
glmQ3m<-update(glm3m, family=quasipoisson(link="log"))
glmQ4m<-update(glm4m, family=quasipoisson(link="log"))

# Summary of the model allow me to interpret the estimates of the parameters (e.g effect sizes) in the model and the difference from cero and between them.
# The summary give me the effect size. 
summary(glm1m)
summary(glmQ1m)

#####overdisperssion parameters 
summary(glmQ1m)

#fusca=2.26
#cerulea=5.05
#peregrina=3.19
#canadensis= 4
#pitiayumi=3.45
#guira=
#chrysops=4

###Model used  (glm1Q,glm2Q,glm3Q,glm4Q), for all the species which have overdispersed datta, remember to use the glm models instead of glmQ
aictab(list(glm1m,glm2m,glm3m,glm4m),
       modnames=c("Socialcontext+Dayofseason",
                  "Socialcontext",
                  "Dayofseason",
                  "Intercept"),
       c.hat=3.45)

summary (glmQ2)

######################
######SUMMARY TABLE
#A beautiful table in html for the parameters of glm models
#For the parameter stimates we need to consider the models with the Quasipoisson distribution, the points estimates are identical
# to the model with poisson distribution but the standard error and confidence intervals are wider!
sjt.glm(glmQ1m,glmQ2m,glmQ3m,
        depvar.labels = c("Model1: Socialcontext+Dayofseason","Model2:Socialcontext","Model3:Day of season"),
        pred.labels = NULL,
        show.aicc=TRUE, 
        show.family=TRUE, 
        group.pred = FALSE,
        exp.coef =TRUE,  # if true Regression and cof.intervals are exponentiaded st.error are not in the unstransformed scale
        p.numeric=TRUE,
        robust = TRUE,
        show.se = TRUE,
        show.r2=TRUE,
        show.dev=TRUE,
        show.chi2=TRUE,
        cell.spacing = 0.001,
        sep.column = FALSE)

#########Anova test###########################################################
## Anova allow me to do hyphothesis testing, note that I used the best selected model for each species 
# The ANOVA atable compares the fit of two models, the null model foraging~1 vs foraging~sociality+....
## f test appropiate for quasipoisson distributions, and type3 anova order orthe term does not matter
# test the null H that there is no differences in the foraging rate  among flocking and non-flocking individuals
#Note: Use anova(test = "F"), ie an F test, if testing hypotheses that use gaussian, quasibinomial or quasipoisson link functions.
#This is because the quasi-likelihood is not a real likelihood and the generalized log-likelihood ratio test is not accurate FROM dolph class
#Anova type=3 rom the car package, in that case order of the terms in the model does not matter and neither does hirarchy

#anova(model, test="Chi") using the "best model" I guess i can use any of the two qasi or poisson since the stimates are the same?

### For the migrants 
glm1m<-glm(movement~sociality+dayofseason, data=fusca,family=poisson(link="log"), offset=log(minutes))
glmQ1m<-update(glm1m, family=quasipoisson(link="log"))
glm2m<-update(glm1m, . ~ . - dayofseason) #social context
glmQ2m<-update(glm2m, family=quasipoisson(link="log"))
Anova(glmQ2m, type = 3, test="F")

glm1m<-glm(movement~sociality+dayofseason, data=peregrina1,family=poisson(link="log"), offset=log(minutes))
glmQ1m<-update(glm1m, family=quasipoisson(link="log"))
glm2m<-update(glm1m, . ~ . - dayofseason) #social context
glmQ2m<-update(glm2m, family=quasipoisson(link="log"))
Anova(glmQ2m, type = 3, test="F")

glm1m<-glm(movement~sociality+dayofseason, data=canadensis1,family=poisson(link="log"), offset=log(minutes))
glmQ1m<-update(glm1m, family=quasipoisson(link="log"))
glm2m<-update(glm1m, . ~ . - dayofseason) #social context
glmQ2m<-update(glm2m, family=quasipoisson(link="log"))
Anova(glmQ2m, type = 3, test="F")

#For the residents
glm1m<-glm(movement~sociality+dayofseason, data=pitiayumi1,family=poisson(link="log"), offset=log(minutes))
glmQ1m<-update(glm1m, family=quasipoisson(link="log"))
Anova(glmQ1m, type = 3, test="F") 
summary(glmQ1m)

glm1m<-glm(movement~sociality+dayofseason, data=chrysops1,family=poisson(link="log"), offset=log(minutes))
glmQ1m<-update(glm1m, family=quasipoisson(link="log"))
glm2m<-update(glm1m, . ~ . - dayofseason) #social context
glmQ2m<-update(glm2m, family=quasipoisson(link="log"))
Anova(glmQ2m, type = 3, test="F")
summary (glmQ2m)

###################################################################################################
#3 The effect of group size -----------------------------------------------
##################################################################################################
#### 
library(dplyr)
flocks<-filter(foraging,sociality=="Flock")
flocks

#Model selection for flock size
glm1fs<-glm(capture~flocksizespecies+I(flocksizespecies^2), data=fusca,family=poisson(link="log"), offset=log(minutes))
glm1fs<-glm(capture~flocksizespecies+I(flocksizespecies^2), data=canadensis,family=poisson(link="log"), offset=log(minutes))
glm1fs<-glm(capture~flocksizespecies+I(flocksizespecies^2), data=cerulea,family=poisson(link="log"), offset=log(minutes)) #Overdisspersion to high model innacurat, small sample size
glm1fs<-glm(capture~flocksizespecies+I(flocksizespecies^2), data=peregrina,family=poisson(link="log"), offset=log(minutes))
glm1fs<-glm(capture~flocksizespecies+I(flocksizespecies^2), data=pitiayumi,family=poisson(link="log"), offset=log(minutes))
glm1fs<-glm(capture~flocksizespecies+I(flocksizespecies^2), data=guira,family=poisson(link="log"), offset=log(minutes))
glm1fs<-glm(capture~flocksizespecies+I(flocksizespecies^2), data=chrysops,family=poisson(link="log"), offset=log(minutes))

glm2fs<-glm(capture~flocksizespecies, data=chrysops,family=poisson(link="log"), offset=log(minutes))
glm3fs<-update(glm2fs, . ~ . - flocksizespecies)
glm4fs<-glm(capture~I(flocksizespecies^2), data=chrysops,family=poisson(link="log"), offset=log(minutes))
glmQ1fs<-update(glm1fs, family=quasipoisson(link="log"))
glmQ2fs<-update(glm2fs, family=quasipoisson(link="log"))
glmQ3fs<-update(glm3fs, family=quasipoisson(link="log"))
glmQ4fs<-update(glm4fs, family=quasipoisson(link="log"))

# Summary of the model allow me to interpret the estimates of the parameters (e.g effect sizes) in the model and the difference from cero and between them.
# The summary give me the effect size. 
summary(glm1fs)
summary(glmQ1fs)

#####overdisperssion parameters 
dfun<-function(glmQ1fs){with(glmQ1fs,sum((weights * residuals^2)[weights > 0])/df.residual)}
dfun(glmQ1fs)
summary(glmQ4fs)

#fusca=2.50
#cerulea=4
#peregrina=4
#canadensis=1.76
#pitiayumi=2.75
#guira=3.17
#chrysops=

#### Model selection using  AICcmodavg
library(AICcmodavg)

###Model used  (glm1,glm2,glm3,glm4), for all the species which have overdispersed datta, remember to use the glm models instead of glmQ
aictab(list(glm1fs,glm2fs,glm3fs,glm4fs),
       modnames=c("FlockSizeSpecies+POLY",
                  "FlockSizeSpecies",
                  "Intercept",
                  "poly"),
       c.hat=2.50)

aictab(list(glm1fs,glm2fs,glm3fs),
       modnames=c("FlockSizeSpecies+POLY",
                  "FlockSizeSpecies",
                  "Intercept"),
       c.hat=1.05)


#### To make the table for AIC cfor species that are not overdisperssed

aictab (list(glm1,glm2fs,glm3fs,glm4fs,glm5fs),
        modnames=c("FlockSizeSpecies+Dayofseason",
                   "FlockSizeSpecies",
                   "Dayofseason",
                   "Intercept",
                   "poly"))

######SUMMARY TABLE
#A beautiful table in html for the parameters of glm models
#For the parameter stimates we need to consider the models with the Quasipoisson distribution, the points estimates are identical
# to the model with poisson distribution but the standard error and confidence intervals are wider!
sjt.glm(glmQ1fs,glmQ2fs,glmQ4fs,
        depvar.labels = c("Model1: Socialcontext+Dayofseason","Model2:Socialcontext","Model3:Day of season"),
        pred.labels = NULL,
        show.aicc=TRUE, 
        show.family=TRUE, 
        group.pred = FALSE,
        exp.coef = FALSE,  # if true Regression and cof.intervals are exponentiaded st.error are not in the unstransformed scale
        p.numeric=TRUE,
        robust = TRUE,
        show.se = TRUE,
        show.r2=TRUE,
        show.dev=TRUE,
        show.chi2=TRUE,
        cell.spacing = 0.001,
        sep.column = FALSE)

########################################################################
# Graphs for manuscript ---------------------------------------------------
########################################################################
########################################################################
# Read the data from the file, View the first few lines
foraging<-read.csv(file.choose("foraging_all_ species_offset"), stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
str(foraging)

####1. Plot foraging rate vs social context
foraging$plotMigrants<- as.factor(foraging$plotMigrants)
foraging$plotResidents<- as.factor(foraging$plotResidents)

#For the migrants, useDingbats=FALSE allow me to import this ch to illustrator without problems
pdf("Fig1Socialcontexmigrants.pdf",width=14,height=8,useDingbats=FALSE)
stripchart(foraging$foragingRate~foraging$plotMigrants, method='jitter', jitter=0.3, cex.axis=1.3,vertical=TRUE, pch=19, cex=1.5, ylim=c(0,25), col=c("black","gray"),  ylab='Foraging rate ( #captures / min)', xlab='Social context')
points(c(1,2,3,4,5,6,7,8), tapply (foraging$foragingRate,foraging$plotMigrants,mean), pch="----", col=c("black","darkgray"),cex=2  )
dev.off()


# For the residents
pdf("Fig1SocialcontexResidents.pdf",width=14,height=8, useDingbats=FALSE)
stripchart(foraging$foragingRate~foraging$plotResidents, method='jitter', jitter=0.3, cex.axis=1.3,vertical=TRUE, pch=19, cex=1.5,cex.lab=1.3, ylim=c(0,15), col=c("black","gray"),  ylab='Foraging rate (#captures / min)', xlab='Social context')
points(c(1,2,3,4,5,6), tapply (foraging$foragingRate,foraging$plotResidents,mean), pch="----", col=c("black","darkgray"),cex=2  )
dev.off()

# If I want to include manually the parameter estimates from the model
points(c(1),c(3.66), col="blue", pch="-",cex=4)
# if I want to include confidence intervals
lines( c(1,1), c(3,10)) 
lines( c(2,2), c(5,10))
lines ()
coefficients(glm1)
#add the means##

##########################################################
####2. Plot movement rate vs social context
foraging$plotMigrants<- as.factor(foraging$plotMigrants)
foraging$plotResidents<- as.factor(foraging$plotResidents)
#For the migrants
foraging
pdf("Fig2SocialcontexmigrantsMVT.pdf",width=14,height=8,useDingbats=FALSE)
stripchart(foraging$movementRate~foraging$plotMigrants, method='jitter', jitter=0.3, cex.axis=1.3,vertical=TRUE, pch=19, cex=1.5, ylim=c(0,55), col=c("black","gray"),  ylab='Foraging rate (#captures / min)', xlab='Social context')
points(c(1,2,3,4,5,6,7,8), tapply (foraging$movementRate,foraging$plotMigrants,mean), pch="----", col=c("black","grey"),cex=2  )
dev.off()


# For the residents
pdf("Fig2bSocialcontexResidentsMvt.pdf",width=14,height=8,useDingbats=FALSE)
stripchart(foraging$movementRate~foraging$plotResidents, method='jitter', jitter=0.3, cex.axis=1.3,vertical=TRUE, pch=19, cex=1.5, ylim=c(0,50), col=c("black","gray"),  ylab='Movement rate (#movements / min)', xlab='Social context')
points(c(1,2,3,4,5,6), tapply (foraging$movementRate,foraging$plotResidents,mean), pch="----", col=c("black","darkgray"),cex=2)
dev.off()

##############################################################################################
# Foraging maneuvers
###################################################################################################
##############################################################################################
#reading the data
maneuvers<-read.csv(file.choose("maneuver"),stringsAsFactors=FALSE,strip.white=TRUE,na.strings=c("NA",""))

#filtering data
#fuscamaneuver<-filter(maneuvers,species=="Setophaga fusca")
#parulamaneuver<-filter(maneuvers,species=="Parula pitiayumi")

#variables as factors
maneuvers$context<-as.factor(maneuvers$context)

### Barchart for species, I can vary the number of columns and rows, and by using the driver command with dimensions, I can make the graph to fit my desire lenght###
#barchart(proportion~maneuver|species,data=maneuvers, layout=c(3,2),scale=list(x=list(rot=0,cex=0.7),y=list(cex=0.7)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,1),box.ratio=1.5,col="gray")
# using context as a grouping category
#barchart(proportion~maneuver|species, groups=context, data=maneuvers, layout=c(3,2),scale=list(x=list(rot=0,cex=0.7),y=list(cex=0.7)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,1),box.ratio=1.5, auto.key=TRUE)
# in gay colors
barchart(proportion~maneuver|species, groups=context, data=maneuvers, layout=c(3,2),scale=list(x=list(rot=0,cex=0.7),y=list(cex=0.7)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,1),box.ratio=1.5, col=c("grey","white"))
barchart(percentage~maneuver|species, groups=context, data=maneuvers, layout=c(3,2),scale=list(x=list(rot=0,cex=0.7),y=list(cex=0.7)),ylab="Percentage",xlab="Foraging maneuver",ylim=c(0,100),box.ratio=1.5, col=c("grey","white"))

pdf(file="Fig4ManeuversResidentsMigrants_Proportion.pdf",width=14,height=10)
barchart(proportion~maneuver|species, groups=context, data=maneuvers, layout=c(3,2),scale=list(x=list(rot=0,cex=1),y=list(cex=1)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,1),box.ratio=1.5, col=c("grey","white"))
dev.off()

pdf(file="Fig4ManeuversResidentsMigrants_Percentage.pdf",width=14,height=10)
barchart(percentage~maneuver|species, groups=context, data=maneuvers, layout=c(3,2),scale=list(x=list(rot=0,cex=1),y=list(cex=1)),ylab="Percentage",xlab="Foraging maneuver",ylim=c(0,100),box.ratio=1.5, col=c("grey","white"))
dev.off()

tiff("Plotmaneuvers.tiff", res = 300, width =1000 , height =1015 )
barchart(proportion~maneuver|species, groups=context, data=maneuvers, layout=c(3,2),scale=list(x=list(rot=0,cex=1),y=list(cex=1)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,1),box.ratio=1.5, col=c("grey","white"))
dev.off()

##############################################################################################
################Substrate########################################################################
##################################################################################################################################

library(lattice)

#reading the data
substrate<-read.csv(file.choose("substrate"),stringsAsFactors=FALSE,strip.white=TRUE,na.strings=c("NA",""))

#variables as factors
substrate$context<-as.factor(substrate$context)
str(substrate)

# in gay colors
barchart(proportion~substrate|species, groups=context, data=substrate, layout=c(3,2),scale=list(x=list(rot=0,cex=0.7),y=list(cex=0.7)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,1),box.ratio=1.5, auto.key=TRUE)
barchart(proportion~substrate|species, groups=context, data=substrate, layout=c(3,2),scale=list(x=list(rot=0,cex=1),y=list(cex=1)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,1),box.ratio=1.5, col=c("grey","white"))
barchart(percentage~substrate|species, groups=context, data=substrate, layout=c(3,2),scale=list(x=list(rot=0,cex=1),y=list(cex=1)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,100),box.ratio=1.5, col=c("grey","white"))

pdf(file="Fig5Substrate_Percentage.pdf",width=14,height=8)
barchart(percentage~substrate|species, groups=context, data=substrate, layout=c(3,2),scale=list(x=list(rot=0,cex=0.75),y=list(cex=1)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,100),box.ratio=1.5, col=c("grey","white"))
dev.off()

pdf(file="Fig5Substrate_Proportion.pdf",width=14,height=8)
barchart(proportion~substrate|species, groups=context, data=substrate, layout=c(3,2),scale=list(x=list(rot=0,cex=0.75),y=list(cex=1)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,1),box.ratio=1.5, col=c("grey","white"))
dev.off()

#Note remember to order by aerial and non aerial for easier analyises

##############################################################################################
#Group size effect
########################################
#### Visualizing model fits
##############################################################################################
#Sethophaga fusca
#Basic glm plots
quartz(title="flock size vs. foraging rate")                                         # creates a quartz window with title
plot(jitter(capture,amount=0.9)~flocksizespecies,data=fusca,xlab="",ylab="", xlim=c(0,40))    # plot 
symbols(fusca2$foragingRate~fusca2$flocksizespecies, circles=fusca2$dayofseason, inches=0.15,pch=11, col= "black", bg="black",fg="black")
modelfusca<-glm(capture~flocksizespecies+I(flocksizespecies^2), data=fusca,family=poisson(link="log"), offset=log(minutes))
#model<-glm(capture~flocksizespecies+I(flocksizespecies^2)+offset(log(minutes)), data=fusca2,family=poisson(link="log")
curve(predict(modelfusca,data.frame(flocksizespecies=x,minutes=1),type="resp"), add=TRUE) # draws a curve based on prediction from regression model

### Using visreg

#You can use visreg to visualize the model fit on the transformed scale (the function uses predict(z) to generate the result). 
# The glm method fits a linear model on the transformed scale, and this is what you will visualize. The dots are not the transformed data, however. 
# They are "working values" obtained by transforming the residuals of the fitted model on the original scale. glm repeatedly recalculates the working values and 
#the fitted model as it converges on the maximum likelihood estimate. visreg shows you the results from the final iteration.

# If want to change the size of the points
#visreg(model,"flocksizespecies", type = "conditional", xlim=c(0,40), ylim=c(0,20), scale = "response", ylab="",  xlab=NA, rug=0) 
#par(new = TRUE)
#symbols(fusca$foragingRate~fusca$flocksizespecies, circles=fusca$dayofseason, inches=0.15,pch=11, col= "black", bg="black",fg="black",xlim=c(0,40), ylim=c(0,20))

visreg(modelfusca,xvar="flocksizespecies",trans=exp)
visreg(modelfusca,xvar="flocksizespecies")
visreg(modelfusca,xvar="flocksizespecies",scale="response")

pdf(file="Fig6afuscaflocksize.pdf",width=14,height=10,useDingbats=FALSE)
modelfusca<-glm(capture~flocksizespecies+I(flocksizespecies^2)+offset(log(minutes)), data=fusca,family=poisson(link="log"))
visreg(modelfusca,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,20), scale = "response", ylab=NA,  xlab=NA, rug=0, main="Setophaga fusca",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizespecies,1.5), data=fusca, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,40), ylim=c(0,20),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
dev.off()

#####################################################
#Cardellina canadensis
#####################################################
modelcanadensis<-glm(capture~flocksizespecies+I(flocksizespecies^2)+offset(log(minutes)), data=canadensis,family=poisson(link="log"))

pdf(file="Fig6acanadensisflocksize.pdf",width=14,height=10,useDingbats=FALSE)
visreg(modelcanadensis,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,20), ylab="", scale ="response", xlab=NA,ylab=NA, rug=0, main="Cardellina canadensis",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizespecies,1.5), data=canadensis, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,40), ylim=c(0,20),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
dev.off()

#####################################################
#Cerulea
#####################################################
modelcerulea<-glm(capture~flocksizespecies+offset(log(minutes)), data=cerulea,family=poisson(link="log"))

pdf(file="Fig6bceruleaflocksize.pdf",width=14,height=10,useDingbats=FALSE)
visreg(modelcerulea,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,20), scale = "response", ylab=NA,  xlab=NA, rug=0,main="Setophaga cerulea",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizespecies,1.5), data=cerulea, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,40), ylim=c(0,20),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
dev.off()
#####################################################
#peregrina
#####################################################

modelperegrina<-glm(capture~flocksizespecies+I(flocksizespecies^2)+offset(log(minutes)), data=peregrina,family=poisson(link="log"))
pdf(file="Fig6cperegrinaflocksize.pdf",width=14,height=10,useDingbats=FALSE)
visreg(modelperegrina,"flocksizespecies",cond=list(minutes=1), xlim=c(0,35), ylim=c(0,30), scale = "response", ylab=NA,  xlab=NA, rug=0, main="Oreothlypis peregrina",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizespecies,1.5), data=peregrina, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,35), ylim=c(0,30),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
dev.off()

#####################################################
#pitiayumi
#####################################################
###Note:Remove the influencial point to see it apply ylim=c(0,20)
modelpitiayumi<-glm(capture~flocksizespecies+I(flocksizespecies^2)+offset(log(minutes)), data=pitiayumi,family=poisson(link="log"))

pdf(file="Fig6apitiayumiflocksize.pdf",width=14,height=10,useDingbats=FALSE)
visreg(modelpitiayumi,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,15), scale = "response", ylab=NA,  xlab=NA, rug=0,main="Parula pitiayumi",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizespecies,1.5), data=pitiayumi, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,40), ylim=c(0,15),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
dev.off()

#####################################################
#guira
#####################################################

modelguira<-glm(capture~flocksizespecies+offset(log(minutes)), data=guira,family=poisson(link="log"))
pdf(file="Fig6eguiraflocksize.pdf",width=14,height=10,useDingbats=FALSE)
visreg(modelguira,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,15), scale = "response", rug=0,main="Hemithraupis guira",cex=1.5,cex.axis=1.5,ylab=NA,  xlab=NA) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizespecies,1.5), data=guira, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,40), ylim=c(0,15),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
dev.off()
#####################################################
#chrysops
#####################################################
modelchrysops<-glm(capture~flocksizespecies+I(flocksizespecies^2)+offset(log(minutes)), data=chrysops,family=poisson(link="log"))
pdf(file="Fig6fchrysopsflocksize.pdf",width=14,height=10,useDingbats=FALSE)
visreg(modelchrysops,"flocksizespecies",cond=list(minutes=1),xlim=c(0,40), ylim=c(0,15), scale = "response",rug=0,main="Zimmerius chrysops ",cex=1.5,cex.axis=1.5)
par(new = TRUE)
plot(foragingRate~jitter(flocksizespecies,1.5), data=chrysops, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,40), ylim=c(0,15),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
dev.off()



#USING GGPLOT
ggplot(fusca2, aes(x = flocksizespecies, y = foragingRate)) + #shoudl it be capture?
  geom_point(aes(size = dayofseason)) 
stat_smooth(method = "glm", family=quasipoisson, formula=y~x+I(x^2))  + 
  theme_bw() 

# other summary graphs no used in the manuscript ----------------------------------

# visualizing the data for each species
## boxplot with stripchart on the back for foraging data
#Migrants
#Setophaga fusca
boxplot(fusca$foragingrate~fusca$sociality, main='Setophaga fusca',ylab='Capture rate/min',col="white",ylim=c(0,20))
stripchart(fusca$foragingrate~fusca$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,20))
#Setophaga cerulea
boxplot(cerulea$foragingrate~cerulea$sociality, main='Setophaga cerulea',ylab='Capture rate/min',col="white",ylim=c(0,20), width=c(1.0, 1.0))
stripchart(cerulea$foragingrate~cerulea$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,20))
#Cardelina canadensis
boxplot(canadensis$foragingrate~canadensis$sociality, main='Cardelina canadensis',ylab='Capture rate/min',col="white",ylim=c(0,15))
stripchart(canadensis$foragingrate~canadensis$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,15))
#Oreothlypis peregrina
boxplot(peregrina$foragingrate~peregrina$sociality, main='Oreothlypis peregrina',ylab='Capture rate/min',col="white",ylim=c(0,25))
stripchart(peregrina$foragingrate~peregrina$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,25))
#Residents
#Hemithraupis guira
boxplot(guira$foragingrate~guira$sociality, main='Hemithraupis guira',ylab='Capture rate/min',col="white",ylim=c(0,15))
stripchart(guira$foragingrate~guira$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,15))
#Zimmerius chrysops
boxplot(chrysops$foragingrate~chrysops$sociality, main='Zimmerius chrysops',ylab='Capture rate/min',col="white",ylim=c(0,10))
stripchart(chrysops$foragingrate~chrysops$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,10))
#Parula pitiayumi
boxplot(pitiayumi$foragingrate~pitiayumi$sociality, main='Parula pitiayumi',ylab='Capture rate/min',col="white",ylim=c(0,20))
stripchart(pitiayumi$foragingrate~pitiayumi$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,20))

#Effect size
effect<-read.csv(file.choose("effect"), stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))  

str(effect)

plot(Benefit ~ Ocurrence, data = effect )
cor.test(effect$Benefit,effect$Ocurrence, method="kendall", alternative="greater")
model<-lm(Benefit~Ocurrence, data=effect)

plot(Benefit ~ Ocurrence, data = effect,xlab="Ocurrence in flocks", ylab="Increase in capture rate", xlim=c(0,10),ylim=c(0,15),pch=16, cex=1.5, cex.axis=1.5, cex.lab=1.5)
abline(model)
#Kendall's rank correlation tau

#data:  effect$Benefit and effect$Ocurrence
#T = 3, p-value = 0.8833
#alternative hypothesis: true tau is greater than 0
#sample estimates:
# tau 
#-0.4 



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
############Reviews from journal of Animal behavior###################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

foraging<-read.csv(file.choose("foraging_all_species_offset"), stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))  

str(foraging)
foraging$flocksizeind<- as.numeric(foraging$flocksizeind)

library(dplyr)

fusca<-filter(foraging, species=="Sethophaga fusca")
canadensis<-filter(foraging, species=="Cardelina canadensis")
cerulea<-filter(foraging,species=="Setophaga cerulea")
peregrina<-filter(foraging,species=="Oreothlypis peregrina")
guira<-filter(foraging,species=="Hemithraupis guira")
chrysops<-filter(foraging,species=="Zimmerius chrysops")
pitiayumi<-filter(foraging,species=="Parula pitiayumi")
canadensis<-filter(foraging, species=="Cardelina canadensis")

str(canadensis)

#Model selection for flock size
glm1fsi<-glm(capture~flocksizeind+I(flocksizeind^2), data=fusca,family=poisson(link="log"), offset=log(minutes))
glm1fsi<-glm(capture~flocksizeind+I(flocksizeind^2), data=canadensis,family=poisson(link="log"), offset=log(minutes))
glm1fsi<-glm(capture~flocksizeind+I(flocksizeind^2), data=cerulea,family=poisson(link="log"), offset=log(minutes)) #Overdisspersion to high model innacurat, small sample size
glm1fsi<-glm(capture~flocksizeind+I(flocksizeind^2), data=peregrina,family=poisson(link="log"), offset=log(minutes))
glm1fsi<-glm(capture~flocksizeind+I(flocksizeind^2), data=pitiayumi,family=poisson(link="log"), offset=log(minutes))
glm1fsi<-glm(capture~flocksizeind+I(flocksizeind^2), data=guira,family=poisson(link="log"), offset=log(minutes))
glm1fsi<-glm(capture~flocksizeind+I(flocksizeind^2), data=chrysops,family=poisson(link="log"), offset=log(minutes))

glm2fsi<-glm(capture~flocksizeind, data=canadensis,family=poisson(link="log"), offset=log(minutes))
glm3fsi<-update(glm2fs, . ~ . - flocksizeind)
glm4fsi<-glm(capture~I(flocksizeind^2), data=canadensis,family=poisson(link="log"), offset=log(minutes))
glmQ1fsi<-update(glm1fs, family=quasipoisson(link="log"))
glmQ2fsi<-update(glm2fs, family=quasipoisson(link="log"))
glmQ3fsi<-update(glm3fs, family=quasipoisson(link="log"))
glmQ4fsi<-update(glm4fs, family=quasipoisson(link="log"))

# Summary of the model allow me to interpret the estimates of the parameters (e.g effect sizes) in the model and the difference from cero and between them.
# The summary give me the effect size. 
summary(glm1fsi)
summary(glmQ1fsi)


#####overdisperssion parameters 
dfun<-function(glmQ1fsi){with(glmQ1fsi,sum((weights * residuals^2)[weights > 0])/df.residual)}
dfun(glmQ1fsi)


#fusca=2.50
#cerulea=4
#peregrina=4
#canadensis=1.76
#pitiayumi=2.75
#guira=3.17
#chrysops=

#### Model selection using  AICcmodavg
library(AICcmodavg)

###Model used  (glm1,glm2,glm3,glm4), for all the species which have overdispersed datta, remember to use the glm models instead of glmQ
aictab(list(glm1fs,glm2fs,glm3fs,glm4fs),
       modnames=c("FlockSizeSpecies+POLY",
                  "FlockSizeSpecies",
                  "Intercept",
                  "poly"),
       c.hat=2.50)

aictab(list(glm1fs,glm2fs,glm3fs),
       modnames=c("FlockSizeSpecies+POLY",
                  "FlockSizeSpecies",
                  "Intercept"),
       c.hat=1.05)


#### To make the table for AIC cfor species that are not overdisperssed

aictab (list(glm1,glm2fs,glm3fs,glm4fs,glm5fs),
        modnames=c("FlockSizeSpecies+Dayofseason",
                   "FlockSizeSpecies",
                   "Dayofseason",
                   "Intercept",
                   "poly"))

######SUMMARY TABLE
#A beautiful table in html for the parameters of glm models
#For the parameter stimates we need to consider the models with the Quasipoisson distribution, the points estimates are identical
# to the model with poisson distribution but the standard error and confidence intervals are wider!
sjt.glm(glmQ1fs,glmQ2fs,glmQ4fs,
        depvar.labels = c("Model1: Socialcontext+Dayofseason","Model2:Socialcontext","Model3:Day of season"),
        pred.labels = NULL,
        show.aicc=TRUE, 
        show.family=TRUE, 
        group.pred = FALSE,
        exp.coef = FALSE,  # if true Regression and cof.intervals are exponentiaded st.error are not in the unstransformed scale
        p.numeric=TRUE,
        robust = TRUE,
        show.se = TRUE,
        show.r2=TRUE,
        show.dev=TRUE,
        show.chi2=TRUE,
        cell.spacing = 0.001,
        sep.column = FALSE)






#Group size effect using 
########################################
#### Visualizing model fits
##############################################################################################
#Sethophaga fusca
#Basic glm plots
quartz(title="flock size vs. foraging rate")                                         # creates a quartz window with title
plot(jitter(capture,amount=0.9)~flocksizeind,data=fusca,xlab="",ylab="", xlim=c(0,40))    # plot 
symbols(fusca$foragingRate~fusca$flocksizeind, circles=fusca$dayofseason, inches=0.15,pch=11, col= "black", bg="black",fg="black")
modelfusca<-glm(capture~flocksizeind+I(flocksizeind^2), data=fusca,family=poisson(link="log"), offset=log(minutes))
model<-glm(capture~flocksizeind+I(flocksizeind^2)+offset(log(minutes)), data=fusca2,family=poisson(link="log"),curve(predict(modelfusca,data.frame(flocksizeind=x,minutes=1),type="resp"), add=TRUE)) # draws a curve based on prediction from regression model

### Using visreg

#You can use visreg to visualize the model fit on the transformed scale (the function uses predict(z) to generate the result). 
# The glm method fits a linear model on the transformed scale, and this is what you will visualize. The dots are not the transformed data, however. 
# They are "working values" obtained by transforming the residuals of the fitted model on the original scale. glm repeatedly recalculates the working values and 
#the fitted model as it converges on the maximum likelihood estimate. visreg shows you the results from the final iteration.

# If want to change the size of the points
#visreg(model,"flocksizeind", type = "conditional", xlim=c(0,40), ylim=c(0,20), scale = "response", ylab="",  xlab=NA, rug=0) 
#par(new = TRUE)
#symbols(fusca$foragingRate~fusca$flocksizeind, circles=fusca$dayofseason, inches=0.15,pch=11, col= "black", bg="black",fg="black",xlim=c(0,40), ylim=c(0,20))

visreg(modelfusca,xvar="flocksizeind",trans=exp)
visreg(modelfusca,xvar="flocksizeind")
visreg(modelfusca,xvar="flocksizeind",scale="response")

pdf(file="Fig11.fuscaflocksizeind.pdf",width=14,height=10,useDingbats=FALSE)
modelfusca<-glm(capture~flocksizeind+I(flocksizeind^2)+offset(log(minutes)), data=fusca,family=poisson(link="log"))
visreg(modelfusca,"flocksizeind",cond=list(minutes=1), xlim=c(0,60), ylim=c(0,20), scale = "response", ylab=NA,  xlab=NA, rug=0, main="Setophaga fusca",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizeind,1.5), data=fusca, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,60), ylim=c(0,20),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
dev.off()

#####################################################
#Cardellina canadensis
#####################################################
modelcanadensis<-glm(capture~flocksizeind+I(flocksizeind^2)+offset(log(minutes)), data=canadensis,family=poisson(link="log"))

pdf(file="Fig11.canadensisflocksize.pdf",width=14,height=10,useDingbats=FALSE)
visreg(modelcanadensis,"flocksizeind",cond=list(minutes=1), xlim=c(0,60), ylim=c(0,20), ylab="", scale ="response", xlab=NA,ylab=NA, rug=0, main="Cardellina canadensis",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizeind,1.5), data=canadensis, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,60), ylim=c(0,20),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
dev.off()

#####################################################
#Cerulea
#####################################################
modelcerulea<-glm(capture~flocksizeind+offset(log(minutes)), data=cerulea,family=poisson(link="log"))

pdf(file="Fig11.ceruleaflocksize.pdf",width=14,height=10,useDingbats=FALSE)
visreg(modelcerulea,"flocksizeind",cond=list(minutes=1), xlim=c(0,70), ylim=c(0,20), scale = "response", ylab=NA,  xlab=NA, rug=0,main="Setophaga cerulea",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizeind,1.5), data=cerulea, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,70), ylim=c(0,20),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
dev.off()
#####################################################
#peregrina
#####################################################

modelperegrina<-glm(capture~flocksizeind+I(flocksizeind^2)+offset(log(minutes)), data=peregrina,family=poisson(link="log"))
pdf(file="Fig11peregrinaflocksize.pdf",width=14,height=10,useDingbats=FALSE)
visreg(modelperegrina,"flocksizeind",cond=list(minutes=1), xlim=c(0,60), ylim=c(0,30), scale = "response", ylab=NA,  xlab=NA, rug=0, main="Oreothlypis peregrina",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizeind,1.5), data=peregrina, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,60), ylim=c(0,30),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
dev.off()

#####################################################
#pitiayumi
#####################################################
###Note:Remove the influencial point to see it apply ylim=c(0,20)
modelpitiayumi<-glm(capture~flocksizeind+I(flocksizeind^2)+offset(log(minutes)), data=pitiayumi,family=poisson(link="log"))

pdf(file="Fig11pitiayumiflocksize.pdf",width=14,height=10,useDingbats=FALSE)
visreg(modelpitiayumi,"flocksizeind",cond=list(minutes=1), xlim=c(0,60), ylim=c(0,15), scale = "response", ylab=NA,  xlab=NA, rug=0,main="Parula pitiayumi",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizeind,1.5), data=pitiayumi, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,60), ylim=c(0,15),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
dev.off()

#####################################################
#guira
#####################################################

modelguira<-glm(capture~flocksizeind+offset(log(minutes)), data=guira,family=poisson(link="log"))
pdf(file="Fig11guiraflocksize.pdf",width=14,height=10,useDingbats=FALSE)
visreg(modelguira,"flocksizeind",cond=list(minutes=1), xlim=c(0,60), ylim=c(0,15), scale = "response", rug=0,main="Hemithraupis guira",cex=1.5,cex.axis=1.5,ylab=NA,  xlab=NA) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizeind,1.5), data=guira, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,60), ylim=c(0,15),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
dev.off()
#####################################################
#chrysops
#####################################################
modelchrysops<-glm(capture~flocksizeind+I(flocksizeind^2)+offset(log(minutes)), data=chrysops,family=poisson(link="log"))
pdf(file="Fig11.chrysopsflocksize.pdf",width=14,height=10,useDingbats=FALSE)
visreg(modelchrysops,"flocksizeind",cond=list(minutes=1),xlim=c(0,50), ylim=c(0,15), scale = "response",rug=0,main="Zimmerius chrysops ",cex=1.5,cex.axis=1.5)
par(new = TRUE)
plot(foragingRate~jitter(flocksizeind,1.5), data=chrysops, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,50), ylim=c(0,15),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
dev.off()



#USING GGPLOT
ggplot(fusca2, aes(x = flocksizeind, y = foragingRate)) + #shoudl it be capture?
  geom_point(aes(size = dayofseason)) 
stat_smooth(method = "glm", family=quasipoisson, formula=y~x+I(x^2))  + 
  theme_bw() 

# other summary graphs no used in the manuscript ----------------------------------

# visualizing the data for each species
## boxplot with stripchart on the back for foraging data
#Migrants
#Setophaga fusca
boxplot(fusca$foragingrate~fusca$sociality, main='Setophaga fusca',ylab='Capture rate/min',col="white",ylim=c(0,20))
stripchart(fusca$foragingrate~fusca$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,20))
#Setophaga cerulea
boxplot(cerulea$foragingrate~cerulea$sociality, main='Setophaga cerulea',ylab='Capture rate/min',col="white",ylim=c(0,20), width=c(1.0, 1.0))
stripchart(cerulea$foragingrate~cerulea$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,20))
#Cardelina canadensis
boxplot(canadensis$foragingrate~canadensis$sociality, main='Cardelina canadensis',ylab='Capture rate/min',col="white",ylim=c(0,15))
stripchart(canadensis$foragingrate~canadensis$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,15))
#Oreothlypis peregrina
boxplot(peregrina$foragingrate~peregrina$sociality, main='Oreothlypis peregrina',ylab='Capture rate/min',col="white",ylim=c(0,25))
stripchart(peregrina$foragingrate~peregrina$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,25))
#Residents
#Hemithraupis guira
boxplot(guira$foragingrate~guira$sociality, main='Hemithraupis guira',ylab='Capture rate/min',col="white",ylim=c(0,15))
stripchart(guira$foragingrate~guira$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,15))
#Zimmerius chrysops
boxplot(chrysops$foragingrate~chrysops$sociality, main='Zimmerius chrysops',ylab='Capture rate/min',col="white",ylim=c(0,10))
stripchart(chrysops$foragingrate~chrysops$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,10))
#Parula pitiayumi
boxplot(pitiayumi$foragingrate~pitiayumi$sociality, main='Parula pitiayumi',ylab='Capture rate/min',col="white",ylim=c(0,20))
stripchart(pitiayumi$foragingrate~pitiayumi$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,20))


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
############Reviews for JFO###################
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


foraging<-read.csv(file.choose("foraging_all_species_offset"), stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))  

str(foraging)
foraging$flocksizeind<- as.numeric(foraging$flocksizeind)

library(dplyr)

fusca<-filter(foraging, species=="Sethophaga fusca")
canadensis<-filter(foraging, species=="Cardelina canadensis")
cerulea<-filter(foraging,species=="Setophaga cerulea")
peregrina<-filter(foraging,species=="Oreothlypis peregrina")
guira<-filter(foraging,species=="Hemithraupis guira")
chrysops<-filter(foraging,species=="Zimmerius chrysops")
pitiayumi<-filter(foraging,species=="Parula pitiayumi")
canadensis<-filter(foraging, species=="Cardelina canadensis")


#Calculate the average foraging rate 
flocks<-filter(foraging, sociality=="Flock")
mean(flocks$foragingRate)
#4.878172
mean(flocks$movementRate)
#21.83333908

solitary<-filter(foraging, sociality=="Solitary")
mean(solitary$foragingRate)
#1.286309
mean(solitary$movementRate)
#

#a)#####Relationship of foraging and seasonality

canadensisplot<-ggplot(data=canadensis, aes(x=dayofseason, y=foragingRate))+
   geom_point(data=canadensis, aes(x=dayofseason, y=foragingRate),size=4)+
    geom_smooth(method = "glm")+
    ggtitle("Canada Warbler")+
    labs(y = " Attack rate (number of attacks/min", x = "Day of the season", main="canadensis")+
  xlim(0, 80)+
  ylim(0,20)+
   theme_classic()+
  My_theme

fuscaplot<-ggplot(data=fusca, aes(x=dayofseason, y=foragingRate))+
     geom_point(data=fusca, aes(x=dayofseason, y=foragingRate),size=4)+
     geom_smooth(method = "glm")+
     ggtitle("Bluckburnian Warbler")+
    labs(y = " Attack rate (number of attacks/min", x = "Day of the season")+
  xlim(0, 80)+
  ylim(0,20)+
     theme_classic()+ My_theme
  
ceruleaplot<-ggplot(data=cerulea, aes(x=dayofseason, y=foragingRate))+
    geom_point(data=cerulea, aes(x=dayofseason, y=foragingRate),size=4)+
    geom_smooth(method = "glm")+
    ggtitle("Cerulean Warbler")+
    labs(y = " Attack rate (number of attacks/min)", x = "Day of the season")+
  xlim(0, 80)+
  ylim(0,20)+
   theme_classic()+ My_theme
peregrinaplot<-ggplot(data=peregrina, aes(x=dayofseason, y=foragingRate))+
   geom_point(data=peregrina, aes(x=dayofseason, y=foragingRate),size=4)+
   geom_smooth(method = "glm")+
   ggtitle("Tennessee Warbler")+
   labs(y = " Attack rate (number of attacks/min)", x = "Day of the season")+
  xlim(0, 80)+
  ylim(0,30)+
   theme_classic()+ My_theme
guiraplot<-ggplot(data=guira, aes(x=dayofseason, y=foragingRate))+
    geom_point(data=guira, aes(x=dayofseason, y=foragingRate),size=4)+
   geom_smooth(method = "glm")+
   ggtitle("Guira Tanager")+
   labs(y = " Attack rate (number of attacks/min)", x = "Day of the season")+
  xlim(0, 80)+
  ylim(0,15)+
   theme_classic()+ My_theme
chrysopsplot<-ggplot(data=chrysops, aes(x=dayofseason, y=foragingRate))+
   geom_point(data=chrysops, aes(x=dayofseason, y=foragingRate),size=4)+
   geom_smooth(method = "glm")+
   ggtitle("Golden-faced tyrannulet")+
   labs(y = " Attack rate (number of attacks/min)", x = "Day of the season")+
  xlim(0, 80)+
  ylim(0,15)+
   theme_classic()+ My_theme
pitiayumiplot<-ggplot(data=pitiayumi, aes(x=dayofseason, y=foragingRate))+
  geom_point(data=pitiayumi, aes(x=dayofseason, y=foragingRate),size=4)+
    geom_smooth(method = "glm")+
   ggtitle("Tropical Parula")+
   labs(y = " Attack rate (number of attacks/min)", x = "Day of the season")+
  xlim(0, 80)+
  ylim(0,15)+
   theme_classic()+ My_theme
grid.arrange(fuscaplot, canadensisplot,peregrinaplot,ceruleaplot,guiraplot,chrysopsplot,pitiayumiplot,  nrow = 2, top="Seasonality and attack rate")

grid.arrange

#b)####Relationship of movement and foraging
My_theme = theme( 
  legend.position="none",
  axis.title.x = element_text(size = 30),
  axis.text.y = element_text(size = 30, color="black",angle=0),
  axis.text.x =  element_text(size = 30, color="black",angle=0),
  axis.title.y = element_text(size = 30))

canadensisplot_mvt_foraging<-ggplot(data=canadensis, aes(x=movementRate , y=foragingRate, color=sociality))+
  scale_color_manual(values=c("black","#999999"))+
  geom_point(data=canadensis, aes(x=movementRate , y=foragingRate),size=4)+
  geom_smooth(method = "glm")+
  ggtitle("Canada Warbler")+
  labs(y = " Attack rate (number of attacks/min)", x = "Movement rate")+
  xlim(0, 40)+
  ylim(0,20)+
  theme_classic()+
  My_theme

fuscaplot_mvt_foraging<-ggplot(data=fusca, aes(x=movementRate , y=foragingRate, color=sociality))+
  scale_color_manual(values=c("black","#999999"))+
  geom_point(data=fusca, aes(x=movementRate , y=foragingRate),size=4)+
  geom_smooth(method = "glm")+
  ggtitle("Blackburnian Warbler")+
  labs(y = " Attack rate (number of attacks/min)", x = "Movement rate")+
  xlim(0, 40)+
  ylim(0,20)+
  theme_classic()+
  My_theme

ceruleaplot_mvt_foraging<-ggplot(data=cerulea, aes(x=movementRate , y=foragingRate, color=sociality))+
  scale_color_manual(values=c("black","#999999"))+
  geom_point(data=cerulea, aes(x=movementRate , y=foragingRate),size=4)+
  geom_smooth(method = "glm")+
  ggtitle("Cerulean Warbler")+
  labs(y = " Attack rate (number of attacks/min)", x = "Movement rate")+
  xlim(0, 40)+
  ylim(0,20)+
  theme_classic()+
  My_theme

 peregrinaplot_mvt_foraging<-ggplot(data=peregrina, aes(x=movementRate , y=foragingRate, color=sociality))+
  scale_color_manual(values=c("black","#999999"))+
  geom_point(data=peregrina, aes(x=movementRate , y=foragingRate),size=4)+
  geom_smooth(method = "glm")+
  ggtitle("Tennessee Warbler")+
  labs(y = " Attack rate (number of attacks/min)", x = "Movement rate")+
   xlim(0, 40)+
   ylim(0,30)+
  theme_classic()+
  My_theme
 
guiraplot_mvt_foraging<-ggplot(data=guira, aes(x=movementRate , y=foragingRate, color=sociality))+
   scale_color_manual(values=c("black","#999999"))+
   geom_point(data=guira, aes(x=movementRate , y=foragingRate),size=4)+
   geom_smooth(method = "glm")+
   ggtitle("Guira Tanager")+
   labs(y = " Attack rate (number of attacks/min)", x = "Movement rate")+
   xlim(0, 40)+
  ylim(0,15)+
   theme_classic()+
   My_theme

chrysopsplot_mvt_foraging<-ggplot(data=chrysops, aes(x=movementRate , y=foragingRate, color=sociality))+
  scale_color_manual(values=c("black","#999999"))+
  geom_point(data=chrysops, aes(x=movementRate , y=foragingRate),size=4)+
  geom_smooth(method = "glm")+
  ggtitle("Golden-faced tyrannuletr")+
  labs(y = " Attack rate (number of attacks/min)", x = "Movement rate")+
  xlim(0, 25)+
  ylim(0,15)+
  theme_classic()+
  My_theme

pitiayumiplot_mvt_foraging<-ggplot(data= pitiayumi, aes(x=movementRate , y=foragingRate, color=sociality))+
  scale_color_manual(values=c("black","#999999"))+
  geom_point(data= pitiayumi, aes(x=movementRate , y=foragingRate),size=4)+
  geom_smooth(method = "glm")+
  ggtitle("Tropical Parula")+
  labs(y = " Attack rate (number of attacks/min)", x = "Movement rate")+
  xlim(0, 40)+
  ylim(0,15)+
  theme_classic()+
  My_theme

grid.arrange(fuscaplot_mvt_foraging, canadensisplot_mvt_foraging,peregrinaplot_mvt_foraging,ceruleaplot_mvt_foraging,guiraplot_mvt_foraging,chrysopsplot_mvt_foraging,pitiayumiplot_mvt_foraging,  nrow = 2, top="Seasonality and attack rate")


#c)#####Combining the two movement and foraging rate models for each species

#####################################################
#Cardellina canadensis
#####################################################

str(canadensis)
modelcanadensis<-glm(capture~flocksizespecies+I(flocksizespecies^2)+offset(log(minutes)), data=canadensis,family=poisson(link="log"))

pdf(file="Fig6acanadensisflocksize.pdf",width=14,height=10,useDingbats=FALSE)
visreg(modelcanadensis,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,20), ylab="", scale ="response", xlab=NA,ylab=NA, rug=0, main="Cardellina canadensis",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizespecies,1.5), data=canadensis, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,40), ylim=c(0,40),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
dev.off()

modelcanadensis1<-glm(movement~flocksizespecies+I(flocksizespecies^2)+offset(log(minutes)), data=canadensis,family=poisson(link="log"))
visreg(modelcanadensis1,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,40), ylab="", scale ="response", xlab=NA,ylab=NA, rug=0, main="Cardellina canadensis",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(movementRate~jitter(flocksizespecies,1.5), data=canadensis, ylab="Movement rate (number of moves / min)", xlab="Flock size (Number of species)",xlim=c(0,40), ylim=c(0,40),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
par(new = TRUE)


#Combined plots

visreg(modelcanadensis,"flocksizespecies",cond=list(minutes=1), xlim=c(0,35), ylim=c(0,40), ylab="", scale ="response", xlab=NA,ylab=NA, rug=0, main="",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
visreg(modelcanadensis1,"flocksizespecies",cond=list(minutes=1), xlim=c(0,35), ylim=c(0,40), ylab="", scale ="response", xlab=NA,ylab=NA, rug=0, main="Canada Warbler",cex=1.5,cex.axis=1.5, xlab="Flock diversity (number of species)",ylab="Attack rate (number of attacks / min)") 

#####################################################
#Cerulea
#####################################################
modelcerulea<-glm(capture~flocksizespecies+offset(log(minutes)), data=cerulea,family=poisson(link="log"))
visreg(modelcerulea,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,40), scale = "response", ylab=NA,  xlab=NA, rug=0,main="Setophaga cerulea",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizespecies,1.5), data=cerulea, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,40), ylim=c(0,40),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)

modelcerulea1<-glm(movement~flocksizespecies+offset(log(minutes)), data=cerulea,family=poisson(link="log"))
visreg(modelcerulea1,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,40), scale = "response", ylab=NA,  xlab=NA, rug=0,main="Setophaga cerulea",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(movementRate~jitter(flocksizespecies,1.5), data=cerulea, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,40), ylim=c(0,40),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)

#Combined plots
visreg(modelcerulea,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,35), ylab="", scale ="response", xlab=NA,ylab=NA, rug=0, main="Cerulean Warbler",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
visreg(modelcerulea1,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,35), scale ="response", rug=0, main="Cerulean Warbler",cex=1.5,cex.axis=1.5, xlab="Flock diversity (number of species)", ylab="Attack rate (number of attacks / min)")


#####################################################
#peregrina
#####################################################

modelperegrina<-glm(capture~flocksizespecies+I(flocksizespecies^2)+offset(log(minutes)), data=peregrina,family=poisson(link="log"))
visreg(modelperegrina,"flocksizespecies",cond=list(minutes=1), xlim=c(0,35), ylim=c(0,30), scale = "response", ylab=NA,  xlab=NA, rug=0, main="Oreothlypis peregrina",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizespecies,1.5), data=peregrina, ylab="Foraging rate (#captures / min)", xlab="Flock diversity (Number of species)",xlim=c(0,35), ylim=c(0,30),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)

modelperegrina1<-glm(movement~flocksizespecies+offset(log(minutes)), data=peregrina,family=poisson(link="log"))
visreg(modelperegrina1,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,40), scale = "response", ylab=NA,  xlab=NA, rug=0,main="Tennessee Warbler",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(movementRate~jitter(flocksizespecies,1.5), data=peregrina, ylab="Foraging rate (number of attack/min)", xlab="Flock diversity (umber of species)",xlim=c(0,40), ylim=c(0,40),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)

#Combined plots
visreg(modelperegrina,"flocksizespecies",cond=list(minutes=1), xlim=c(0,35), ylim=c(0,30), scale = "response", ylab=NA,  xlab=NA, rug=0,cex=1.5,cex.axis=1.5)
par(new = TRUE)
visreg(modelperegrina1,"flocksizespecies",cond=list(minutes=1), xlim=c(0,35), ylim=c(0,30), scale = "response", ylab=NA,  xlab=NA, rug=0,main="Tennessee Warbler",cex=1.5,cex.axis=1.5,xlab="Flock diversity (number of species)",ylab="Attack rate (number of attacks / min)") 


#####################################################
#fusca
#####################################################

modelfusca<-glm(capture~flocksizespecies+I(flocksizespecies^2)+offset(log(minutes)), data=fusca,family=poisson(link="log"))
visreg(modelfusca,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,40), scale = "response", ylab=NA,  xlab=NA, rug=0, main="Blackburnian Warbler",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizespecies,1.5), data=fusca, ylab="Foraging rate (#captures / min)", xlab="Flock diversity (Number of species)",xlim=c(0,40), ylim=c(0,40),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)

modelfusca1<-glm(movement~flocksizespecies+offset(log(minutes)), data=fusca,family=poisson(link="log"))
visreg(modelfusca1,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,40), scale = "response", ylab=NA,  xlab=NA, rug=0,main="Blackburnian Warbler",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(movementRate~jitter(flocksizespecies,1.5), data=fusca, ylab="Attack rate (number of attacks/min)", xlab="Flock diversity (number of species)",xlim=c(0,40), ylim=c(0,40),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)

#Combined plots
visreg(modelfusca,"flocksizespecies",cond=list(minutes=1), xlim=c(0,35), ylim=c(0,30), scale = "response", ylab=NA,  xlab=NA, rug=0,cex=1.5,cex.axis=1.5)
par(new = TRUE)
visreg(modelfusca1,"flocksizespecies",cond=list(minutes=1), xlim=c(0,35), ylim=c(0,30), scale = "response", rug=0,main="Blackburnian Warbler",cex=1.5,cex.axis=1.5, xlab="Flock diversity (number of species)",ylab="Attack rate (number of attacks / min)") 


#####################################################
#pitiayumi
#####################################################
###Note:Remove the influencial point to see it apply ylim=c(0,20)
modelpitiayumi<-glm(capture~flocksizespecies+I(flocksizespecies^2)+offset(log(minutes)), data=pitiayumi,family=poisson(link="log"))
visreg(modelpitiayumi,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,15), scale = "response", ylab=NA,  xlab=NA, rug=0,main="Parula pitiayumi",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizespecies,1.5), data=pitiayumi, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,40), ylim=c(0,15),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)

modelpitiayumi1<-glm(movement~flocksizespecies+offset(log(minutes)), data=pitiayumi,family=poisson(link="log"))
visreg(modelpitiayumi1,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,40), scale = "response", ylab=NA,  xlab=NA, rug=0,main="Parula pitiayumi",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(movementRate~jitter(flocksizespecies,1.5), data=pitiayumi, ylab="", xlab="Flock diversity (number of species)",xlim=c(0,40), ylim=c(0,40),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)

visreg(modelpitiayumi,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,40), scale = "response", ylab=NA,  xlab=NA, rug=0,main="Parula pitiayumi",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
visreg(modelpitiayumi1,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,40), scale = "response", ylab=NA,  xlab="Flock diversity (number of species)", ylab="Attack rate (number of attacks / min)",rug=0,main="Parula pitiayumi",cex=1.5,cex.axis=1.5) 



#####################################################
#guira
#####################################################

modelguira<-glm(capture~flocksizespecies+offset(log(minutes)), data=guira,family=poisson(link="log"))
visreg(modelguira,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,15), scale = "response", rug=0,main="Hemithraupis guira",cex=1.5,cex.axis=1.5,ylab=NA,  xlab=NA) 
par(new = TRUE)
plot(foragingRate~jitter(flocksizespecies,1.5), data=guira, ylab="Attack rate (number of attacks/min)", xlab="Flock size (Number of species)",xlim=c(0,40), ylim=c(0,15),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)

modelguira1<-glm(movement~flocksizespecies+offset(log(minutes)), data=guira,family=poisson(link="log"))
visreg(modelguira1,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,40), scale = "response", ylab=NA,  xlab=NA, rug=0,main="Hemithraupis guira",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(movementRate~jitter(flocksizespecies,1.5), data=guira, ylab="Attack rate (number of attacks/min)", xlab="Flock diversity (number of species)",xlim=c(0,40), ylim=c(0,40),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)

visreg(modelguira1,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,40), scale = "response", ylab=NA,  xlab=NA, rug=0,main="",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
visreg(modelguira,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,40), scale = "response", rug=0,main="Guira Tanager",cex=1.5,cex.axis=1.5,ylab=NA,xlab="Flock diversity (number of species)", ylab="Attack rate (number of attacks / min)") 

Fig10.Guira_combined

#####################################################
#chrysops
#####################################################
modelchrysops<-glm(capture~flocksizespecies+I(flocksizespecies^2)+offset(log(minutes)), data=chrysops,family=poisson(link="log"))
visreg(modelchrysops,"flocksizespecies",cond=list(minutes=1),xlim=c(0,40), ylim=c(0,15), scale = "response",rug=0,main="Zimmerius chrysops ",cex=1.5,cex.axis=1.5)
par(new = TRUE)
plot(foragingRate~jitter(flocksizespecies,1.5), data=chrysops, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,40), ylim=c(0,15),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)

modelchrysops1<-glm(movement~flocksizespecies+offset(log(minutes)), data=chrysops,family=poisson(link="log"))
visreg(modelchrysops1,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,40), scale = "response", ylab=NA,  xlab=NA, rug=0,main="Zimmerius chrysops",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(movementRate~jitter(flocksizespecies,1.5), data=chrysops, ylab="Attack rate (number of attacks/min)", xlab="Flock diversity (number of species)",xlim=c(0,40), ylim=c(0,40),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)

visreg(modelchrysops,"flocksizespecies",cond=list(minutes=1),xlim=c(0,40), ylim=c(0,15), scale = "response",rug=0,main="",cex=1.5,cex.axis=1.5,xlab="Flock diversity (number of species)", ylab="Attack rate (number of attacks / min)")
par(new = TRUE)
visreg(modelchrysops1,"flocksizespecies",cond=list(minutes=1), xlim=c(0,40), ylim=c(0,15), scale = "response", ylab=NA,  xlab=NA, rug=0,main="Golden-faced tyrannulet",cex=1.5,cex.axis=1.5) 






flocks_per<-filter(flocks, sociality=="")

str (peregrina)

peregrina_flocks<-filter(peregrina,sociality=="Flocks")


modelperegrina<-glm(movement~flocksizespecies+I(flocksizespecies^2)+offset(log(minutes)), data=peregrina_flocks,family=poisson(link="log"))
pdf(file="Fig1peregrinamvt.pdf",width=14,height=10,useDingbats=FALSE)
visreg(modelperegrina,"flocksizespecies",cond=list(minutes=1), xlim=c(0,60), ylim=c(0,50), scale = "response", ylab=NA,  xlab=NA, rug=0, main="Oreothlypis peregrina",cex=1.5,cex.axis=1.5) 
par(new = TRUE)
plot(movementRate ~jitter(flocksizespecies,1.5), data=peregrina_flocks, ylab="Foraging rate (#captures / min)", xlab="Flock size (Number of species)",xlim=c(0,60), ylim=c(0,50),pch=16,cex=1.5, cex.axis=1.5,cex.lab=1.5)
dev.off()



flocksizespecies
#############################END OF SCRIPT ##########################################
########################### END OF SCRIPt FOR MANUSCRIPT############################


## Using glm flocksize #species including flock an d solitary individuals includeit in the paper
ggplot(fusca2, aes(x = flocksizespecies, y = foragingRate)) + #shoudl it be capture?
geom_point(aes(size = dayofseason)) 
stat_smooth(method = "glm", family=quasipoisson, formula=y~x+I(x^2))  + 
theme_bw() 
ggplot(fusca2, aes(x = flocksizespecies, y = capture)) + #shoul it be capture?
geom_point(aes(size = dayofseason)) + 
geom_smooth(method = "glm", formula=y~x+I(x^2))  + 
theme_bw()

# individuals
ggplot(foraging, aes(x = flocksizeind, y = foragingRate), color="grey") +
geom_point(aes(size = dayofseason)) + 
geom_smooth(method = "glm") + 
facet_wrap(~species, scales = "free") + 
theme_bw()

ggplot(canadensis2, aes(x = flocksizespecies, y = foragingRate)) + #shoudl it be capture?
geom_point(aes(size = dayofseason)) + 
stat_smooth(method = "glm", formula=y~x+I(x^2))+ 
facet_wrap(~species, scales = "free")+
theme_bw()

#including only individuals in flocks
ggplot(flocks, aes(x = flocksizespecies, y = capture)) +
geom_point(aes(size = dayofseason)) + 
geom_smooth(method = "glm", formula=y~x+I(x^2))+ 
facet_wrap(~species, scales = "free")+ 
theme_bw()
### Using say of the saseon in the x axis
ggplot(foraging, aes(x = dayofseason, y = foragingRate)) +
geom_point(aes(size = flocksizeind)) + 
geom_smooth(method = "glm")+ 
facet_wrap(~species, scales = "free") + 
theme_bw()
# free smooth
ggplot(foraging, aes(x =flocksizespecies , y = foragingRate)) +
geom_point(aes(size = dayofseason), alpha = 0.6) + 
# geom_smooth(method = "lm") + 
geom_smooth() + 
facet_wrap(~species, scales = "free") + 
theme_bw()

##############################################################################################


#############################Previous analyses##########################################


##Prepring the da Variables as factor
### convert variables 
foraging$foragingrate<- as.numeric(foraging$foragingrate)
##foraging$foragingrate<- as.integer(foraging$foragingrate) # because I will need data to be discrete to use a poisson distribution for count data 
foraging$movementrate<- as.numeric(foraging$movementrate)
foraging$sociality<- as.factor(foraging$sociality)
foraging$flocksizespecies<- as.numeric(foraging$flocksizespecies)
foraging$flocksizeind<- as.numeric(foraging$flocksizeind)

#Filter the data using dplyr for each species
foraging %>%
  filter(species == "Setophaga fusca")

fusca<-filter(foraging,species=="Setophaga fusca")
canadensis<-filter(foraging, species=="Cardelina canadensis")
cerulea<-filter(foraging,species=="Setophaga cerulea")
peregrina<-filter(foraging,species=="Oreothlypis peregrina")
guira<-filter(foraging,species=="Hemithraupis guira")
chrysops<-filter(foraging,species=="Zimmerius chrysops")
pitiayumi<-filter(foraging,species=="Parula pitiayumi")


# Visualizing the data ----------------------------------------------------

# Combining the plots in a multiplot graph?

ggplot(foraging, aes(flocksizespecies, foragingrate)) + xlab("X variable name") + ylab("Y variable name") +
  geom_point(col = "black", size = I(2)) +
  geom_smooth(method = lm, size = I(1), se = FALSE, col = "black") +
  facet_wrap(~species, ncol = 0)



# Variables as factor
### convert variables 
foraging$foragingrate<- as.numeric(foraging$foragingrate)
foraging$sociality<- as.factor(foraging$sociality)
foraging$flocksizespecies<- as.numeric(foraging$flocksizespecies)


##### Evaluate the fit of Models

help(family)
stripchart (foragingrate~sociality, vertical=TRUE, data=foraging,pch=16,col=c("seagreen4","seagreen2","turquoise","skyblue1","royalblue","darkblue"),method="jitter",jitter=0.1,ylab="",xlab="")
points(c(1,2,3,4,5,6),mean, pch="-",col="red",cex=5)
stripchart(fitted(glm1)~foraging$foraging, vertical=TRUE,add=TRUE, pch="-",cex=2,method="jitter",col="black",na.strings=c("NA",""))
stripchart(fitted(glm2)~foraging$sociality, vertical=TRUE,add=TRUE, pch="-",cex=2, method="jitter",col="tan1")
stripchart(fitted(glm3)~foraging$sociality, vertical=TRUE,add=TRUE, pch="-",cex=2,method="jitter",col="gray")
legend("bottomright", legend=c("additive model","dominance model","genotype model","observed mean"),bty="n",lwd=2,cex=0.8, col=c("black","tan1","gray","red"), lty=c(1,1,1))
dev.off()

pdf(file="figure2.pdf")
stripchart (foragingrate~sociality, vertical=TRUE, data=foraging,method="jitter",jitter=0.1,ylab="foraging rate",xlab="sociality", pch=19,cex=0.8, ylim=c(0,20))
stripchart(fitted(glm1)~foraging$sociality, vertical=TRUE,add=TRUE, pch="-",cex=2,method="jitter",col="black")
points(c(1,2),mean, pch="-",col="red",cex=5)
legend("bottomright", legend=c("fitted","observed mean"),bty="n",lwd=2,cex=0.8, col=c("black","red"), lty=c(1,1,1))
dev.off()

boxplot(fusca$foragingrate~fusca$sociality, main='Setophaga fusca',ylab='Capture rate/min',col="white",ylim=c(0,20))
stripchart(fusca$foragingrate~fusca$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,20), method="jitter", col="black")

# Calculating teh explained deviance
pseudo.R2<-(glm1$null.deviance-glm1$deviance)/glm1$null.deviance
pseudo.R2
#pseudoR2() explained deviance)=0.3489

##Visualizing the data using ggplot!!####################################################################################################
#####!!

library(ggplot2)

# For fusca
ggplot(fusca, aes(x = dayofseason, y = foragingrate, colour = sociality)) +
  geom_point() + 
  geom_smooth(method = "lm")

head(fusca)

ggplot(fusca, aes(x = dayofseason, y = foragingrate, colour = sociality)) +
  geom_point(aes(size = flocksizeind)) + 
  geom_smooth(method = "lm")

## The best fitted model 
ggplot(fusca, aes(x = dayofseason, y = foragingrate, colour = sociality)) +
  geom_point(aes(size = flocksizespecies)) + 
  geom_smooth(method = "glm",fullrange=TRUE, method.args=list(family="quasipoisson")) + 
  facet_wrap(~species, scales = "free") + 
  theme_bw()
#ggplot(fusca, aes(x = dayofseason, y = foragingrate, z= sociality, colour = sociality)) +
  geom_point(aes(size = flocksizespecies)) + 
  #geom_smooth(method= "glm", formula = foragingrate ~ sociality + dayofseason, family = quasipoisson(link = "log") + 
  #facet_wrap(~species, scales = "free") 
  
  
ggplot(fusca, aes(x = dayofseason, y = foragingrate, colour=sociality)) +
  geom_point(aes(size = flocksizespecies)) + 
  geom_smooth(method = "glm", family="quasipoisson", formula= foraging~sociality + dayofseason ) + 
  facet_wrap(~species, scales = "free") + 
theme_bw()

geom_smooth(method = "glm", formula = foragingrate~ sociality + dayofseason, method.args=list(family="quasipoisson"),data=fusca) +

plot (fusca$dayofseason,fusca$foragingrate, pch=19,
      ylab="Number of species",xlab="dayofseason", xlim=c(0,40), col=as.numeric(fusca$sociality))
lm2<-lm(foragingrate~dayofseason+sociality, data= foraging) 
abline(lm2)

plot(Numspecies~Mean_canopy_H, data=flocksd)
model<-lm(Numspecies~Mean_canopy_H, data=flocksd)
abline(model, col="red", cex=7)
summary(model)

groups <- levels(as.factor(fusca$sociality))     
for(i in 1:length(groups)){
  xi <- fusca$dayofseason[fusca$sociality==groups[i]]                  
  yhati <- fitted(glm1)[fusca$sociality==groups[i]]       
  lines(xi[order(xi)], yhati[order(xi)],col=as.numeric(i))}


### For all the species
# Smooth linear
ggplot(foraging, aes(x = dayofseason, y = foragingRate, colour = sociality)) +
  geom_point(aes(size = flocksizeind)) + 
  geom_smooth(method = "glm")+ 
  facet_wrap(~species, scales = "free") + 
  theme_bw()

#changing alfa
ggplot(foraging, aes(x = dayofseason, y = foragingRate, colour = sociality)) +
  geom_point(aes(size = flocksizeind), alpha = 0.6) + 
   geom_smooth(method = "lm") + 
  facet_wrap(~species, scales = "free") + 
  theme_bw()

### FREE SMOOTH
ggplot(foraging, aes(x = dayofseason, y = foragingRate, colour = sociality)) +
  geom_point(aes(size = flocksizeind), alpha = 0.6) + 
  # geom_smooth(method = "lm") + 
  geom_smooth() + 
  facet_wrap(~species, scales = "free") + 
  theme_bw()

# With number of species 
ggplot(foraging, aes(x = dayofseason, y = foragingrate, colour = sociality)) +
  geom_point(aes(size = flocksizespecies), alpha = 0.6) + 
  # geom_smooth(method = "lm") + 
  geom_smooth() + 
  facet_wrap(~species, scales = "free") + 
  theme_bw()

############EXTRACT THE OVERDISPERSION PARAMETER, also given when run the model, it is different, but it is recommendend the one given in the model
theta<-glmQ1$deviance/glm1$df.residual
theta

# In the analysis I incorporate the following calculation from Balker 2016, and the overdisppersion number that i gave me
#is the same that the summary of the models gave me too.
dfun<-function(glmQ1){with(glmQ1,sum((weights * residuals^2)[weights > 0])/df.residual)}
dfun(glmQ1)
summary(glmQ1)    #also give me the overdispersion parameter

#
### visualizing the  best model fit,  how to interprete this?
visreg(glmQ1, type = "conditional") #using the median for the ther predictors
visreg(glmQ1,"sociality", type = "conditional", xlim=c(0,40), ylim=c(0,30), scale = "response", ylab="",  xlab=NA) 
par(new = TRUE)
plot(foragingrate~sociality, data=fusca, ylab="foraging rate ", xlab="Sociality", main="a) Dendroica fusca ",pch=16, cex=1.3, ylim=c(0,30), xlim=c(0,40))

#visreg(glm1, type="conditional", scale="response", ylim=c(0,40), ylab="Foraging rate", pch=1,rug=0, whitespace=0.6, par=TRUE
#par(new = TRUE),
stripchart(fusca$foragingrate ~ fusca$sociality, method='jitter', vertical=TRUE, pch=19, cex=0.8, ylim=c(0,40), ylab="", whitespace=0.6)
str(fusca)

# RESULTS FOR MANUSCRIPT Hypothesis testing ------------------------------------------------------

#####################################################################
#######Testing Hypothesis USING THE BEST MODEL
# Summary of the model allow me to interpret the estimates parameters (e.g effect sizes) in the model and the difference from cero and between them.
# The summary give me the effect sizes. For example when I do summary (glm2, the model including only sociality), I can interprete the Intercept as the log transformation of  the mean of the first level of the factr sociality, 
# Inthis case the intercept is the log (mean capture rate for birds in flocks),  the socialitysolitary is the log Difference of the mean on birds in flocks with the mean of solitary. 
# And the significance p-value tell me if they are different form cero (that is not really interesting0), THE SUMMARY (model) IS MOSTLY TO SEE THE EEFFECT SIZE. From this I can say or example that birds in flock forage n given times, more when they are in flocks.
#significant difference among them. From this I can say or example that birds in flock forage n given times, more when they are in flocks.
# Be aware that interpretation is a bit more complex for the model including sociality and day of the season
# In the summary the Estimate are the stimate parameters in the log scale, 
#the first one is the mean of the first group type, the others are the difference between the first group mean and the other groups,
summary(glm2)

###################################################################################################
###################################################################################################
##########################################################################
#                Results for manuscript
###########################################################################
###################################################################################################


###################################################################################################
#1. Mixed species flocks features
###################################################################################################

#Filter only the flock data
flocks<-filter(foraging,sociality=="Flock")
flockfeatures<-na.omit(flocks)                             #omit NA values
summarise_each(flockfeatures,funs(mean))                   #Summarise all the columns
summarise(flockfeatures,mean(flocksizespecies))            #Summarise flock size 
summarise(flockfeatures,mean(flocksizeind))                #mean
summarise(flockfeatures,sd(flocksizespecies))              #standard deviation
summarise(flockfeatures,sd(flocksizeind)) 

summarise(flockfeatures,max(flocksizespecies))             #max and mins
summarise(flockfeatures,min(flocksizespecies))
summarise(flockfeatures,max(flocksizeind))
summarise(flockfeatures,min(flocksizeind))

flocksfusca<-filter(fusca,sociality=="Flock")
a<-summarise(flocksfusca,mean(foragingrate))
a

#Normality of the variables
hist(flockfeatures$flocksizespecies)
hist(flockfeatures$flocksizeind)
qqnorm(flockfeatures$flocksizeind) #normality
qqline(flockfeatures$flocksizeind, lty=2)
qqnorm(flockfeatures$flocksizespecies) #normality
qqline(flockfeatures$flocksizespecies, lty=2)

#Correlation test
#The alternative hypothesis of interest is that the flocksize is positively associated with the flockind.
cor.test(flockfeatures$flocksizespecies, flockfeatures$flocksizeind, method="kendall", alternative="greater")

###################################################################################################
################################################
#2.  Influence of social context in capture rate
# Model fit, assumptions and model selection ------------------------------
################################################
###################################################################################################

######Fit the model
### Using two predictors
###### GLM Generalized linear models ( warning: Sociality need to be a  factor)
#quasipoisson used for over disspersed count data, log link used for count data
# Because there is not AIC for quasipoisson distributions I followed the recommendtions fron (Bolker 2016)
#For migrant species I included sociality as variable, for modl selection
glm1<-glm(foragingrate~sociality+dayofseason, data =fusca,family=poisson(link="log")) 
glm2<-update(glm1, . ~ . - dayofseason)
glm3<-update(glm1, . ~ . - sociality)
glm4<-update(glm3, . ~ . - dayofseason) # Null model same than glm9<-glm(foragingrate~1, data =fusca, family=poisson(link="log"))
glmQ1<-update(glm1, family=quasipoisson(link="log"))
glmQ2<-update(glm2, family=quasipoisson(link="log"))
glmQ3<-update(glm3, family=quasipoisson(link="log"))
glmQ4<-update(glm4, family=quasipoisson(link="log"))
warnings()
#glm5<-glm(foragingrate~1, data =fusca, family=poisson(link="log")) Null model
summary(glmQ1)
summary(glm3)
summary(glm4)
summary(glm9)
# Checking the assumptions of the model################################################
#1 Residuals are normally distributed
fv<-fitted(glm1)    # predicted (fitted) values on original scale
re<-residuals(glm1,type="response")  # residuals
plot(fv,re)     # For constant variance, if homogeneous distrbuted in the graph
qqnorm(re) # to check the normality of the data
qqline(re) # to add the line of normalite
hist(re)         # to see if the residuals distribute normaly

###2 the residuals are  homogeneus  (homosedatic), if they are disperse equal acroos the graph ###oops they are not 
plot(re~glm1$fitted.values)

#then, calculate theta to see if there are overdisoersion, if >1 is overdispersion, the use a quiasipoisson, 
theta<-glm1$deviance/glm1$df.residual
theta
#theta=1.98

#3### not autocorrelation ( correlation between succesive data points)
durbinWatsonTest(glm1)

#Results of the test, if pvalue less than 0.05  there is a significan serial autocorrelation
#durbinWatsonTest(gamdiversity1)
#lag Autocorrelation D-W Statistic p-value
#1      0.09337523       1.80187    0.16
#Alternative hypothesis: rho != 0

#4# the model is not biased by influencial observations
#leverage #not sure how to calculate it
#172 is an influencial point
# removing the influencial point
flocksd$Number_of_species[172]
gamdiversity1updated<-update(gamdiversity1,subset=(flocksd$Number_of_species!=6))
##### ploting the new model with out the influencial point
#plot(gamdiversity1updated, pch=19, cex=0.25, col='#FF8000',shade=TRUE, shade.col="grey",se=TRUE)
#plot(gamdiversity1updated, residuals=TRUE, pch=19, cex=0.6, col='navy',shade=TRUE, shade.col="grey",se=TRUE)
#summary(gamdiversity1) 
#summary(gamdiversity1updated) 

################################################################
# Model selection ---------------------------------------------------------
######################################################################
# With two predictors
glm1<-glm(foragingrate~sociality+dayofseason, data =fusca,family=poisson(link="log"))
glm1<-glm(foragingrate~sociality+dayofseason, data =peregrina,family=poisson(link="log"))
glm1<-glm(foragingrate~sociality+dayofseason, data =canadensis,family=poisson(link="log"))
glm1<-glm(foragingrate~sociality+dayofseason, data =cerulea,family=poisson(link="log"))

glm1<-glm(foragingrate~sociality+dayofseason, data =chrysops,family=poisson(link="log")) 
glm1<-glm(foragingrate~sociality+dayofseason, data =guira,family=poisson(link="log"))
glm1<-glm(foragingrate~sociality+dayofseason, data =pitiayumi,family=poisson(link="log"))

glm2<-update(glm1, . ~ . - dayofseason) #social context
glm3<-update(glm1, . ~ . - sociality) #day of season
glm4<-update(glm3, . ~ . - dayofseason)
glmQ1<-update(glm1, family=quasipoisson(link="log"))
glmQ2<-update(glm2, family=quasipoisson(link="log"))
glmQ3<-update(glm3, family=quasipoisson(link="log"))
glmQ4<-update(glm4, family=quasipoisson(link="log"))

warnings()

glm5<-glm(foragingrate~1, data =fusca, family=poisson(link="log")) #Null model equivalent to glm2 long version

#Interesting models
glm6<-glm(foragingrate~sociality+dayofseason-1, data=fusca, family=poisson(link="log"))# Give me the log estimates of each group directly
glm7<-update(glm1, . ~ . + dayofseason:sociality) # including an interaction in the model

# Summary of the model allow me to interpret the estimates of the parameters (e.g effect sizes) in the model and the difference from cero and between them.
# The summary give me the effect size. 
summary(glm2)
summary (glm6)
summary(glmQ1)

#####overdisperssion parameters 
dfun<-function(glmQ1q){with(glmQ1,sum((weights * residuals^2)[weights > 0])/df.residual)}
dfun(glmQ1)

summary(glmQ1)

#fusca=1.97
#cerulea=2.22
#peregrina=2.27
#canadensis=1.56
#pitiayumi=1.67
#guira=2.92
#chrysops= na


#### Model selection using  AICcmodavg
library(AICcmodavg)

###Model used  (glm1Q,glm2Q,glm3Q,glm4Q), for all the species which have overdispersed datta, remember to use the glm models instead of glmQ
aictab(list(glm1,glm2,glm3,glm4),
       modnames=c("Socialcontext+Dayofseason",
                  "Socialcontext",
                  "Dayofseason",
                  "Intercept"),
       c.hat=1.97)

summary (glmQ2)

#### To make the table for AIC cfor species that are not overdisperssed

aictab(list(glm1,glm2,glm3,glm4),
       modnames=c("Socialcontext+Dayofseason",
                  "Socialcontext",
                  "Dayofseason",
                  "Intercept"))


#4.   ######SUMMARY TABLE
#A beautiful table in html for the parameters of glm models
#For the parameter stimates we need to consider the models with the Quasipoisson distribution, the points estimates are identical
# to the model with poisson distribution but the standard error and confidence intervals are wider!
sjt.glm(glmQ1,glmQ2,glmQ3,
        depvar.labels = c("Model1: Socialcontext+Dayofseason","Model2:Socialcontext","Model3:Day of season"),
        pred.labels = NULL,
        show.aicc=TRUE, 
        show.family=TRUE, 
        group.pred = FALSE,
        exp.coef = TRUE,  # if true Regression and cof.intervals are exponentiaded st.error are not in the unstransformed scale
        p.numeric=TRUE,
        robust = TRUE,
        show.se = TRUE,
        show.r2=TRUE,
        show.dev=TRUE,
        show.chi2=TRUE,
        cell.spacing = 0.001,
        sep.column = FALSE)

# A pvalue less than 0.05 indicates a good model fit ?
# Summary of any model 
summary (glm2)
###To underestand the summary(glm) and the table with the estimate for sethophaga fusca.
#Let's start form the model that only include sociality as the only predictor (Flock and Solitaryare the leves). In this one the Intercept is the mean of the first level of the factor, in flocks, but it is in the log scale mean 5.25 log(5.25 =1.66) 
# next predictor is sociality(solitary) and is calculated as the difference of the logs of flock and solitary. log flock mean(5.25)log - LOG solitary mean =(1.65). In other words 1.65-0.5=1.15. If we exponentiate that value exp(1.15) that means that a fusca individual in flocks
# will eat 3.15 times more than solitary * 0.03 times the day of the season. Or what is the same a  solitary individual of fusca will eat exp(-1.33)= 0.26 times less what it will eat when in flocks~ 26% less

#The model  will be for fusca f.rate= exp [1.27]* exp[-1.33]*exp[0.03 * #day of the season]
# To extract the parameters manually
glmQ1<-glm(foragingrate~sociality+dayofseason, data =fusca,family=poisson(link="log"))
summary(glmQ1)
coef(glmQ1)
exp(coef(glmQ1))
exp(confint.default(glmQ1)) 
# to calculate confidence intervals manually
library(MASS)
exp (confint(glmQ1, level = 0.95))


#Hypothesis testing using ANOVA, are the difference in foraging rate statistically significant? USING
glm1<-glm(foragingrate~sociality+dayofseason, data =fusca,family=poisson(link="log")) 
glmQ1<-update(glm1, family=quasipoisson(link="log"))

### Anova allow me to do hyphothesis testing
# The ANOVA atable compares the fit of two models, the null model foraging~1 vs foraging~sociality+....
## f test appropiate for quasipoisson distributions, and type3 anova order orthe term does not matter
# test the null H that there is no differences in the foraging rate  among flocking and non-flocking individuals
#Note: Use anova(test = "F"), ie an F test, if testing hypotheses that use gaussian, quasibinomial or quasipoisson link functions.
#This is because the quasi-likelihood is not a real likelihood and the generalized log-likelihood ratio test is not accurate FROM dolph class
#Anova type=3 rom the car package, in that case order of the terms in the model does not matter and neither does hirarchy

#anova(model2, test="Chi") using the "best model" I guess i can use any of the two qasi or poisson since the stimates are the same?
### For the migrants 
glm1<-glm(foragingrate~sociality+dayofseason, data =fusca,family=poisson(link="log"))
glmQ1<-update(glm1, family=quasipoisson(link="log"))
Anova(glmQ1, type = 3, test="F")

glm1<-glm(foragingrate~sociality+dayofseason, data =peregrina,family=poisson(link="log"))
glmQ1<-update(glm1, family=quasipoisson(link="log"))
Anova(glmQ1, type = 3, test="F")

glm1<-glm(foragingrate~sociality+dayofseason, data =canadensis,family=poisson(link="log"))
glmQ1<-update(glm1, family=quasipoisson(link="log"))
Anova(glmQ1, type = 3, test="F")

glm1<-glm(foragingrate~sociality+dayofseason, data =cerulea,family=poisson(link="log"))
glmQ1<-update(glm1, family=quasipoisson(link="log"))
Anova(glmQ1, type = 3, test="F")

#For the residents
glm1<-glm(foragingrate~sociality+dayofseason, data =chrysops,family=poisson(link="log")) 
glm2<-update(glm1, . ~ . - dayofseason) #social context
Anova(glmQ2, type = 3, test="F") 

glm1<-glm(foragingrate~sociality+dayofseason, data =pitiayumi,family=poisson(link="log")) 
glm2<-update(glm1, . ~ . - dayofseason) #social context
glmQ2<-update(glm2, family=quasipoisson(link="log"))
Anova(glmQ2, type = 3, test="F") 

glm1<-glm(foragingrate~sociality+dayofseason, data =guira,family=poisson(link="log")) 
glm2<-update(glm1, . ~ . - dayofseason) #social context
glmQ2<-update(glm2, family=quasipoisson(link="log"))
Anova(glmQ2, type = 3, test="F") 


#Other way to do this no te the "a"
anova (glm1,glm5, test="F") # BEING glm5 the null model 

## f test appropiate for quasipoisson distributions, and type3 anova order orthe term does not matter

# To test a H, using Anova
Anova(glmmodel2, type = 3, test="F") 


#Examining the fit of the best model for each species graphically
#canadensis 
glm1<-glm(foragingrate~sociality+dayofseason, data =canadensis,family=poisson(link="log"))
glm2<-update(glm1, . ~ . - dayofseason) #social context
glmQ2<-update(glm2, family=quasipoisson(link="log"))
plot(glmQ2)
?visreg()
stripchart(canadensis$foragingrate~ canadensis$sociality, method='jitter', vertical=TRUE, pch=19, cex=0.8, ylim=c(0,50), ylab="",whitespace=0.6)
par(new = TRUE)
visreg(glmQ2, type="conditional", scale="response", ylim=c(0,50), ylab="Capture rate", pch=1,rug=0, whitespace=0.6)

#cerulea
glm1<-glm(foragingrate~sociality+dayofseason, data =cerulea,family=poisson(link="log"))
glmQ1<-update(glm1, family=quasipoisson(link="log"))
plot(glmQ1)
?visreg()
stripchart(cerulea$foragingrate~ cerulea$sociality, method='jitter', vertical=TRUE, pch=19, cex=0.8, ylim=c(0,50), ylab="",whitespace=0.6)
par(new = TRUE)
visreg(glmQ1, type="conditional", scale="response", ylim=c(0,50), ylab="Capture rate", pch=1,rug=0, whitespace=0.6)

#fusca
glm1<-glm(foragingrate~sociality+dayofseason, data =fusca,family=poisson(link="log"))
glmQ1<-update(glm1, family=quasipoisson(link="log"))
plot(glmQ1)
?visreg()
stripchart(fusca$foragingrate~ fusca$sociality, method='jitter', vertical=TRUE, pch=19, cex=0.8, ylim=c(0,50), ylab="",whitespace=0.6)
par(new = TRUE)
visreg(glmQ1, type="conditional", scale="response", ylim=c(0,50), ylab="Capture rate", pch=1,rug=0, whitespace=0.6)


#peregrina
glm1<-glm(foragingrate~sociality+dayofseason, data =peregrina,family=poisson(link="log"))
glmQ1<-update(glm1, family=quasipoisson(link="log"))
plot(glmQ1)
?visreg()
stripchart(peregrina$foragingrate~ peregrina$sociality, method='jitter', vertical=TRUE, pch=19, cex=0.8, ylim=c(0,50), ylab="",whitespace=0.6)
par(new = TRUE)
visreg(glmQ1, type="conditional", scale="response", ylim=c(0,50), ylab="Capture rate", pch=1,rug=0, whitespace=0.6)


###################################
# With three predictors not used in the first versin of the paper
###################################
glm1a<-glm(foragingrate~sociality+dayofseason+flocksizespecies, data =fusca,family=poisson(link="log"))
glm1a<-glm(foragingrate~sociality+dayofseason+flocksizespecies, data =peregrina,family=poisson(link="log"))
glm1a<-glm(foragingrate~sociality+dayofseason+flocksizespecies, data =canadensis,family=poisson(link="log"))
glm1a<-glm(foragingrate~sociality+dayofseason+flocksizespecies, data =cerulea,family=poisson(link="log"))

glm1a<-glm(foragingrate~sociality+dayofseason+flocksizespecies, data =chrysops,family=poisson(link="log"))
glm1a<-glm(foragingrate~sociality+dayofseason+flocksizespecies, data =guira,family=poisson(link="log"))
glm1a<-glm(foragingrate~sociality+dayofseason+flocksizespecies, data =pitiayumi,family=poisson(link="log"))

glm2a<-update(glm1a, . ~ . - sociality) # dayofseason + flocksize
glm3a<-update(glm1a, . ~ . - dayofseason) #sociality + flocksize
glm4a<-update(glm1a, . ~ . - flocksizespecies) # sociality+dayofseason
glm5a<-update(glm2a, . ~ . - dayofseason) # flocksize
glm6a<-update(glm2a, . ~ . - flocksizespecies) #dayofseason
glm7a<-update(glm4a, . ~ . - dayofseason )  #sociality
glm8a<-update(glm7a, . ~ . - sociality)  #intercept Equivalent to null model glm(foragingrate~1, data =fusca, family=poisson(link="log"))
glmQ1a<-update(glm1a, family=quasipoisson(link="log"))
glmQ2a<-update(glm2a, family=quasipoisson(link="log"))
glmQ3a<-update(glm3a, family=quasipoisson(link="log"))
glmQ4a<-update(glm4a, family=quasipoisson(link="log"))
glmQ5a<-update(glm5a, family=quasipoisson(link="log"))
glmQ6a<-update(glm6a, family=quasipoisson(link="log"))
glmQ7a<-update(glm7a, family=quasipoisson(link="log"))
glmQ8a<-update(glm8a, family=quasipoisson(link="log"))

############EXTRACT THE OVERDISPERSION PARAMETER, also given when run the model, it is different, but it is recommendend the one given in the model
theta<-glm1a$deviance/glm1a$df.residual
theta
# In the analysis I incorporate the following calculation from Balker 2016, and the overdisppersion number that i gave me
#is the same that the summary of the models gave me too.
dfun<-function(glmQ1q){with(glmQ1a,sum((weights * residuals^2)[weights > 0])/df.residual)}
dfun(glmQ1a)
summary(glmQ1a)    #also give me the overdispersion parameter

#####overdisperssion parameters 
#fusca=1.99
#cerulea=2.36
#peregrina=2.33
#canadensis=1.54
#pitiayumi=1.73
#guira=3.08
#chrysops= na data not overdispersed
###############################################################
##### Model selection using  AICcmodavg
library(AICcmodavg)
###Model used in the thesis
aictab(list(glm1a,glm2a,glm3a,glm4a,glm5a,glm6a,glm7a,glm8a),
       modnames=c("Socialcontext+Dayofseason+Flocksize",
                  "Dayofseason+Flocksize",
                  "Socialcontext+Flocksize",
                  "Socialcontext+dayofseason",
                  "Flocksize",
                  "Dayofseason",
                  "Socialcontext",
                  "Intercept"),
       c.hat=3.08) # need to change for the overdispersion for each species
summary (glmQ5a)
# Influence of sociality in the movement rate ##################################################################################################
###################################################################################################
# Influence of sociality in the movement rate 
###################################################################################################
###################################################################################################

# With two predictors
#Migrants
glm1m<-glm(movementrate~sociality+dayofseason, data =fusca,family=poisson(link="log")) #full model fusca
glm1m<-glm(movementrate~sociality+dayofseason, data =peregrina,family=poisson(link="log"))
glm1m<-glm(movementrate~sociality+dayofseason, data =canadensis,family=poisson(link="log"))
glm1m<-glm(movementrate~sociality+dayofseason, data =cerulea,family=poisson(link="log"))
#Resident
glm1m<-glm(movementrate~sociality+dayofseason, data =chrysops,family=poisson(link="log")) 
glm1m<-glm(movementrate~sociality+dayofseason, data =guira,family=poisson(link="log"))
glm1m<-glm(movementrate~sociality+dayofseason, data =pitiayumi,family=poisson(link="log"))

glm2m<-update(glm1m, . ~ . - dayofseason) #social context
glm3m<-update(glm1m, . ~ . - sociality) #day of season
glm4m<-update(glm3m, . ~ . - dayofseason) #Null or intercept model
glmQ1m<-update(glm1m, family=quasipoisson(link="log"))
glmQ2m<-update(glm2m, family=quasipoisson(link="log"))
glmQ3m<-update(glm3m, family=quasipoisson(link="log"))
glmQ4m<-update(glm4m, family=quasipoisson(link="log"))

#glm5<-glm(foragingrate~1, data =fusca, family=poisson(link="log")) #Null model equivalent to glm2 long version

# To extract summary and coefficients of the data
summary(glmQ1m)
coef(glmQ1m)
exp(coef(glmQ1m))
summary(glmQ2m)
summary(glmQ3m)
summary(glmQ4m)

#####overdisperssion parameters 
dfun<-function(glmQ1qm){with(glmQ1m,sum((weights * residuals^2)[weights > 0])/df.residual)}
dfun(glmQ1m)

#fusca=1.77
#cerulea=2.51
#peregrina=1.89
#canadensis=2.96
#pitiayumi= 2.42
#guira=2.42
#chrysops= 2.44

#### Model selection using  AICcmodavg
library(AICcmodavg)

###Model used  (glm1Qm,glm2Qm,glm3Qm,glm4Qm), for all the species which have overdispersed datta

aictab(list(glm1m,glm2m,glm3m,glm4m),
       modnames = c("Socialcontext+Dayofseason",
                    "Social context",
                    "Day of season",
                    "Intercept"),
       c.hat=2.42)




#### To make the table for AIC cfor species that are not overdisperssed

aictab(list(glm1,glm2,glm3,glm4),
       modnames=c("Socialcontext+Dayofseason",
                  "Socialcontext",
                  "Dayofseason",
                  "Intercept"))


#4.   ######SUMMARY TABLE
#A beautiful table in html for the parameters of glm models
#For the parameter stimates we need to consider the models with the Quasipoisson distribution, the points estimates are identical
# to the model with poisson distribution but the standard error and confidence intervals are wider!
sjt.glm(glmQ1m,glmQ2m,glmQ3m,
        depvar.labels = c("Model1: Socialcontext+Dayofseason","Model2:Socialcontext","Model3:Day of season"),
        pred.labels = NULL,
        show.aicc=TRUE, 
        show.family=TRUE, 
        group.pred = FALSE,
        exp.coef = TRUE,  # if true Regression and cof.intervals are exponentiaded st.error are not in the unstransformed scale
        p.numeric=TRUE,
        robust = TRUE,
        show.se = TRUE,
        show.r2=TRUE,
        show.dev=TRUE,
        show.chi2=TRUE,
        cell.spacing = 0.001,
        sep.column = FALSE)


#Hypothesis testing using Anova
#fusca
glm1m<-glm(movementrate~sociality+dayofseason, data =fusca,family=poisson(link="log")) #full model fusca
glmQ1m<-update(glm1m, family=quasipoisson(link="log"))
glm2m<-update(glm1m, . ~ . - dayofseason) #social context
glm3m<-update(glm1m, . ~ . - sociality) #day of season
glmQ3m<-update(glm3m, family=quasipoisson(link="log"))

Anova(glmQ1m, type = 3, test="F")

#peregrina
glm1m<-glm(movementrate~sociality+dayofseason, data =peregrina,family=poisson(link="log")) #full model fusca
glmQ1m<-update(glm1m, family=quasipoisson(link="log"))
summary(glmQ1m)
Anova(glmQ1m, type = 3, test="F")

#canadensis
glm1m<-glm(movementrate~sociality+dayofseason, data =canadensis,family=poisson(link="log")) #full model fusca
glmQ1m<-update(glm1m, family=quasipoisson(link="log"))
glm2m<-update(glm1m, . ~ . - dayofseason) #social context
glm3m<-update(glm1m, . ~ . - sociality) #day of season
glmQ3m<-update(glm3m, family=quasipoisson(link="log"))

Anova(glmQ1m, type = 3, test="F")

#cerulea the best model is the null model!!!

#chrysops
glm1m<-glm(movementrate~sociality+dayofseason, data =chrysops,family=poisson(link="log")) #full model fusca
glmQ1m<-update(glm1m, family=quasipoisson(link="log"))
summary(glmQ1m)
Anova(glmQ1m, type = 3, test="F")

#parula
glm1m<-glm(movementrate~sociality+dayofseason, data =pitiayumi,family=poisson(link="log"))
glm2m<-update(glm1m, . ~ . - dayofseason) #social context
glmQ2m<-update(glm2m, family=quasipoisson(link="log"))

Anova(glmQ2m, type = 3, test="F")


# With three predictors not used in the first versin of the paper
###################################

glm1a<-glm(movementrate~sociality+dayofseason+flocksizespecies, data =fusca,family=poisson(link="log"))
glm1a<-glm(movementrate~sociality+dayofseason+flocksizespecies, data =peregrina,family=poisson(link="log"))
glm1a<-glm(movementrate~sociality+dayofseason+flocksizespecies, data =canadensis,family=poisson(link="log"))
glm1a<-glm(movementrate~sociality+dayofseason+flocksizespecies, data =cerulea,family=poisson(link="log"))

glm1a<-glm(movementrate~sociality+dayofseason+flocksizespecies, data =chrysops,family=poisson(link="log"))
glm1a<-glm(movementrate~sociality+dayofseason+flocksizespecies, data =guira,family=poisson(link="log"))
glm1a<-glm(movementrate~sociality+dayofseason+flocksizespecies, data =pitiayumi,family=poisson(link="log"))

glm2a<-update(glm1a, . ~ . - sociality) # dayofseason + flocksize
glm3a<-update(glm1a, . ~ . - dayofseason) #sociality + flocksize
glm4a<-update(glm1a, . ~ . - flocksizespecies) # sociality+dayofseason
glm5a<-update(glm2a, . ~ . - dayofseason) # flocksize
glm6a<-update(glm2a, . ~ . - flocksizespecies) #dayofseason
glm7a<-update(glm4a, . ~ . - dayofseason )  #sociality
glm8a<-update(glm7a, . ~ . - sociality)  #intercept Equivalent to null model glm(movementrate~1, data =fusca, family=poisson(link="log"))
glmQ1a<-update(glm1a, family=quasipoisson(link="log"))
glmQ2a<-update(glm2a, family=quasipoisson(link="log"))
glmQ3a<-update(glm3a, family=quasipoisson(link="log"))
glmQ4a<-update(glm4a, family=quasipoisson(link="log"))
glmQ5a<-update(glm5a, family=quasipoisson(link="log"))
glmQ6a<-update(glm6a, family=quasipoisson(link="log"))
glmQ7a<-update(glm7a, family=quasipoisson(link="log"))
glmQ8a<-update(glm8a, family=quasipoisson(link="log"))


aictab(list(glm1a,glm2a,glm3a,glm4a,glm5a,glm6a,glm7a,glm8a),
       modnames=c("Socialcontext+Dayofseason+Flocksize",
                  "Dayofseason+Flocksize",
                  "Socialcontext+Flocksize",
                  "Socialcontext+dayofseason",
                  "Flocksize",
                  "Dayofseason",
                  "Socialcontext",
                  "Intercept"),
       c.hat=1.92)

summary(glmQ1a)

# Foraging maneuvers##################################################################################################
###################################################################################################
# Foraging maneuvers
###################################################################################################
##############################################################################################

###Load Lattice ###
library(lattice)
library(dplyr)
library(dev2bitmap)
#reading the data
maneuvers<-read.csv(file.choose("maneuver"),stringsAsFactors=FALSE,strip.white=TRUE,na.strings=c("NA",""))

#filtering data

#fuscamaneuver<-filter(maneuvers,species=="Setophaga fusca")
#parulamaneuver<-filter(maneuvers,species=="Parula pitiayumi")

#variables as factors
maneuvers$context<-as.factor(maneuvers$context)

### Barchart for species, I can vary the number of columns and rows, and by using the driver command with dimensions, I can make the graph to fit my desire lenght###
barchart(proportion~maneuver|species,data=maneuvers, layout=c(3,2),scale=list(x=list(rot=0,cex=0.7),y=list(cex=0.7)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,1),box.ratio=1.5,col="gray")

# using context as a grouping category
barchart(proportion~maneuver|species, groups=context, data=maneuvers, layout=c(3,2),scale=list(x=list(rot=0,cex=0.7),y=list(cex=0.7)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,1),box.ratio=1.5, auto.key=TRUE)
# in gay colors
barchart(proportion~maneuver|species, groups=context, data=maneuvers, layout=c(3,2),scale=list(x=list(rot=0,cex=0.7),y=list(cex=0.7)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,1),box.ratio=1.5, col=c("grey","white"))
barchart(percentage~maneuver|species, groups=context, data=maneuvers, layout=c(3,2),scale=list(x=list(rot=0,cex=0.7),y=list(cex=0.7)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,100),box.ratio=1.5, col=c("grey","white"))

pdf(file="Figure maneuvers.pdf",width=8,height=10)
barchart(percentage~maneuver|species, groups=context, data=maneuvers, layout=c(3,2),scale=list(x=list(rot=0,cex=0.7),y=list(cex=0.7)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,100),box.ratio=1.5, col=c("grey","white"))
dev.off()

tiff("Plotmaneuvers.tiff", res = 300, width =1000 , height =1015 )
barchart(proportion~maneuver|species, groups=context, data=maneuvers, layout=c(3,2),scale=list(x=list(rot=0,cex=0.7),y=list(cex=0.7)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,1),box.ratio=1.5, col=c("grey","white"))
dev.off()

####################################
#Substrate
####################################

library(lattice)

#reading the data
substrate<-read.csv(file.choose("substrate"),stringsAsFactors=FALSE,strip.white=TRUE,na.strings=c("NA",""))

#variables as factors
substrate$context<-as.factor(substrate$context)
str(substrate)

# in gay colors
barchart(proportion~substrate|species, groups=context, data=substrate, layout=c(3,2),scale=list(x=list(rot=0,cex=0.7),y=list(cex=0.7)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,1),box.ratio=1.5, auto.key=TRUE)
barchart(proportion~substrate|species, groups=context, data=substrate, layout=c(3,2),scale=list(x=list(rot=0,cex=0.7),y=list(cex=0.7)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,1),box.ratio=1.5, col=c("grey","white"))
barchart(percentage~substrate|species, groups=context, data=substrate, layout=c(3,2),scale=list(x=list(rot=0,cex=0.7),y=list(cex=0.7)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,100),box.ratio=1.5, col=c("grey","white"))

pdf(file="Figure substrate2.pdf",width=12,height=8)
barchart(percentage~substrate|species, groups=context, data=substrate, layout=c(3,2),scale=list(x=list(rot=0,cex=0.7),y=list(cex=0.7)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,100),box.ratio=1.5, col=c("grey","white"))
dev.off()

#Note remember to order by aerial and non aerial for easier analyises

# RESULTS F test if p<0.005 differences between groups are significant
#Analysis of Deviance Table (Type III tests)

Response: Similarity
#SS    Df      F Pr(>F)
#Type.flock 0.17824  2 2.1055 0.1342
#Residuals  1.82006 43 

#### years
boxplot(stability2$Similarity ~ stability2$Type.flock, ylab="Stability over years (1-Jaccard similarity index)",xlab="Elevation", xlab='Type of flock', col="light blue",, ylim=c(0,1))
stripchart(stability2$Similarity ~ stability2$Type.flock, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,1))

#### days small sample size
boxplot(stability1$Similarity ~ stability1$Type.flock, ylab='Flock similarity over time', xlab='Type of flock', col="light blue",, ylim=c(0,1))
stripchart(stability1$Similarity ~ stability1$Type.flock, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,1))

#Linearfixed model
lm<-lm(Foraging.rate ~ Flocksize+Gender+ Flocksize:Gender, data=foraging)
plot(lm)

#Transforming the response variable to met assumptions
foraging$Logforaging <- log((foraging$Foraging.rate)+1)

# visualizing the data
plot(foraging$Flocksize,foraging$Logforaging, 
     pch = as.numeric(foraging$Gender), col = c(foraging$Gender))


#Lineal fixed model including  gender, and the interaction
linear2<-lm(Logforaging ~ Gender+Flocksize + Flocksize:Gender, data=foraging)
plot(linear2)
#Lineal fixed model including the factor gender
linear1<-lm(Logforaging ~ Gender+Flocksize, data=foraging)
#Lineal simple fixed model
linear<-lm(Logforaging ~ Flocksize, data=foraging)

#Exploring the explanatory power of both models
# Anova type =2 where hierarchy matter but order do not.

library(car)

Anova(linear2,type=2)

#The parameter estimates for the model including the interaction

summary (linear2)

# But after the Anova I knew that
#the simpler model (linear) does not have signicicantly lower explanatoty power
summary (linear)

#Plot the model including both factors and the interaction term.
plot (foraging$Flocksize,foraging$Logforaging, 
      pch = as.numeric(foraging$Gender), col = c(foraging$Gender),
      ylab="Foraging rate (cap/min)",xlab="Flock size (Number of species)")

groups <- levels(foraging$Gender)                     
for(i in 1:length(groups)){
  xi <- foraging$Flocksize[foraging$Gender==groups[i]]                  
  yhati <- fitted(linear2)[foraging$Gender==groups[i]]       
  lines(xi[order(xi)], yhati[order(xi)],col=as.numeric(i)) }


# Plot the simpler model

plot (foraging$Flocksize,foraging$Logforaging, 
      pch = as.numeric(foraging$Gender), col = c(foraging$Gender),
      ylab="Foraging rate (cap/min)",xlab="Flock size (Number of species)")

groups <- levels(foraging$Gender)                     
for(i in 1:length(groups)){
  xi <- foraging$Flocksize[foraging$Gender==groups[i]]                  
  yhati <- fitted(linear)[foraging$Gender==groups[i]]       
  lines(xi[order(xi)], yhati[order(xi)],col=as.numeric(i)) }

### plots foraging rate and flock size
plot(fusca$flocksizespecies,fusca$foragingrate, 
     pch = as.numeric(fusca$gender), col = c(fusca$gender))

plot(fusca$foragingrate,fusca$flocksizespecies)

#####For movement rate

## boxplot with stripchart on the back for foraging data
#Migrants
#Setophaga fusca
boxplot(fusca$movementrate~fusca$sociality, main='fusca',ylab='Capture rate/min',col="white",ylim=c(0,20))
stripchart(fusca$movementrate~fusca$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,20))
#Setophaga cerulea
boxplot(cerulea$movementrate~cerulea$sociality, main='cerulea',ylab='Capture rate/min',col="white",ylim=c(0,20))
stripchart(cerulea$movementrate~cerulea$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,20))
#Cardelina canadensis
boxplot(canadensis$movementrate~canadensis$sociality, main='canadensis',ylab='Capture rate/min',col="white",ylim=c(0,15))
stripchart(canadensis$movementrate~canadensis$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,15))
#Oreothlypis peregrina
boxplot(peregrina$movementrate~peregrina$sociality, main='peregrina',ylab='Capture rate/min',col="white",ylim=c(0,25))
stripchart(peregrina$movementrate~peregrina$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,25))
#Residents
#Hemithraupis guira
boxplot(guira$movementrate~guira$sociality, main='guira',ylab='Capture rate/min',col="white",ylim=c(0,15))
stripchart(guira$movementrate~guira$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,15))
#Zimmerius chrysops
boxplot(chrysops$movementrate~chrysops$sociality, main='chrysops',ylab='Capture rate/min',col="white",ylim=c(0,10))
stripchart(chrysops$movementrate~chrysops$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,10))
#Parula pitiayumi
boxplot(pitiayumi$movementrate~pitiayumi$sociality, main='pitiayumi',ylab='Capture rate/min',col="white",ylim=c(0,20))
stripchart(pitiayumi$movementrate~pitiayumi$sociality, method='jitter', add=TRUE, vertical=TRUE, pch=19, cex=0.8, ylim=c(0,20))

# Foraging rate vs group size
# Combining the plots in a multiplot graph, be aware of the linear model

ggplot(foraging, aes(flocksizespecies, foragingRate)) + xlab("X variable name") + ylab("Y variable name") +
  geom_point(col = "black", size = I(2)) +
  geom_smooth(method = glm, size = I(1), se = FALSE, col = "black") +
  facet_wrap(~species, ncol = 0)

Models
fuscaflocks<-read.csv(file.choose("fusca"), stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))

glmfusca<-glm(foragingrate~flocksizespecies, data =fuscaflocks,family=poisson(link="log"))
glmfusca1<-glm(foragingrate ~ poly(flocksizespecies,2), family=poisson(link="log"), data=fuscaflocks) 


aictab(list(glmfusca,glmfusca1),
       modnames=c("glm",
                  "polynomial"))
####

#note neeed to exclude solitary individuals
ggplot(flocks, aes(flocksizeind, foragingRate)) + xlab("X variable name") + ylab("Y variable name") +
  geom_point(col = "black", size = I(2)) +
  geom_smooth(method = glm, size = I(1), se = FALSE, col = "black") +
  facet_wrap(~species, ncol = 0)

ggplot(foraging, aes(flocksizeind, movementrate)) + xlab("X variable name") + ylab("Y variable name") +
  geom_point(col = "black", size = I(2)) +
  geom_smooth(method = glm, size = I(1), se = FALSE, col = "black") +
  facet_wrap(~species, ncol = 2)

ggplot(flock, aes(flocksizeind, foragingrate)) + xlab("X variable name") + ylab("Y variable name") +
  geom_point(col = "black", size = I(2)) +
  geom_smooth(method = glm, size = I(1), se = FALSE, col = "black") +
  facet_wrap(~species, ncol = 2)

ggplot(flock, aes(flocksizeind, foragingrate)) + xlab("X variable name") + ylab("Y variable name") +
  geom_point(col = "black", size = I(2)) +
  geom_smooth( size = I(1), se = FALSE, col = "black") +
  facet_wrap(~species, ncol = 2)

# graph
ggplot(data=foraging) +
  geom_point(mapping= aes(x=flocksizeind, y=foragingRate, col=species))+
  geom_smooth(mapping = aes(x=flocksizeind, y=foragingRate, linetype=species))
geom_point(mapping= aes(x=flocksizeind, y=movementrate))

ggplot(flocks, aes(flocksizeind, movementrate)) +
  geom_point( col="black", size=I(2)) +
  geom_smooth(method = glm, size=I(1), se=FALSE, col="black")+
  geom_point(mapping= aes(x=flocksizeind, y=movementrate))

ffusca<-filter(flock,species=="Setophaga fusca")

ggplot(data=fusca) +
  geom_point(mapping= aes(x=flocksizeind, y=foragingrate), col="blue")+
  stat_smooth(mapping = aes(x=flocksizeind, y=foragingrate), col="blue",method = "glm", formula = y~ (x^2))

facet_wrap(~species,ncol = 2)
+
  geom_point(mapping= aes(x=flocksizeind, y=movementrate), col="black")+
  geom_smooth(mapping = aes(x=flocksizeind, y=movementrate), col="black")


####Using lattice
barchart(proportion~substrate|species, groups=context, data=substrate, layout=c(3,2),scale=list(x=list(rot=0,cex=0.7),y=list(cex=0.7)),ylab="Relative proportions",xlab="Foraging maneuver",ylim=c(0,1),box.ratio=1.5, auto.key=TRUE)
xyplot(foragingrate~flocksizespecies|species, data=flock, cex=0.5)

#Using visreg
library(visreg)
library(mgcv)''
glmfusca<-glm(foragingRate~flocksizespecies, data =fusca2,family=quasipoisson(link="log"))
glm(y~x)

glmfusca1<-glm(foragingRate ~ poly(flocksizespecies,2), family=quasipoisson(link="log"), data=fusca2) 
glmfusca<-glm(flocksizespecies ~ poly(foragingrate,2), data=ffusca) 

summary(glmfusca1)
#fitting a polynomial function

visreg(glmfusca,"foragingrate",xlab="Flock size",ylab="Foragingrate",type = "conditional", scale="response",ylim=c(0,20),xlim=c(0,40),shade=NULL)
par(new = TRUE)
plot(foragingrate~flocksizespecies, data=ffusca, ylab="Foraging rate", xlab="Flock size (Number of species)",xlim=c(0,40), ylim=c(0,20))

aictab(list(glmfusca,glmfusca1),
       modnames=c("glm",
                  "polynomial"))

# for the number of individuals

visreg(glmfusca,"flocksizespecies", type = "conditional", scale = "response", xlim=c(0,40), ylim=c(0,20), ylab=NA, xlab=NA)
par(new = TRUE)
plot(foragingrate~flocksizespecies, data=ffusca, ylab="Flock size (Number of species)", xlab="Foraging rate",xlim=c(0,40), ylim=c(0,20))

gamfusca<-gam(foragingrate~s(flocksizespecies), data=ffusca, family=quasipoisson,link=("log")) 
plot(gamfusca, residuals=FALSE, pch=19, cex=0.6, col='navy',shade=TRUE, shade.col="grey",se=TRUE)
vis.gam(gamfusca, type = "response", plot.type = "contour")
visreg(gamfusca,"flocksizespecies",xlab="Flock size",ylab="Foragingrate",scale="response")









