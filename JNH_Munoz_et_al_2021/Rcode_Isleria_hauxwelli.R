###############################JENNY MUNOZ#####################################
###############################################################################
##  Manuscript Nest description Isleria hauxwelli
## Journal format Wilson journal of ornithology
## R-code 
## Jenny Munoz
#### last update: February 2021
################################################################################
R.Version()
citation()
citation("incR")

# Assumptions": 

# Loading packages --------------------------------------------------------
install.packages("cooccur")
install.packages("car")
install.packages("vegan")
install.packages("AICcmodavg")
install.packages("stats")
install.packages("ggplot2")
install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
install.packages("lme4")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("visreg")
install.packages("skimr")
install.packages("ggthemes")

install.packages("incR")
library(incR)
packageVersion("incR")

##general packages  and visualization
library(lattice)
library(ggplot2); theme_set(theme_bw())
library(reshape2)
library(car)
library(visreg)
library (lme4)
library(tidyverse)
library(skimr)
library(ggthemes)

#Packages for biparental care
install.packages("tidyverse") #includer dplyr, tidyr, and ggplot2 
install.packages("dplyr")
install.packages("tidyr")

install.packages("ggplot2")
install.packages("ggthemes")
install.packages("extrafont")

# Load libraries
library(tidyverse)
library (dplyr)
library (tidyr)
library(ggplot2)
library(ggthemes)
library(extrafont)

#####Loading packages for Incubation pattern analyses
# Instaling package from v1.1.0.9000. https://github/PabloCapilla/incR
#install.packages("devtools")
#devtools::install_github(repo = "PabloCapilla/incR")
### Or for CRAN V 1.1.0 stable version
install.packages("incR")
library(incR)
packageVersion("incR")

# #Loading packages for Incubation pattern analyses -------------------
# install and load this packages to run this vignette
install.packages ("ggplot2")
install.packages ("dplyr")
install.packages ("data.table")
install.packages ("extrafont")
install.packages("tidyr")#I added this libraries for easier manipulation of data
install.packages("tidyverse") #I added this libraries for easier manipulation of data

library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(extrafont)
library(visreg)

loadfonts()


# Read the data ---------------------------------------------------------------
#growthrate<-read.csv(file="Isleria_hauxwelli_growth_rate.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
#growthrate<-read.csv(file="Isleria_hauxwelli_growth_rate1.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))

head(growthrate)
summary(growthrate)
str(growthrate)

#Exploring for missing data
which(is.na(growthrate)==T)

#Exploring the data
summary(growthrate) ## summary statistics for each column
unique(growthrate$date) ## unique values in a column

## the `skimr` package provides a really useful function for doing many of 
## these data exploration steps (and so much more!) in one simple command
skimr::skim(growthrate)


#Manuscript 

# Description egss, nestlings,nest  etc ----------------------------------------
####Eggs description
#Note that I excluded teh eggs number 3 that was found in one of the nest because I am not sure it represent the general eggmeasurements.
eggs<-read.csv ("Isleria_eggs_measurements.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
str(eggs)
##Lenght,width and mass

#Count the number of observations
eggs %>% tally(mass)
eggs %>% count (mass)

#Sumarize measurements
mean(eggs$egg_lenght)
sd(eggs$egg_lenght)
n(eggs$egg_lenght)


#####Nest description
nest<-read.csv(file="Isleria_nest_measurements.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
str(nest)

mean(nest$inner_length_nest)
sd(nest$inner_length_nest)

mean(nest$inner_width_nest)
sd(nest$inner_width_nest)

mean(nest$tickness)
sd(nest$tickness)

mean(nest$cup_depth)
sd(nest$cup_depth)

mean(nest$outher_length_nest)
sd(nest$outher_length_nest)
mean(nest$outher_width_nest)
sd(nest$outher_width_nest)

mean(nest$length)
sd(nest$length )

#nest height from ground
mean(nest$height)
sd(nest$height)

mean(nest$mass,na.rm=TRUE)
sd(nest$mass,na.rm=TRUE)

mean(nest$mass_layer_1,na.rm=TRUE)
sd(nest$mass_layer_1,na.rm=TRUE)

mean(nest$mass_layer_2,na.rm=TRUE)
sd(nest$mass_layer_2,na.rm=TRUE)

#####Nestling period
dayone<-filter(growthrate,True_date=="1")
mean(dayone$mass)
sd(dayone$mass)

fledgingday<-filter(growthrate,True_date=="9")
mean(fledgingday$mass)
sd(fledgingday$mass)

fledgingday<-filter(growthrate,True_date=="9")
mean(fledgingday$wing)
sd(fledgingday$wing)

mean(fledgingday$tail)
sd(fledgingday$tail)

mean(fledgingday$tarsus)
sd(fledgingday$tarsus)

#Weight change over time
str(growthrate)
mean(growthrate$mass_gained,na.rm=TRUE)
sd(growthrate$mass_gained,na.rm=TRUE)

#Adults size and measurements

adults<-read.csv(file.choose("Isleria_adult_measurements.csv"), stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))

str(adults)

mean(adults$Mass, na.rm=TRUE)
sd(adults$Mass, na.rm=TRUE)
mean(adults$Tail, na.rm=TRUE)
sd(adults$Tail, na.rm=TRUE)
mean(adults$Wing, na.rm=TRUE)
sd(adults$Wing, na.rm=TRUE)


# Incubation pattern analysis ---------------------------------------------
# Reading and prep the data (data need to be processed by incRprep and incRscan)-----------------------------------------------
#rawdata<-read.csv(file.choose("Myrmotherula_hauxwelli_P02a_AJRC13.csv"), stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
#rawdata<-read.csv(file="Isleria_selected.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
rawdata<-read.csv(file="Myrmotherula_hauxwelli_P02a_AJRC13.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
rawdata<-read.csv(file="Myrmotherulla_hauxwelli_P100a_MASM13.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))

str(rawdata)
max(rawdata$t_nest)
min(rawdata$t_nest)
mean(rawdata$t_nest)
sd(rawdata$t_nest)

max(rawdata$t_env)
min(rawdata$t_env)
mean(rawdata$t_env)
sd(rawdata$t_env)

hist(rawdata$t_nest)
hist(rawdata$t_env)

#all code together for practicity
### Format preparation, create an new data with the right format
##Note that some of the column content is missing, I ma not sure if this is right?
inc_sp<-incRprep(data=rawdata,
                 date.name="DATE", 
                 date.format="%d/%m/%Y %H:%M", 
                 timezone="GMT", 
                 temperature.name="t_nest")

head(inc_sp,3)

incubation.analysis <- incRscan (data=inc_sp,
                                 temp.name="t_nest",
                                 lower.time=18,
                                 upper.time=5,
                                 sensitivity=0.1,
                                 temp.diff.threshold=8,
                                 maxNightVariation=3,
                                 env.temp="t_env")
inc_scores_sp <- incubation.analysis[[1]]
inc.thresholds_sp <- incubation.analysis[[2]]

write.csv(inc_scores_sp, file = "inc_scores_isleria.csv")
inc_scores_sp<-read.csv(file="inc_scores_Isleria.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))

#subset
dayone<-inc_scores_sp%>%
  filter(date ==as.Date('2013-09-14') )

my_plot <- incRplot(data = dayone[complete.cases(dayone$t_nest),],
                    time.var = "dec_time",
                    day.var = "date",
                    inc.temperature.var = "t_nest",
                    env.temperature.var = "t_env",
                    vector.incubation = "incR_score")
my_plot

max(dayone$t_nest)
mean(dayone$t_env)
mean(dayone$t_nest)
#incubation
incRatt (data=inc_scores_sp,
         vector.incubation="incR_score")
incRact (data=inc_scores_sp,
         time_column="time",
         vector.incubation="incR_score")

incRbouts(data=inc_scores_sp,
          vector.incubation="incR_score",
          dec_time="dec_time",
          temp="temperature", #or temp1
          sampling.rate=240) # sampling rate in seconds.


#day only
day_inc_scores_sp<-inc_scores_sp %>% 
  filter(hour>=5 & hour<=18 )
incRatt (data=day_inc_scores_sp,
         vector.incubation="incR_score")

my_plot <- incRplot(data = day_inc_scores_sp[complete.cases(day_inc_scores_sp$t_nest),],
                    time.var = "dec_time",
                    day.var = "date",
                    inc.temperature.var = "t_nest",
                    env.temperature.var = "t_env",
                    vector.incubation = "incR_score")
my_plot

write.csv(day_inc_scores_sp, file = "inc_scores_Isleria.csv")


# Nest attentiveness ------------------------------------------------------
#Extracting nest attentiveness from sensors

inc_scores_Isleria_sensor<-read.csv("inc_scores_Isleria.csv", stringsAsFactors = FALSE, header=TRUE, strip.white=TRUE, na.strings = c("NA",""))
str(inc_scores_Isleria_sensor)

daily_nest_attentiveness<-inc_scores_Isleria_sensor%>%
  filter(between (hour_minute,"5.31","17.30"))%>%
  group_by(date)%>%
  summarize(sumincR_score=sum(incR_score))%>% 
  mutate(daily_inc_percentage=((sumincR_score-6)*100)/ 720)
  
#Note 720 is the total minutes in 12 hours (530-1730_), will be 100 percent of incubating all the day
#but because there is an error in the filter of the data I have 6 minutes more each day, 5.4, 5.5, 5.6 , 5.7, 5.8, 5.9 so I -6min per day todatl 14 days 

write.csv(daily_nest_attentiveness, file = "daily_nest_attentiveness_sensors_.csv")

#Nest attentiveness from camera traps

inc_scores_isleria_camera<-read.csv("nest_attentiveness_isleria_eggs.csv", stringsAsFactors = FALSE, header=TRUE, strip.white=TRUE, na.strings = c("NA",""))

daily_nest_attentiveness_camera<-inc_scores_isleria_camera%>%
  filter(nest_id=="P02_AJRC13")%>%
  filter(event=="on")%>%
  filter(nocturnal_incubation=="N")%>%
  filter(date!="29-Sep-13")%>% 
  filter(date!="23-Sep-13")%>% 
  filter(date!="08-Sep-13")%>% 
  group_by(date)%>%
  summarise(inc_total_daily_min=sum(duration_min))

nest_attentiveness_camera_night<-inc_scores_isleria_camera%>%
  filter(nest_id=="P02_AJRC13")%>%
  filter(event=="on")%>%
  filter(nocturnal_incubation=="Y")%>%
  filter(date!="29-Sep-13")%>% 
  filter(date!="23-Sep-13")%>% 
  filter(date!="08-Sep-13")%>% 
  group_by(date)%>%
  summarise(inc_total_night_min=sum(duration_min))


nest_attentiveness_camera_calculation<-full_join(daily_nest_attentiveness_camera,nest_attentiveness_camera_night, by="date")

nest_attentiveness_camera_calculation<-nest_attentiveness_camera_calculation%>%
mutate(total_daytime_lenght=(1440-inc_total_night_min))

nest_attentiveness_camera_calculation<-nest_attentiveness_camera_calculation%>%
  mutate(daily_nest_attentiveness=((inc_total_daily_min*100)/total_daytime_lenght))


#write.csv(nest_attentiveness_camera_calculation, file = "daily_nest_attentiveness_camera_calculation.csv")
                              

# Biparental care from camera traps ---------------------------------------
# Graphs and analysis on biparental care#####
########read the data
biparental<-read.csv("nest_attentiveness_isleria_eggs.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
biparental1<-read.csv("nest_attentiveness_nestlings.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
biparentalhours<-read.csv(file.choose("nest_attentiveness_isleria_hours"),stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
str(biparentalhours)

#Tidying the data
#select specific columns and filter
select(biparental, sex) %>%    
  filter (sex=="f")

#Lets calculate new variables that summarise info
#ANALYSING BY DAY # First we need to exclude dates when the info is incomplete sept 8, sept 23, and sept 29 2014 from nest P02_AJRC13 are incomplete, maybe camera was off or not batteries
#calculate total duration on "incubation time" per day for male and female.
#a######INCUBATION TIME ANALYSES
#Trying to get info for the other nest

#incubation nest P02_IGZ14
inc_time_female_total<-biparental%>%
  filter(nest_id=="P02_IGZ14")%>%
  filter(sex=="f")%>%
  filter(event=="on")%>% 
  group_by(date)%>%
  summarise(inc_female_total_h=sum(duration_h))

inc_time_female_daily<-biparental%>%
  filter(sex=="f")%>%
  filter(event=="on")%>%
  filter(nocturnal_incubation=="N")%>%
  filter(nest_id=="P02_IGZ14")%>%
  group_by(date)%>%
  summarise(inc_female_daily_h=sum(duration_h))

inc_time_male_total<-biparental%>%
  filter(nest_id=="P02_IGZ14")%>%
  filter(sex=="m")%>%
  filter(event=="on")%>% 
  group_by(date)%>%
  summarise(inc_male_total_h=sum(duration_h))

inc_time_male_daily<-biparental%>%
  filter(sex=="m")%>%
  filter(event=="on")%>%
  filter(nocturnal_incubation=="N")%>%
  filter(nest_id=="P02_IGZ14")%>%
  group_by(date)%>%
  summarise(inc_male_daily_h=sum(duration_h))

inc_time_total_daily<-biparental%>%
  filter(event=="on")%>%
  filter(nocturnal_incubation=="N")%>%
  filter(nest_id=="P02_IGZ14")%>%
  group_by(date)%>%
  summarise(inc_total_daily_h=sum(duration_h))
#Combine DATA SETS

inc_male_female<-full_join(inc_time_female_daily,inc_time_male_daily,by="date")
incubation_summary<-full_join(inc_male_female,inc_time_total_daily,by="date")
incubation_summary<-full_join(incubation_summary,inc_time_female_total,by="date")

true_inc_day<-biparental%>%
  filter(nest_id=="P02_IGZ14")%>%
  select(true_inc_day,date)

true_inc_day<-distinct(true_inc_day) #  removeduplicate rows
incubation_summary<-full_join(true_inc_day, incubation_summary, by="date")
str(incubation_summary)

write.csv(incubation_summary, "incubation_summary_P02_IGZ14.csv")

#Initial nest

#female incubation
inc_time_female_total<-biparental%>%
  filter(sex=="f")%>%
  filter(event=="on")%>%
  filter(date!="23-Sep-13")%>% 
  filter(date!="29-Sep-13")%>% 
  filter(date!="08-Sep-13")%>% 
  filter(nest_id=="P02_AJRC13")%>%
  group_by(date)%>%
  summarise(inc_female_total_min=sum(duration_min))

#female incubation nest P02_AJRC13

inc_time_female_daily<-biparental%>%
  filter(sex=="f")%>%
  filter(event=="on")%>%
  filter(nocturnal_incubation=="N")%>%
  filter(date!="23-Sep-13")%>% 
  filter(date!="29-Sep-13")%>% 
  filter(date!="08-Sep-13")%>% 
  filter(nest_id=="P02_AJRC13")%>%
  group_by(date)%>%
  summarise(inc_female_daily_min=sum(duration_min))
  
#male incubation time
inc_time_male_daily<-biparental%>%
  filter(sex=="m")%>%
  filter(event=="on")%>%
  filter(nocturnal_incubation=="N")%>%
  filter(date!="23-Sep-13")%>% 
  filter(date!="29-Sep-13")%>% 
  filter(date!="08-Sep-13")%>% 
  filter(nest_id=="P02_AJRC13")%>%
  group_by(date)%>%
  summarise(inc_male_daily_min=sum(duration_min))

#total incubation daily
inc_time_total_daily<-biparental%>%
  filter(event=="on")%>%
  filter(nocturnal_incubation=="N")%>%
  filter(date!="23-Sep-13")%>% 
  filter(date!="29-Sep-13")%>% 
  filter(date!="08-Sep-13")%>% 
  filter(nest_id=="P02_AJRC13")%>%
  group_by(date)%>%
  summarise(inc_total_daily_min=sum(duration_min))

#combine data sets
inc_male_female<-full_join(inc_time_female_daily,inc_time_male_daily,by="date")
incubation_summary<-full_join(inc_male_female,inc_time_total_daily,by="date")
incubation_summary<-full_join(incubation_summary,inc_time_female_total,by="date")

# change minutes to hours, and create new variable that sumarize femaletotalnight+male

incubation_summary<-incubation_summary%>%
  filter(date!="13-Sep-13")%>% #this day was excluded because seems to be an errr on the analyses, only 3 min of incubation this day
  mutate(inc_female_daily_h=inc_female_daily_min/60)%>%
  mutate(inc_male_daily_h=inc_male_daily_min/60)%>%
  mutate(inc_total_daily_h=inc_total_daily_min/60)%>%
  mutate(inc_female_total_h=inc_female_total_min/60)%>%
  mutate(inc_total_h=inc_female_total_h+inc_male_daily_h)
#create a new colum with true inc day and combine datasets
true_inc_day<-biparental%>%
  filter(nest_id=="P02_AJRC13")%>%
  filter(date!="23-Sep-13")%>% 
  filter(date!="29-Sep-13")%>% 
  filter(date!="08-Sep-13")%>% 
  filter(date!="13-Sep-13")%>% 
  select(true_inc_day,date)

true_inc_day<-distinct(true_inc_day) #  removeduplicate rows
incubation_summary<-full_join(true_inc_day, incubation_summary, by="date")
str(incubation_summary)

incubation_vertical<-gather(incubation_summary, "inc_female_daily_H","inc_male_daily_h","inc_total_daily_h",3:5 )




#export the summary as a csv
write.csv(incubation_summary, "incubation_summary.csv")

# organize all the info in one column 
incubation_vertical<-gather(incubation_summary, "inc_female_daily_min","inc_male_daily_min","inc_total_daily_min",3:5 )

incubation_vertical<-incubation_vertical%>%
  rename (incubation=inc_female_daily_min)%>%
  rename (incubation_time=inc_male_daily_min)

#let's plot
incubation_summary<-read.csv("incubation_summary.csv")
str(incubation_summary)
# Incubation daily
pdf("Fig6a.Dailyincubationbiparental.pdf",width=10,height=10, useDingbats=FALSE)
plot(inc_female_daily_h~true_inc_day, data=incubation_summary, type="o", col="black", ylim=c(0,10), pch=19, cex.lab=2, cex=1.5, xlab='Day', ylab='')
axis(1, seq(0, 23, by = 1), las=0) 
title(ylab="Daylight Incubation (h/day)", line=2, cex.lab=2)
lines(inc_male_daily_h~true_inc_day, data=incubation_summary, type="o", pch=19, lty=1, col="blue", cex=1.5)
lines(inc_total_daily_h~true_inc_day, data=incubation_summary, type="l", pch=22, lty=3, col="black")
dev.off()

#Incubation including night data

pdf("Fig6b.totalincubationbiparental.pdf",width=10,height=10, useDingbats=FALSE)
plot(inc_female_total_h~true_inc_day, data=incubation_summary, type="o", col="black", ylim=c(0,25), pch=19, cex.lab=2, cex=1.5, xlab='Day', ylab='')
axis(1, seq(0, 23, by = 1), las=0) 
title(ylab="Incubation(h/day)", line=2, cex.lab=2)
lines(inc_male_daily_h~true_inc_day, data=incubation_summary, type="o", pch=19, lty=1, col="blue", cex=1.5)
lines(inc_total_h~true_inc_day, data=incubation_summary, type="l", pch=22, lty=3, col="black")
dev.off()

#b FORAGING TRIPS
#female
foraging_trips_female<-biparental%>% 
  filter(date!="23-Sep-13")%>% 
  filter(date!="29-Sep-13")%>% 
  filter(date!="08-Sep-13")%>% 
  filter(sex=="f")%>%
  filter(event=="off")%>%
  filter(nest_id=="P02_AJRC13")%>%
  group_by(date)%>%
  count(event)%>%
  rename (trips_female=n)
  
#male
foraging_trips_male<-biparental%>% 
  filter(date!="23-Sep-13")%>% 
  filter(date!="29-Sep-13")%>% 
  filter(date!="08-Sep-13")%>% 
  filter(sex=="m")%>%
  filter(event=="off")%>%
  filter(nest_id=="P02_AJRC13")%>%
  group_by(date)%>%
  count(event)%>%
  rename (trips_male=n)

#total foraging trips
foraging_trips_total<-biparental%>% 
  filter(nocturnal_incubation=="N")%>%
  filter(date!="23-Sep-13")%>% 
  filter(date!="29-Sep-13")%>% 
  filter(date!="08-Sep-13")%>% 
  filter(event=="off")%>%
  filter(nest_id=="P02_AJRC13")%>%
  group_by(date)%>%
  count(event)%>%
  rename (trips_total=n)

#combine data sets
trips_female_male<-full_join(foraging_trips_female,foraging_trips_male,by="date")
trips_summary<-full_join(trips_female_male,foraging_trips_total,by="date")

true_inc_day<-biparental%>%
  filter(nest_id=="P02_AJRC13")%>%
  filter(date!="23-Sep-13")%>% 
  filter(date!="29-Sep-13")%>% 
  filter(date!="08-Sep-13")%>% 
  filter(date!="13-Sep-13")%>% 
  select(true_inc_day,date)

true_inc_day<-distinct(true_inc_day) #  removeduplicate rows
trips_summary<-full_join(true_inc_day, trips_summary, by="date")

#export summary as csv
#write.csv(trips_summary, "foraging_trips_summary.csv")

#plot number of trips
pdf("Fig6c.totaltrips.pdf",width=10,height=10, useDingbats=FALSE)
plot(trips_female~true_inc_day, data=trips_summary, type="o", col="black", ylim=c(0,17), pch=19, cex.lab=2, cex=1.5, xlab='Day', ylab='')
axis(1, seq(0, 24, by = 1), las=0) 
title(ylab="Foraging trips (events/day)", line=2, cex.lab=2)
lines(trips_male~true_inc_day, data=trips_summary, type="o", pch=19, lty=1, col="blue", cex=1.5)
lines(trips_total~true_inc_day, data=trips_summary, type="l", pch=22, lty=3, col="black")
dev.off()

# VISIT TRIPS
#female
visits_female<-biparental%>% 
  filter(date!="23-Sep-13")%>% 
  filter(date!="29-Sep-13")%>% 
  filter(date!="08-Sep-13")%>% 
  filter(sex=="f")%>%
  filter(event=="visita")%>%
  filter(nest_id=="P02_AJRC13")%>%
  group_by(date)%>%
  count(event)%>%
  rename (female_visits=n)

#male
visits_male<-biparental%>% 
  filter(date!="23-Sep-13")%>% 
  filter(date!="29-Sep-13")%>% 
  filter(date!="08-Sep-13")%>% 
  filter(sex=="m")%>%
  filter(event=="visita")%>%
  filter(nest_id=="P02_AJRC13")%>%
  group_by(date)%>%
  count(event)%>%
  rename(male_visits=n)
#total
visits_total<-biparental%>% 
  filter(date!="23-Sep-13")%>% 
  filter(date!="29-Sep-13")%>% 
  filter(date!="08-Sep-13")%>% 
  filter(event=="visita")%>%
  filter(nest_id=="P02_AJRC13")%>%
  group_by(date)%>%
  count(event)%>%
  rename (total_visits=n)

#combine data sets
visits_female_male<-full_join(visits_female,visits_male,by="date")
visits_summary<-full_join(visits_female_male,visits_total,by="date")

true_inc_day<-biparental%>%
  filter(nest_id=="P02_AJRC13")%>%
  filter(date!="23-Sep-13")%>% 
  filter(date!="29-Sep-13")%>% 
  filter(date!="08-Sep-13")%>% 
  filter(date!="13-Sep-13")%>% 
  select(true_inc_day,date)

true_inc_day<-distinct(true_inc_day) #  removeduplicate rows
visits_summary<-full_join(true_inc_day, visits_summary, by="date")

#write.csv(visits_summary, "nest_visits_summary.csv")

#plot number of visits
pdf("Fig6d.totalvisits.pdf",width=10,height=10, useDingbats=FALSE)
plot(female_visits~true_inc_day, data=visits_summary, type="o", col="black", ylim=c(0,17), pch=19, cex.lab=2, cex=1.5, xlab='Day', ylab='')
axis(1, seq(0, 24, by = 1), las=0) 
title(ylab="Nest visits (events/day)", line=2, cex.lab=2)
lines(male_visits~true_inc_day, data=visits_summary, type="o", pch=19, lty=1, col="blue", cex=1.5)
lines(total_visits~true_inc_day, data=visits_summary, type="l", pch=22, lty=3, col="black")
dev.off()

######FEEDING
#female
feeding_female<-biparental1%>% 
  filter(sex=="f")%>%
  filter(event=="Feeding")%>%
  filter(nest_id=="P02_IGZ14")%>%
  group_by(date)%>%
  count(event)%>%
  rename (feeding_female=n)

#male
feeding_male<-biparental1%>% 
  filter(sex=="m")%>%
  filter(event=="Feeding")%>%
  filter(nest_id=="P02_IGZ14")%>%
  group_by(date)%>%
  count(event)%>%
  rename (feeding_male=n)

#total foraging trips
feeding_total<-biparental1%>% 
  filter(event=="Feeding")%>%
  filter(nest_id=="P02_IGZ14")%>%
  group_by(date)%>%
  count(event)%>%
  rename (feeding_total=n)

#combine data sets
feeding_female_male<-full_join(feeding_female,feeding_male,by="date")
feeding_summary<-full_join(feeding_female_male,feeding_total,by="date")

true_day<-biparental1%>%
  filter(nest_id=="P02_IGZ14")%>%
  select(true_day,date)

true_day<-distinct(true_day) #  removeduplicate rows
feeding_summary<-full_join(true_day, feeding_summary, by="date")

#write.csv(feeding_summary, "feeding_events_summary_isleria.csv")

feeding_summary<-read_csv("feeding_events_summary_isleria.csv")

feeding_summary
#plot feeding frequency per day
pdf("Fig6e.Feeding frequency per day.pdf",width=10,height=10, useDingbats=FALSE)
plot(feeding_female~true_day, data=feeding_summary, type="o", col="black", ylim=c(0,30), pch=19, cex.lab=2, cex=1.5, xlab='Day', ylab='')
axis(1, seq(0, 15, by = 1), las=0) 
title(ylab="Nestling provisioning (events/day)",line=2, cex.lab=2)
lines(feeding_male~true_day, data=feeding_summary, type="o", pch=19, lty=1, col="blue", cex=1.5)
lines(feeding_total~true_day, data=feeding_summary, type="l", pch=22, lty=3, col="black")
dev.off()

####Brooding
#female 
brooding_time_female_total<-biparental1%>%
  filter(sex=="f")%>%
  filter(event=="brooding")%>%
  filter(nest_id=="P02_IGZ14")%>%
  group_by(date)%>%
  summarise(brooding_female_total_h=sum(duration_h))

brooding_time_female_daily<-biparental1%>%
  filter(sex=="f")%>%
  filter(event=="brooding")%>%
  filter(nocturnal_incubation=="N")%>%
  filter(nest_id=="P02_IGZ14")%>%
  group_by(date)%>%
  summarise(brooding_female_daily_h=sum(duration_h))

#male incubation time


brooding_time_male_daily<-biparental1%>%
  filter(sex=="m")%>%
  filter(event=="brooding")%>%
  filter(nocturnal_incubation=="N")%>%
  filter(nest_id=="P02_IGZ14")%>%
  group_by(date)%>%
  summarise(brooding_male_daily_h=sum(duration_h))

#total incubation daily
brooding_time_total_daily<-biparental1%>%
  filter(event=="brooding")%>%
  filter(nocturnal_incubation=="N")%>%
  group_by(date)%>%
  summarise(brooding_total_daily_h=sum(duration_h))

#combine data sets
brooding_male_female<-full_join(brooding_time_female_daily,brooding_time_male_daily,by="date")
brooding_summary<-full_join(brooding_male_female,brooding_time_total_daily,by="date")
brooding_summary<-full_join(brooding_summary,brooding_time_female_total,by="date")

#true_day_column
true_day<-biparental1%>%
  filter(nest_id=="P02_IGZ14")%>%
  select(true_day,date)

true_day<-distinct(true_day) #  removeduplicate rows
brooding_summary<-full_join(true_day, brooding_summary, by="date")

#create the variable 

brooding_summary<-brooding_summary%>%
  mutate (brooding_total_h=brooding_female_total_h+brooding_male_daily_h)

#write.csv(brooding_summary, "brooding_events_summary_isleria.csv")

#plots

pdf("Fig6f.Broodingdaylight.pdf",width=10,height=10, useDingbats=FALSE)
plot(brooding_female_daily_h~true_day, data=brooding_summary, type="o", col="black", ylim=c(0,10), pch=19, cex.lab=2, cex=1.5, xlab='Day', ylab='')
axis(1, seq(0, 15, by = 1), las=0) 
title(ylab="Daylight brooding (h/day)",line=2, cex.lab=2)
lines(brooding_male_daily_h~true_day, data=brooding_summary, type="o", pch=19, lty=1, col="blue", cex=1.5)
lines(brooding_total_daily_h~true_day, data=brooding_summary, type="l", pch=22, lty=3, col="black")
dev.off()

pdf("Fig6g.Broodingtotal.pdf",width=10,height=10, useDingbats=FALSE)
plot(brooding_female_total_h~true_day, data=brooding_summary, type="o", col="black", ylim=c(0,25), pch=19, cex.lab=2, cex=1.5, xlab='Day', ylab='')
axis(1, seq(0, 15, by = 1), las=0) 
title(ylab="Brooding (h/day)",line=2, cex.lab=2)
lines(brooding_male_daily_h~true_day, data=brooding_summary, type="o", pch=19, lty=1, col="blue", cex=1.5)
lines(brooding_total_h~true_day, data=brooding_summary, type="l", pch=22, lty=3, col="black")
dev.off()

--
# Biparental care behavior (Summary stats) ---------------------------------------------------------

#Incubation
biparental_incubation<-read.csv(file="incubation_summary_isleria.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
str(biparental_incubation)
View(biparental_incubation)

#get some descriptive stats
incubation_time<-filter(biparental_incubation,sex=="total_f&m")
incubation_time 
mean_inc_h<-mean(incubation_time$inc_daily_h)
mean_inc_h
percentage_incubation<-(mean_inc_h*100)/11
percentage_incubation
max(incubation_time$inc_daily_h)
min(incubation_time$inc_daily_h)
sd(incubation_time$inc_daily_h)

mean_inc_h_overall<-mean(incubation_time$inc_daynight_h)
mean_inc_h_overall
(mean_inc_h_overall*100)/24
max(incubation_time$inc_daynight_h)
min(incubation_time$inc_daynight_h)
sd(incubation_time$inc_daynight_h)

#Female
incubation_time<-filter(biparental_incubation,sex=="f")
incubation_time 
mean_inc_h<-mean(incubation_time$inc_daily_h)
mean_inc_h
max(incubation_time$inc_daily_h)
min(incubation_time$inc_daily_h)
sd(incubation_time$inc_daily_h)

mean_inc_h_overall<-mean(incubation_time$inc_daynight_h)
mean_inc_h_overall
max(incubation_time$inc_daynight_h)
min(incubation_time$inc_daynight_h)
sd(incubation_time$inc_daynight_h)

#male
incubation_time<-filter(biparental_incubation,sex=="m")
incubation_time 
mean_inc_h<-mean(incubation_time$inc_daily_h)
mean_inc_h
max(incubation_time$inc_daily_h)
min(incubation_time$inc_daily_h)
sd(incubation_time$inc_daily_h)

mean_inc_h_overall<-mean(incubation_time$inc_daynight_h)
mean_inc_h_overall
max(incubation_time$inc_daynight_h)
min(incubation_time$inc_daynight_h)
sd(incubation_time$inc_daynight_h)


#standard_error <- function(x) sd(x) / sqrt(length(x))
#standard_error(x)
#x<-(incubation_time$inc_daynight_h)

#Incubation graphs

pdf("Fig.Daytime incubation.pdf",width=14,height=8, useDingbats=FALSE) #export as pdf
boxplot(inc_daily_h~sex*stage_eggs, data=biparental_incubation,main="Daytime Incubation", ylab="Daytime incubation (h/day)",xlab="Incubation stage", ylim=c(0,10), boxwex=0.5, cex.axis=1.5,cex=0.5,col=(c("white","darkslategray4", "lightgrey")))
stripchart(inc_daily_h~code_day,data=biparental_incubation,  add=TRUE,vertical=TRUE, method='jitter', jitter=0.12, pch=19, cex=2, ylim=c(0,10))
dev.off()
#method='jitter', jitter=0.03,

pdf("Fig.Incubation.pdf",width=14,height=8, useDingbats=FALSE) #export as pdf
boxplot(inc_daynight_h~sex_night*stage_eggs, data=biparental_incubation,main=" Incubation", ylab="Incubation time (h/day)",xlab="Incubation stage", ylim=c(0,20), boxwex=0.5, col=(c("white","darkslategray4")))
stripchart(inc_daynight_h~code_night,data=biparental_incubation, method='jitter', jitter=0.12, add=TRUE,vertical=TRUE, pch=19, cex=2, ylim=c(0,20))
dev.off()

####
###Brooding
####
biparental_brooding<-read.csv(file="brooding_events_summary_isleria.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
biparental_brooding<-filter(biparental_brooding,date!="02-Sep-14") #sept 2 was filterd cause this was probably the hatching day,and the ndata is difficult to

str(biparental_brooding)
View(biparental_brooding)

#brooding graphs
pdf("Fig.Daytime brooding.pdf",width=14,height=8, useDingbats=FALSE) #export as pdf
boxplot(brooding_daily_h~sex*nestlings_stage, data=biparental_brooding,main="Daytime brooding", ylab="Daylight brooding time (h/day)",xlab="Incubation stage", ylim=c(0,10), boxwex=0.5, col=(c("white","darkslategray4", "lightgrey")))
stripchart(brooding_daily_h~code_day,data=biparental_brooding, method='jitter', jitter=0.12, add=TRUE,vertical=TRUE, pch=19, cex=2, ylim=c(0,10))
dev.off()

pdf("Fig.Brooding.pdf",width=14,height=8, useDingbats=FALSE) #export as pdf
boxplot(brooding_daynight_h~sex_night*nestlings_stage, data=biparental_brooding, main="Brooding",ylab="Total brooding time(h/day)",xlab="Nestlings stage", ylim=c(0,20), boxwex=0.5, col=(c("white","darkslategray4")))
stripchart(brooding_daynight_h~code_night,data=biparental_brooding, method='jitter', jitter=0.12, add=TRUE,vertical=TRUE, pch=19, cex=2, ylim=c(0,20))
dev.off()

# descriptive stats brooding
biparental_brooding<-filter(biparental_brooding,date!="02-Sep-14") #sept 2 was filterd cause this was probably the hatching day,and the ndata is difficult to
brooding_time<-filter(biparental_brooding,sex=="total_f&m")
brooding_time

mean_brood_h<-mean(brooding_time$brooding_daily_h)
mean_brood_h
max(brooding_time$brooding_daily_h)
min(brooding_time$brooding_daily_h)
sd(brooding_time$brooding_daily_h)

mean_brood_overall<-mean(brooding_time$brooding_daynight_h)
mean_brood_overall
max(brooding_time$brooding_daynight_h)
min(brooding_time$brooding_daynight_h)
sd(brooding_time$brooding_daynight_h)

#female
biparental_brooding<-filter(biparental_brooding,date!="02-Sep-14") #sept 2 was filterd cause this was probably the hatching day,and the ndata is difficult to
brooding_time<-filter(biparental_brooding,sex=="f")
brooding_time

mean_brood_h<-mean(brooding_time$brooding_daily_h)
mean_brood_h
max(brooding_time$brooding_daily_h)
min(brooding_time$brooding_daily_h)
sd(brooding_time$brooding_daily_h)

mean_brood_overall<-mean(brooding_time$brooding_daynight_h)
mean_brood_overall
max(brooding_time$brooding_daynight_h)
min(brooding_time$brooding_daynight_h)
sd(brooding_time$brooding_daynight_h)

#male
biparental_brooding<-filter(biparental_brooding,date!="02-Sep-14") #sept 2 was filterd cause this was probably the hatching day,and the ndata is difficult to
brooding_time<-filter(biparental_brooding,sex=="m")
brooding_time

mean_brood_h<-mean(brooding_time$brooding_daily_h)
mean_brood_h
max(brooding_time$brooding_daily_h)
min(brooding_time$brooding_daily_h)
sd(brooding_time$brooding_daily_h)

mean_brood_overall<-mean(brooding_time$brooding_daynight_h)
mean_brood_overall
max(brooding_time$brooding_daynight_h)
min(brooding_time$brooding_daynight_h)
sd(brooding_time$brooding_daynight_h)


#percentage_brooding<-(mean_brood_h*100)/10
#percentage_brooding(mean_inc_h_overall*100)/24

#provisioning rate
### important to note in the grph the total number of events in a given day 
#can be higher than m+female cause sometimes it was not posssible to id the sex of the bird. but it was documented as a feeding event. This is the case for example in the last boxplot.
biparental_feeding<-read.csv(file="feeding_events_summary_isleria.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
str(biparental_feeding)
View (biparental_feeding)

pdf("Fig.Feeding.pdf",width=14,height=8, useDingbats=FALSE) #export as pdf
boxplot(feeding_events~sex*nestling_stage, data=biparental_feeding,main="Nestlings provisioning", ylab="Feeding frequency (n/day)",xlab="Incubation stage", ylim=c(0,30), boxwex=0.5, col=(c("white","darkslategray4", "lightgrey")))
stripchart(feeding_events~code,data=biparental_feeding, method='jitter', jitter=0.12, add=TRUE,vertical=TRUE, pch=19, cex=2, ylim=c(0,30))
dev.off()

plot(true_day~feeding_events_1,data=biparental_feeding, method='jitter', jitter=0.12, pch=19, cex=0.8, ylim=c(0,20), col=(c("white","darkslategray4", "lightgrey"))) 

#get some descriptive stats
feeding_frequency<-filter(biparental_feeding,sex=="total_f&m")
feeding_frequency
mean_feeding<-mean(feeding_frequency$feeding_events)
mean_feeding

max(feeding_frequency$feeding_events)
min(feeding_frequency$feeding_events)
sd(feeding_frequency$feeding_events)

#female
feeding_frequency<-filter(biparental_feeding,sex=="f")
feeding_frequency
mean_feeding<-mean(feeding_frequency$feeding_events)
mean_feeding

max(feeding_frequency$feeding_events)
min(feeding_frequency$feeding_events)
sd(feeding_frequency$feeding_events)

#male
feeding_frequency<-filter(biparental_feeding,sex=="m")
feeding_frequency
mean_feeding<-mean(feeding_frequency$feeding_events)
mean_feeding

max(feeding_frequency$feeding_events)
min(feeding_frequency$feeding_events)
sd(feeding_frequency$feeding_events)


#foraging trips
#number of foraging trips.
biparental_trips<-read.csv(file="foraging_trips_summary_isleria.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
str(biparental_trips)
#during incubation
pdf("Fig.Trips.pdf",width=14,height=8, useDingbats=FALSE) #export as pdf
boxplot(trips~sex*stage_eggs, data=biparental_trips,main="Foraging trips", ylab=" Trips frequency (trips/day)",xlab="", ylim=c(0,15), boxwex=0.5, col=(c("white","darkslategray4", "lightgrey")))
stripchart(trips~code,data=biparental_trips, method='jitter', jitter=0.05, add=TRUE,vertical=TRUE, pch=19, cex=2, ylim=c(0,30))
dev.off()

#get some descriptive stats
trips_frequency<-filter(biparental_trips,sex=="total_f&m")
trips_frequency
mean_trips<-mean(trips_frequency$trips)
mean_trips
sd(trips_frequency$trips)

max(trips_frequency$trips)
min(trips_frequency$trips)

#female
trips_frequency<-filter(biparental_trips,sex=="f")
trips_frequency
mean_trips<-mean(trips_frequency$trips)
mean_trips

max(trips_frequency$trips)
min(trips_frequency$trips)
sd(trips_frequency$trips)

#male
trips_frequency<-filter(biparental_trips,sex=="m")
trips_frequency
mean_trips<-mean(trips_frequency$trips)
mean_trips
sd(trips_frequency$trips)

max(trips_frequency$trips)
min(trips_frequency$trips)

##duration of foraging trips
biparental_trips_duration<-read.csv(file="foraging_trips_duration_isleria.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
str(biparental_trips_duration)


pdf("Fig.Foraging_trips_duration.pdf",width=14,height=8, useDingbats=FALSE) #export as pdf
boxplot(duration_min~stage_eggs, data=biparental_trips_duration,main="Foraging trips duration", ylab="Duration(min)",xlab="Incubation stage", ylim=c(0,300), boxwex=0.5)
stripchart(duration_min~stage_eggs, data=biparental_trips_duration,  add=TRUE,vertical=TRUE, method='jitter', jitter=0.03, pch=19, cex=1.5, ylim=c(0,300))
dev.off()

lm <- lm(duration_min~true_inc_day, data = biparental_trips_duration)
plot(duration_min~true_inc_day, data=biparental_trips_duration, ylim=c(0,300), xlim=c(0,30))
par(new = TRUE)
visreg(lm, ylab=NA, xlab=NA, color="black", overlay=TRUE,width=14,height=8, useDingbats=FALSE, ylim=c(0,300),  xlim=c(0,30))

##summary statistics

mean(biparental_trips_duration$duration_min)
sd(biparental_trips_duration$duration_min)
max(biparental_trips_duration$duration_min)
min(biparental_trips_duration$duration_min)


# Plots for final manuscript ----------------------------------------------

biparental_incubation<-read.csv("incubation_summary_isleria.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
#nest_attentiveness<-read.csv("nest_attentiveness_by_method.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
biparental_brooding<-read.csv("brooding_events_summary_isleria.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
biparental_trips_duration<-read.csv("foraging_trips_duration_isleria.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
biparental_trips<-read.csv("foraging_trips_summary_isleria.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
biparental_feeding<-read.csv("feeding_events_summary_isleria.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))

#Total incubation 
biparental_incubation<-biparental_incubation%>%
  filter(sex!="total_f&m") 

str(biparental_incubation)
View(biparental_incubation)

qqnorm(biparental_incubation$inc_daynight_min)
qqline(biparental_incubation$inc_daynight_min, lty=2)
shapiro.test(biparental_incubation$inc_daynight_min) # iF P<0.05 data deviate from normality

mean (biparental_incubation$inc_daynight_min)
var(biparental_incubation$inc_daynight_min)

ggplot(biparental_incubation, aes(y=inc_daynight_min , x=true_inc_day, colour=factor(sex),shape=factor(sex)))+
  geom_point(aes(y=inc_daynight_min , x=true_inc_day), size=4)+
  scale_shape_manual(values=c(1,19))+
  geom_smooth(method = "glm",method.args=list(family ="quasipoisson" (link = "log")), fill="grey")+
  scale_color_manual(breaks = c("male", "female"),values=c("black", "darkgrey"))+
  geom_rug()+
  labs(title="", y = "Total incubation (min/day)", x = "Incubation day")+
  theme_classic(base_size = 22)+ 
  theme(legend.position = "none")

#Incubation day time
biparental_incubation<-biparental_incubation%>%
  filter(sex!="total_f&m") 

qqnorm(biparental_incubation$inc_daily_min)
qqline(biparental_incubation$inc_daily_min, lty=2)
shapiro.test(biparental_incubation$inc_daily_min) # iF P<0.05 data deviate from normality

mean (biparental_incubation$inc_daily_min)
var(biparental_incubation$inc_daily_min)

str(biparental_incubation)
View(biparental_incubation)

ggplot(biparental_incubation, aes(y=inc_daily_min, x=true_inc_day, colour=factor(sex),shape=factor(sex)))+
  geom_point(aes(y=inc_daily_min, x=true_inc_day), size=4)+
  scale_shape_manual(values=c(1,19))+
  geom_smooth(method = "glm",method.args=list(family ="quasipoisson" (link = "log")), fill="grey")+
  scale_color_manual(breaks = c("male", "female"),values=c("black", "darkgrey"))+
  geom_rug()+
  labs(title="",y = " Daytime incubation (min/day)", x = "Incubation day")+
  theme_classic(base_size = 22)+ 
  theme(legend.position = "none")

#Total brooding
str(biparental_brooding)
qqnorm(biparental_brooding$brooding_daynight_min)
qqline(biparental_brooding$brooding_daynight_min, lty=2)
shapiro.test(biparental_brooding$brooding_daynight_min) # iF P<0.05 data deviate from normality

mean(biparental_trips$trips)
var (biparental_trips$trips)

biparental_brooding<-biparental_brooding%>%
  filter(sex!="total_f&m")%>%
  mutate(brooding_daynight_min=brooding_daynight_h*60)

mean(biparental_brooding$brooding_daynight_min)
var (biparental_brooding$brooding_daynight_min) 

str(biparental_brooding)
View(biparental_brooding)

ggplot(biparental_brooding, aes(y=brooding_daynight_min, x=true_day, colour=factor(sex),shape=factor(sex)))+
  geom_point(aes(y=brooding_daynight_min, x=true_day), size=4)+
  scale_shape_manual(values=c(1,19))+
  geom_smooth(method = "glm",method.args=list(family ="quasipoisson" (link = "log")), fill="grey")+
  scale_color_manual(breaks = c("male", "female"),values=c("black", "darkgrey"))+
  geom_rug()+
  labs(title="",y = " Total brooding (min/day)", x = "Days after hatching")+
  scale_x_continuous (limits = c(0, 11), breaks = seq(0, 10, by = 2))+
  theme_classic(base_size = 22)+
  theme(legend.position = "none")

#Brooding day
biparental_brooding<-biparental_brooding%>%
  filter(sex!="total_f&m")%>%
  mutate(brooding_daily_min=brooding_daily_h*60)

str(biparental_brooding)
qqnorm(biparental_brooding$brooding_daily_min)
qqline(biparental_brooding$brooding_daily_min, lty=2)
shapiro.test(biparental_brooding$brooding_day_min) # iF P<0.05 data deviate from normality

mean(biparental_brooding$brooding_daily_min)
var (biparental_brooding$brooding_daily_min) 

ggplot(biparental_brooding, aes(y=brooding_daily_min, x=true_day, colour=factor(sex),shape=factor(sex)))+
  geom_point(aes(y=brooding_daily_min, x=true_day), size=4)+
  scale_shape_manual(values=c(1,19))+
  geom_smooth(method = "glm",method.args=list(family ="quasipoisson" (link = "log")), fill="grey")+
  scale_color_manual(breaks = c("male", "female"),values=c("black", "darkgrey"))+
  geom_rug()+
  labs(title="",y = " Daytime brooding (min/day)", x = "Days after hatching")+
  scale_x_continuous (limits = c(0, 12), breaks = seq(0, 11, by = 2))+
  theme_classic(base_size = 22)+
  theme(legend.position = "none")


#Foraging trips

qqnorm(biparental_trips$trips)
qqline(biparental_trips$trips, lty=2)
shapiro.test(biparental_trips$trips) # iF P<0.05 data deviate from normality

mean(biparental_trips$trips)
var (biparental_trips$trips)

biparental_trips<-biparental_trips%>%
  filter(sex!="total_f&m") 

str(biparental_trips)
View(biparental_trips)

ggplot(biparental_trips, aes(y=trips, x=true_inc_day, colour=factor(sex),shape=factor(sex)))+
  geom_point(aes(y=trips, x=true_inc_day), size=4)+
  scale_shape_manual(values=c(1,19))+
  geom_smooth(method = "glm",method.args=list(family ="poisson" (link = "log")), fill="grey")+
  scale_color_manual(breaks = c("male", "female"),values=c("black", "darkgrey"))+
  geom_rug()+
  labs(title="", y = "Trips frequency (#trips/day)", x = "Incubation day")+
  theme_classic(base_size = 22)+
  theme(legend.position = "none")

#Feeding
biparental_feeding<-biparental_feeding%>%
  filter(sex!="total_f&m") 

qqnorm(biparental_feeding$feeding_events)
qqline(biparental_feeding$feeding_events, lty=2)
shapiro.test(biparental_feeding$feeding_events) # iF P<0.05 data deviate from normality

mean(biparental_trips$trips)
var (biparental_trips$trips)

biparental_feeding<-biparental_feeding%>%
  filter(sex!="total_f&m") 

str(biparental_feeding)
View(biparental_feeding)

update_geom_defaults("smooth", list(size = 1))


ggplot(biparental_feeding, aes(y=feeding_events, x=true_day, colour=factor(sex),shape=factor(sex)))+
  geom_point(aes(y=feeding_events, x=true_day), size=4)+
  scale_shape_manual(values=c(1,19))+
  stat_smooth(method = "lm",formula=y~x,fill="grey")+
  scale_color_manual(breaks = c("male", "female"),values=c("black", "darkgrey"))+
  geom_rug()+
  labs(title="",y = "Feeding frequency (# events/day)", x = "Days after hatching")+
  scale_x_continuous (limits = c(1, 10), breaks = seq(1, 10, by = 2))+
  theme_classic(base_size = 22)+
  theme(legend.position = "none")


# Growth rate nestlings ---------------------------------------------------
##  Calculating K (Growth rate constant)
## Modified version of Ricklefs from Juan Pablo Gomez
## R-code  with annotations from Jenny Munoz
#### last update: June 26 2018
###############################################################################
#Run the functions 
# Logistic function
exact.logistic<- function(k,b,times,no){
  n.t <- b/(1+((b-no)/no)*exp(-k*times));
  return(n.t)
}

#Function gompertz
exact.gompertz<- function(k,b,times,no){
  n.t <- b*exp(-no*exp(-k*times));
  return(n.t)
}

#Function  von bertalanffy
exact.bert<- function(k,b,times){
  
  n.t <- (b*((1-exp(-k*times)^0.75)));
  return(n.t)
}

#Function to calculate the model probability
model.fit<-function(par,times,weights,model=c("logisitc","gompertz","vonbert")){
  k=par[1]
  b=par[2]
  no=weights[1]
  if(model=="logistic"){
    estimated<-exact.logistic(k,b,times,no)}
  else if (model=="gompertz"){
    estimated	<-	exact.gompertz(k,b,times,no)
  }else if(model=="vonbert"){
    estimated	<-	exact.bert(k,b,times)
  }
  ssq		<-	sum((weights-estimated)^2)
  return(ssq)
}

#  Functionto calculare AIC from  minimos cuadrados

AIC.lssq <- function(n,p,trss){
  
  aic <- n*(log((2*pi*trss)/n)+1) + 2*(p+1)
  return(aic)
}


# Funcion para encontrar b y k a partir de minimos cuadrados
masa  <-	c(1.56,3.21,4.17,4.97,6.26,8,8.26,8.14 ) ## mass data nestling 1
masa1  <-	c(1.31,2.79,3.78,4.53,5.85,6.67,7.57,8.5 ) ## mass data nestling 2
dias 	<-	c(1,2,3,4,5,7,8,9)## days data nestling 1
dias  	<-c(1,2,3,4,5,7,8,9) ## days data nestling 2

SSQmasa		<-	sum((sum(masa)-mean(masa))^2)
logist.mod	<-	optim(par=c(1,1),fn=model.fit,times=dias,weights=masa,model="logistic")
logit.k		<-	logist.mod$par[1]
logit.b		<-	logist.mod$par[2]
r2.logit	<-	1-logist.mod$value/SSQmasa
AIC.logit	<-	AIC.lssq(length(masa),2,logist.mod$value)
gompertz.mod<-	optim(par=c(1,1),fn=model.fit,times=dias,weights=masa,model="gompertz")
gom.k		<-	gompertz.mod$par[1]
gom.b		<-	gompertz.mod$par[2]
r2.gomp		<-	1-gompertz.mod$value/SSQmasa
AIC.gomp	<-	AIC.lssq(length(masa),2,gompertz.mod$value)
bert.mod	<-	optim(par=c(1,1),fn=model.fit,times=dias,weights=masa,model="vonbert")
bert.k		<-	bert.mod$par[1]
bert.b		<-	bert.mod$par[2]
r2.bert		<-	1-bert.mod$value/SSQmasa
AIC.bert	<-	AIC.lssq(length(masa),2,bert.mod$value)

results	<-	data.frame(model=c("Logistic","Gompertz","VonBertalanfy"),
                      b=round(c(logit.b,gom.b,bert.b),3),
                      k=round(c(logit.k,gom.k,bert.k),3),
                      r2=round(c(r2.logit,r2.gomp,r2.bert),3),
                      AIC=round(c(AIC.logit,AIC.gomp,AIC.bert),2))

par(mar=c(4,4,1,1),mgp=c(2.5,0.75,0))
plot(dias,masa,pch=19,cex=1.5,xlab="Times (Days)",ylab="Weight (g)",cex.lab=1.3,bty="l",xaxt="n",yaxt="n")
axis(1,seq(0,10,by=1),lwd=0,lwd.ticks=1.5,cex.axis=1.2)
axis(2,seq(0,10,by=1),lwd=0,lwd.ticks=1.5,cex.axis=1.2)
curve(exact.logistic(logit.k,logit.b,x,1.49),from=1,to=11,col="red",add=T,lwd=2)
curve(exact.gompertz(gom.k,gom.b,x,1.49),from=1,to=11,col="blue",add=T,lwd=2)
curve(exact.bert(bert.k,bert.b,x),from=1,to=11,col="green",add=T,lwd=2)
legend("bottomright",legend=c("Observed Data","Logistic Model","Gompertz","Von Bertalanffy"),col=c("black","red","blue","green")
       ,pch=c(19,-1,-1,-1),lty=c(-1,1,1,1),bty="n",cex=1.3,lwd=2)

results

####### only logistic model

par(mar=c(4,4,1,1),mgp=c(2.5,0.75,0))
plot(dias,masa,pch=19,cex=3,xlab="Times (Days)",ylab="Weight (g)",cex.lab=2,bty="l",xaxt="n",yaxt="n",col="black")
axis(1,seq(0,10,by=1),lwd=0,lwd.ticks=2,cex.axis=2)
axis(2,seq(0,10,by=1),lwd=0,lwd.ticks=2,cex.axis=2)
curve(exact.logistic(logit.k,logit.b,x,1.49),from=1,to=11,col="black",add=T,lwd=2)

par=TRUE

par(mar=c(4,4,1,1),mgp=c(2.5,0.75,0))
plot(dias,masa1,pch=19,cex=3,xlab="Times (Days)",ylab="Weight (g)",cex.lab=2,bty="l",xaxt="n",yaxt="n",col="grey")
axis(1,seq(0,10,by=1),lwd=0,lwd.ticks=2,cex.axis=2)
axis(2,seq(0,10,by=1),lwd=0,lwd.ticks=2,cex.axis=2)
curve(exact.logistic(logit.k,logit.b,x,1.49),from=1,to=11,col="grey",add=T,lwd=2)

plot(dias~masa,pch=19)

plot.default(masa~dias)

################Growth rate plot

#growthrate<-read.csv(file="Isleria_hauxwelli_growth_rate.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
growthrate<-read.csv(file="Isleria_hauxwelli_growth_rate1.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))

head(growthrate)
summary(growthrate)
str(growthrate)

#Exploring for missing data
which(is.na(growthrate)==T)

#Exploring the data
summary(growthrate) ## summary statistics for each column
unique(growthrate$date) ## unique values in a column

## the `skimr` package provides a really useful function for doing many of 
## these data exploration steps (and so much more!) in one simple command
skimr::skim(growthrate)

library(tidyverse)
growthrate<-growthrate%>%
  filter(type!="mass")

str(growthrate$type)
str(growthrate)         


shape=factor(growthrate$Individual)

str(growthrate)
ggplot(growthrate, aes(y=measurements, x=True_date, colour=factor(type)))+
  geom_point (aes(y=measurements, x=True_date,shape=factor(Individual)),size=6)+
  scale_shape_manual(values=c(20,1,0,15))+
  scale_colour_manual(values=c("cornflowerblue","aquamarine4"))+
  stat_smooth(method = "glm", method.args=list(family ="gaussian" (link = "log")),fill="grey60")+
  theme_classic(base_size = 22)+
  xlim(0,10)+
  labs(title="",y = "measurement (mm)", x = "Days after hatching")+
  scale_x_continuous (limits = c(1, 12), breaks = seq(0, 12, by = 2))+
  theme(legend.position = "right")

ggplot(growthrate, aes(y=measurements, x=True_date, colour=factor(type)))+
  geom_point (aes(y=measurements, x=True_date,shape=factor(Individual)),size=6)+
  scale_shape_manual(values=c(20,1,0,15))+
  scale_colour_manual(values=c("grey55","black"))+
  stat_smooth(method = "glm", method.args=list(family ="gaussian" (link = "log")),fill="grey60")+
  theme_classic(base_size = 22)+
  xlim(0,10)+
  labs(title="",y = "Measurement (mm)", x = "Days after hatching")+
  scale_x_continuous (limits = c(0, 10), breaks = seq(0, 10, by = 2))+
  theme(legend.position = "right")


#models
lm <- lm(tarsus ~ True_date +(1|Nestling), data = growthrate)
lm1 <- lm(wing ~ True_date +(1|Nestling), data = growthrate)
lm2 <- lm(mass ~ True_date +(1|Nestling), data = growthrate)

#Fig2_Plot nest attentiveness

nest_attentiveness<-read.csv(file="Nest_attentiveness_by_method.csv", stringsAsFactors=FALSE, strip.white=TRUE, na.strings=c("NA",""))
str(nest_attentiveness)

nest_attentiveness<-nest_attentiveness %>% 
  filter(true_inc_date!=22)

nest_attentiveness(true_inc_date)


ggplot(nest_attentiveness, aes(y=daily_nest_attentiveness, x=true_inc_date, colour=factor(method)))+
  geom_point (aes(y=daily_nest_attentiveness, x=true_inc_date, colour=factor(method)),size=6)+
  scale_shape_manual(values=c(0,16))+
  scale_colour_manual(values=c("black","grey55"))+
  stat_smooth(method= "glm",fill="grey", se=FALSE)+ #se=FALSE removes teh confidence intervals that represend the s.errror
  theme_classic(base_size = 22)+
  scale_x_continuous (limits = c(0, 25), breaks = seq(0, 25, by = 2))+
  scale_y_continuous (limits = c(0, 100), breaks = seq(0, 100, by = 10))+
  labs(title="",y = "Nest attentiveness (%)", x = "Incubation day")+
  theme(legend.position = "none")

#camera grey #sensor black


ggplot(biparental_incubation, aes(y=daily_nest_attentiveness, x=true_inc_date, colour=factor(method)))+
  geom_point (aes(y=daily_nest_attentiveness, x=true_inc_date, colour=factor(method)),size=6)+
  scale_shape_manual(values=c(0,16))+
  scale_colour_manual(values=c("black","coral2"))+
  stat_smooth(method= "glm",fill="grey", se=FALSE)+ #se=FALSE removes teh confidence intervals that represend the s.errror
  theme_classic(base_size = 22)+
  scale_x_continuous (limits = c(0, 25), breaks = seq(0, 25, by = 2))+
  scale_y_continuous (limits = c(0, 100), breaks = seq(0, 100, by = 10))+
  labs(title="",y = "Nest attentiveness (%)", x = "Incubation day")+
  theme(legend.position = "none")



#End of the script



ggsave()


#End of script Manuscript Isleria hauxwelli