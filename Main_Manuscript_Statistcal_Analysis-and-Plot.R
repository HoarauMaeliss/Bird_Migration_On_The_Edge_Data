library(lubridate)
library(tidyverse)
library(nlme) ## Non-linear mied effect
library(lme4)

############### Behavior ################
MSA_raw <- get(load("C:/Users/maeli/OneDrive/Maitrise_Canada/Big_Mammelle/PNAS_submission/Bird-Migration-On-The-Edge_Data/Statistical_Analysis/Donne_accelerometries_IAO_2021_3sec.RData"))

## Regrouping clusters into behaviors:
MSA_behavior <- MSA_raw%>%
  mutate(foraging=case_when(
    EMcluster5==3~1,
    EMcluster5==2~0,
    TRUE~0),
    forandother= case_when(
      EMcluster5==3~1,
      EMcluster5==2~1,
      TRUE~0),
    other= case_when(
      EMcluster5==3~0,
      EMcluster5==2~1,
      TRUE~0))



## Removing pair number 9 because HJ was injured during the capture
MSA_day<-MSA_behavior%>%
  filter(paire!=9)

# # Add number of observation per ID, day, periode
data_MSA<- MSA_day%>%
  group_by(ID,JT,periode, trt, paire,numero_capture)%>%
  summarise(nb_obs= n(),
            prop_foraging = sum(foraging)/n(),
            prop_all= sum(forandother)/n(),
            prop_other = sum(other)/n())


#Filter, relevel and clean dataset
data_MSA$scale_JT <- data_MSA$JT
data_MSA <- data_MSA[,c("prop_foraging","prop_all", "prop_other", "trt", "periode", "scale_JT", "paire", "ID", "nb_obs")]
data_MSA <- data_MSA %>% mutate(trt=factor(trt, levels = c("placebo", "cort")),
                                periode = factor(periode, levels = c("day", "evening", "morning")),
                                paire = factor(paire),
                                ID = factor(ID))


### Table 1. Generalized linear mixed model (GLMM) for the proportion of time spent foraging...
GLM1 <- MASS::glmmPQL(prop_foraging~ periode + trt*scale_JT, random = list(~1|paire, ~1|ID), family = quasibinomial, weight=nb_obs, data=data_MSA[data_MSA$scale_JT>1 ,], na.action = na.omit)
summary(GLM1)


newdata <- data.frame(expand.grid(trt = levels(data_MSA$trt),
                                  periode= "day",
                                  scale_JT = seq(from= 2, to= 10, by=1)))

### Fig.1. Foraging intensity of greater snow goose....
library(ggridges)
library(ggnewscale)

A <- ggplot(data_MSA[data_MSA$periode == "day",], aes(y = as.factor(scale_JT), x = prop_foraging, fill = trt)) +
  stat_summary(aes(col=trt), fun.data = "mean_se", geom = "pointrange", show.legend = F,
               position =  position_dodge(0.25, preserve = 'single'), alpha=1, size=1.2, stroke=2.15, colour="black", pch=21) +  
  scale_fill_manual(values = c("#F3977F", "#A9CDDB"), breaks = c("cort", "placebo"), labels= c("CORT", "PLACEBO"))+
  
  new_scale_fill() +
  geom_density_ridges2(data=data_MSA[data_MSA$periode == "day",], aes(y = as.factor(scale_JT), x = prop_foraging, fill = trt),
                       scale =0.7, show.legend = TRUE, alpha=0.4, bandwidth = 0.045, col= NA)+
  
  scale_fill_manual(values = c("#E73000", "#529CB5"), breaks = c("cort", "placebo"), labels= c("CORT", "PLACEBO"))+
  
  xlab("Proportion of time foraging")+
  ylab("Days after pellet implantation")+
  theme_minimal()+ 
  coord_flip(xlim= c(0,1))+
  
  theme(panel.background = element_rect(fill='white'), #transparent panel bg
        plot.background = element_rect(fill='white'),
        plot.title = element_text(size=22, face="bold", hjust = -0.3),
        legend.title=element_blank(),
        legend.text = element_text(size=14, colour="black"),
        axis.text=element_text(size=17, colour="black"),
        axis.title=element_text(margin =margin(0,20,20,0), size=18,face="bold"),
        axis.line = element_line(color = 'black'))
A


  

############### Departure Date ################
a= read.csv("C:/Users/maeli/OneDrive/Maitrise_Canada/Big_Mammelle/PNAS_submission/Bird-Migration-On-The-Edge_Data/Statistical_Analysis/departure_date.csv", sep=",", header = T)

# Filter pair that are dead killed or with defective collar
'%!in%' <- function(x,y)!('%in%'(x,y))

a<- a[a$paire %!in% c(9,10,11,22,25,26), ]

a$departure_date<- a$julian_departure

a$departure_date<- as.numeric(a$departure_date)

a$date_deploiement<- as.numeric(format(as.POSIXct(strptime(a$Deploy.On.Timestamp , "%Y-%m-%d")),"%j"))


a<- a%>%select(trt,paire,departure_date,date_deploiement)

a<- spread(a,key= trt, value= departure_date)
t.test(a$cort, a$placebo, paired = TRUE)


#Check t-test assumptions (normality of residuals)
d_tmp<-a$cort - a$placebo
d<-d_tmp-mean(d_tmp) #Residuals for paired t-test: differnece between pairwise differnece and mean of those differences
plot(density(d)) #not perfect,
shapiro.test(d)  #but not significantly different than normality
#Another way to visually asses fit to normal distribution
qqnorm(d)
qqline(d, datax = FALSE, distribution = qnorm, probs = c(0.25, 0.75)) #Seems ok
#One more check: similar results obtained through linear regression with pair as random effect?
check_glmer_resid<-lmer(data=a_lmer, departure_date~trt+(1|paire))
#Check normality of residuals and homogeneity of variances
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = check_glmer_resid)
plot(simulationOutput) # not so bad here
summary(check_glmer_resid) #Effect is the same
confint(check_glmer_resid) #95% C.I. for treatment effect [0.08,2.31], fairly equivalent to p=0.04
#We conclude the paired t-test is appropriate

### Fig. 2.Departure date for migration of female greater snow goose ...
a= read.csv("C:/Users/maeli/OneDrive/Maitrise_Canada/Big_Mammelle/PNAS_submission/Bird-Migration-On-The-Edge_Data/Statistical_Analysis/departure_date.csv", sep=",", header = T)

# Filter pair that are dead killed or with defective collar
'%!in%' <- function(x,y)!('%in%'(x,y))

a<- a[a$paire %!in% c(9,10,11,22,25,26), ]

a$departure_date<- a$julian_departure

a<- a%>%select(trt,paire,departure_date,id)
a$departure_date<- as.numeric(a$departure_date)


a<-a%>%
  group_by(trt,paire,id)%>%
  mutate(presence=1)%>%
  pivot_wider(names_from = departure_date, values_from = presence, values_fill = 0)%>%
  mutate('132'=0,
         '133'=0,
         '136'=0,
         '142'=0,
         '143'=0,
         '144'=0)

a<-a%>%
  group_by(trt,paire,id)%>%
  gather(key = "departure_date", value = "depart",
       4:16)
a$trt<- as.factor(a$trt)
a$departure_date <- as.numeric(a$departure_date)

library(lme4)
m1<- lme4::glmer(data=a,depart ~departure_date*trt+(1|paire), family = binomial)
summary(m1)


find_mode <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

a %>% group_by(trt) %>% summarise(find_mode(departure_date))

S<- a%>%
  group_by(trt)%>%
  summarise(m=mean(departure_date),
            SD= sd(departure_date),
            CI_L= m -(SD*1.96)/sqrt(50),
            CI_U= m + (SD*1.96)/sqrt(50))

a<- a%>%mutate(trt=factor(trt, levels = c("placebo", "cort")))

library(ggplot2)
library(ggridges)
### COULEUR 
g<- ggplot()+
  #geom_violin(mapping= aes(x=trt, y=departure_date,fill= trt), data= a,alpha=0.4)+
  geom_density_ridges2(mapping= aes(y=trt, x=departure_date,fill= trt), points=1, data= a, scale=0.9, bandwidth =0.9)+
  geom_dotplot(aes(y=trt, x=departure_date,fill= trt), data= a[a$trt=="placebo",], position = position_nudge(y=1), fill="gray45")+
  geom_dotplot(aes(y=trt, x=departure_date,fill= trt), data= a[a$trt=="cort",], position = position_nudge(y=2), fill="gray45")+
  
  geom_errorbar(mapping = aes(y= trt, xmin= CI_L, xmax= CI_U), data = S, width= 0.1, size=1.2, position= position_nudge(x = 0, y = -0.1))+
  geom_point(mapping = aes(y=trt, x=m, fill= trt),data= S, size=5, position= position_nudge(x = 0, y = -0.1), pch=21, stroke=2.2)+
  
  scale_x_continuous(n.breaks = 10)+
  scale_y_discrete(labels = c("PLACEBO", "CORT"))+
  scale_fill_manual(values = c("#F3977F", "#A9CDDB"), breaks = c("cort", "placebo"))+
  
  labs(y= "", x = "Departure Date in May")+
  
  theme_classic()+
  theme(text = element_text(size = 20),
        legend.position = "none")
g



############### Habitat use ################
library(sf)
library(raster)
library(mapview)

load("C:/Users/maeli/OneDrive/Maitrise_Canada/Big_Mammelle/PNAS_submission/Bird-Migration-On-The-Edge_Data/Statistical_Analysis/gps_with_landclass.RData")


#remove all geometry and spatial structure now that we don't need it
exp_std_all<-habitat_analysis_class_JT

st_geometry(exp_std_all)<-NULL


# Filter all point that are in the St Laurent region
exp_std_stlo<-exp_std_all%>%
  filter(in_stlo=='IN')

#This restricts data for each pair to days when both members of the pair are in the St-Lawrence valley
exp_std_restr<-exp_std_all%>%
  filter(jday<=last_stlo,
         jday>=start_stlo)


#Select which dataset to analyse
dataset<-exp_std_stlo


#### St-Lawrence data #####
#Compute number of points taken every day
fixperday<- dataset%>%
  group_by(ID, jday)%>%
  dplyr::summarize(ntot=n())%>%ungroup()

#Compute proportion of points that are in each habitat 
prop_time<-dataset%>%
  left_join(fixperday, by=c('ID','jday'))%>%
  group_by(ID, jday, JT, class_agri, trt, paire)%>% #Grouping by class allows computing thenumber of points in each habitat class
  dplyr::summarize(pct=(n()/ntot), ntot=mean(ntot))%>%
  print(n=26)

#Returned dataframe contains multiple duplicates (all points are kept when we only want one point per habitat per day with its corresponding percentage)
#remove duplicates
prop_time<-unique(prop_time)

#Keep only proportion of time in cropland, which is the variable we're interested in, and only keep days of treatment 2:10
prop_crop<-prop_time%>%
  filter(class_agri=='cropland')%>%
  filter(JT>=2,JT<=10)

#Run the analysis
mod_pql_all<-MASS::glmmPQL(pct~trt, random = list(~1|paire, ~1|ID), family='quasibinomial', weights=ntot, data=prop_crop)
summary(mod_pql_all)



