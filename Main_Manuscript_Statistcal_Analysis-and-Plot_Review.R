library(lubridate)
library(tidyverse)
library(nlme) ## Non-linear mied effect
library(lme4)
library(ggplot2)

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
data_MSA <- data_MSA[,c("prop_foraging","prop_all", "prop_other", "trt", "periode", "scale_JT", "paire","numero_capture", "ID", "nb_obs")]
data_MSA <- data_MSA %>% mutate(trt=factor(trt, levels = c("placebo", "cort")),
                                periode = factor(periode, levels = c("day", "evening", "morning")),
                                paire = factor(paire),
                                ID = factor(ID))


### Table 1. Generalized linear mixed model (GLMM) for the proportion of time spent foraging...
GLM1 <- MASS::glmmPQL(prop_foraging~ periode + trt*scale_JT, random = list(~1|paire, ~1|ID), family = quasibinomial, weight=nb_obs, data=data_MSA[data_MSA$scale_JT>1 ,], na.action = na.omit)
summary(GLM1)

#### Appendix S1: Table S1
GLM1 <- MASS::glmmPQL(prop_all~ periode + trt*scale_JT, random = list(~1|paire, ~1|ID), family = quasibinomial, weight=nb_obs, data=data_MSA[data_MSA$scale_JT>1 ,], na.action = na.omit)
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
  
  #scale_fill_manual(values = c("#E73000", "#529CB5"), breaks = c("cort", "placebo"), labels= c("CORT", "PLACEBO"))+
  
  xlab("Proportion of time foraging")+
  ylab("Days after pellet implantation")+
  theme_minimal()+ 
  coord_flip(xlim= c(0,1))+
  scale_color_manual(values=c("#E73000","#529CB5"),labels=c('CORT: treated with 90mg of corticosterone in a cholesterol matrix','placebo: treated cholesterol matrix only (contol)'))+
  scale_fill_manual(values=c("#E73000","#529CB5"),labels=c('CORT: treated with 90mg of corticosterone in a cholesterol matrix','placebo: treated cholesterol matrix only (contol)'),
                    guide = guide_legend(override.aes = list(color = NA)))+
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
a= read.csv("C:/Users/maeli/OneDrive/Maitrise_Canada/Big_Mammelle/Ecology_submission/Re_Submission/Bird-Migration-On-The-Edge_Data/Statistical_Analysis/departure_date.csv", sep=",", header = T)

## Calculate mean de parture dates
# Filter pair that are dead killed or with defective collar
'%!in%' <- function(x,y)!('%in%'(x,y))
a<- a[a$paire %!in% c(9,10,11,22,25,26), ]


# formating variable
a$depart_date<- as.numeric(format(as.POSIXct(strptime(a$depart_date , "%Y-%m-%d")),"%d"))
a$departure_date<- as.numeric(a$julian_departure)
a$date_deploiement<- as.numeric(format(as.POSIXct(strptime(a$Deploy.On.Timestamp , "%Y-%m-%d")),"%j"))
a$trt<- factor(a$trt, levels = c("cort", "placebo"))
a$date_capt<- a$date_deploiement-(min(a$date_deploiement)-1)

## Calculate the median departure date since the distribution of the data is an inverse skweness
# Define a bootstrapping function for the median
bootstrap_median <- function(data, indices) {
  sample_data <- data[indices]
  return(median(sample_data))
}

# Apply bootstrapping to calculate median and CI for each group
install.packages("boot")
library(boot)
result <- a %>%
  group_by(trt) %>%
  do({
    boot_result <- boot(
      data = .$depart_date, 
      statistic = bootstrap_median, 
      R = 1000  # Number of bootstrap samples
    )
    boot_ci <- boot.ci(boot_result, type = "perc")  # Percentile CI
    
    data.frame(
      trt = unique(.$trt),
      median = median(.$depart_date),
      CI_L = boot_ci$percent[4],  # 2.5th percentile
      CI_U = boot_ci$percent[5]   # 97.5th percentile
    )
  })



####### Data frame conversion to one line per day with a column is_departed
a<- a %>% dplyr::select(-c(Notes,RECAP.BYLOT,Deploy.Off.Timestamp, julian_departure))
colnames_a <- colnames(a)


# convert dataframe to one line per day with a column is_departed (0/1)
instance_df <- data.frame()

for(i in c(1:nrow(a))){
  current_line <- a[i,]
  time_tracked <- current_line$departure_date-current_line$date_deploiement
  
  all_instance_lines <- current_line[rep(1, each = time_tracked), ]
  all_instance_lines$day_since_trt = c(1:time_tracked)
  all_instance_lines$current_day = all_instance_lines$date_deploiement + all_instance_lines$day_since_trt
  
  all_instance_lines$is_departed <- with(all_instance_lines, case_when(current_day < departure_date ~ 0,
                                                                       current_day == departure_date ~ 1,
                                                                       current_day > departure_date  ~ NA))
  
  instance_df <- rbind(instance_df, all_instance_lines)
}


library(coxme)
cme <- coxme(Surv(day_since_trt, is_departed) ~ trt*date_capt + (1|paire), data= instance_df, x=T)

###### Graphics figure 2 and Annexe #####
#Compute coefficients table for model with interaction
(est_plc <- as.numeric(exp(cme$coefficients)))
(sd <- confint(cme))
coef_table<-as.data.frame(cbind(cme$coefficients,sd))
colnames(coef_table)<-c('Estimate','LowerCI','UpperCI')
coef_table


#for-loop to generate figure, put the capture numbers for which you want to generate the figure in 'selectedcapdates'
selectedcapdates<-sort(unique(instance_df$date_capt))

for(capdat in selectedcapdates){
  # build dummy dataset
  pred_data = data.frame(expand.grid(trt=unique(instance_df$trt),
                                     date_capt=capdat))
  #0.6517440 -0.5920422 
  
  # get predicted values
  # library(devtools)
  # install_github('junkka/ehahelper')
  library(ehahelper)
  preds <- ehahelper::predict_coxme(cme, pred_data, se.fit=TRUE, type = "lp")
  pred_data$est <- preds$fit
  pred_data$upperCI <- preds$fit + 1.96*preds$se.fit
  pred_data$lowerCI <- preds$fit - 1.96*preds$se.fit
  
 
  breslow_est_adj_inter <- function(time, status, X, B, fit_cox){
    data <- data.frame(time,status,X)
    data <- data[order(data$time), ]
    t    <-  unique(data$time)
    k    <- length(t)
    h    <- rep(0,k)
    #[]
    # mean(data$trtplacebo)
    # mean(data$date_capt)
    # mean(data$date_capt)
    
    
    for(i in 1:k) {
      
      LP_sample <- sum(fit_cox$means * coef(fit_cox)) # I'm not sure how this should be different for our interaction model, need to check how this works
      
      #Individual linear predictor for interaction model
      LP_indiv <- (coef(fit_cox)['trtplacebo']*(data$trtplacebo))+(coef(fit_cox)['date_capt']*(data$date_capt))+(coef(fit_cox)['trtplacebo:date_capt']*(data$trtplacebo)*data$date_capt) 
      
      lp_centered <- (LP_indiv - LP_sample)[data$time>=t[i]]
      risk <- exp(lp_centered)
      h[i] <- sum(data$status[data$time==t[i]]) / sum(risk)
    }
    
    res <- cumsum(h)
    return(res)
  }
  
  
  # test the baseline hazard function by comparing the hazard obtained with basehaz, 
  # with breslow estimator for a coxph model and with breslow with a coxme model 
  cph <- coxph(Surv(day_since_trt, is_departed) ~ trt*date_capt, data= instance_df, x=T)
  cme <- coxme(Surv(day_since_trt, is_departed) ~ trt*date_capt + (1|paire), data= instance_df, x=T)

  
  H0_basehaz_cph = basehaz(cph,centered = T)
  H0_basehaz_cph$breslow_cme <- breslow_est_adj_inter(time=instance_df$day_since_trt, status=instance_df$is_departed, X=cme$x, B=cme$coefficients, fit_cox = cme)
  
  
  # building dataset with confidence interval
  bh <- H0_basehaz_cph[, c("breslow_cme", "time")]
  colnames(bh) <- c("hazard","day_since_trt")
  
  #Just a little trick so that CIs appear until last day of treatment on figure
  bh<-rbind(bh,bh[22,])
  bh[23,2]<-23
  
  bh_pred <- rbind(bh, bh)
  bh_pred$trt <- rep(c("placebo", "cort"), each = nrow(bh))
  
  
  ### calculate coxme estimate using the model function: h(t)=h0(t)exp(b1X1) with estimates obtained with coxme models
  bh_pred <- bh_pred %>%
    mutate(est= case_when(
      trt=="placebo" ~ 1-exp(-hazard*exp(as.vector(pred_data[pred_data$trt == "placebo", "est"]))),
      trt=="cort" ~ 1-exp(-hazard*exp(as.vector(pred_data[pred_data$trt == "cort", "est"])))
    ),
    sd_plus= case_when(
      trt=="placebo" ~ 1-exp(-hazard*exp(as.vector(pred_data[pred_data$trt == "placebo", "upperCI"]))),
      trt=="cort" ~ 1-exp(-hazard*exp(as.vector(pred_data[pred_data$trt == "cort", "upperCI"])))
    ),
    sd_moins= case_when(
      trt=="placebo" ~ 1-exp(-hazard*exp(as.vector(pred_data[pred_data$trt == "placebo", "lowerCI"]))),
      trt=="cort" ~ 1-exp(-hazard*exp(as.vector(pred_data[pred_data$trt == "cort", "lowerCI"])))
    ))
  
  
  bh_pred_temp<-bh_pred%>%mutate(capt=capdat)
  
  #Create dataframe for figure with all predictions
  if(capdat==selectedcapdates[1]){
    
    bh_pred_fig<-bh_pred_temp
    
  }else{
    
    bh_pred_fig<-bh_pred_fig%>%
      bind_rows(bh_pred_temp)
    
  }
  
}




#Graphic FIG 2.
library(ggridges)
a<- a%>%mutate(trt=factor(trt, levels = c("placebo", "cort")))

A<- ggplot()+
  #geom_violin(mapping= aes(x=trt, y=departure_date,fill= trt), data= a,alpha=0.4)+
  ggridges::geom_density_ridges2(mapping= aes(y=trt, x=depart_date,fill= trt), data= a, scale=0.9, bandwidth =0.9)+
  geom_dotplot(aes(y=trt, x=depart_date,fill= trt), data= a[a$trt=="placebo",], position = position_nudge(y=1), fill="gray45")+
  geom_dotplot(aes(y=trt, x=depart_date,fill= trt), data= a[a$trt=="cort",], position = position_nudge(y=2), fill="gray45")+
  
  scale_x_continuous(n.breaks = 10)+
  scale_y_discrete(labels = c("PLACEBO", "CORT"))+
  scale_fill_manual(values = c("#F3977F", "#A9CDDB"), breaks = c("cort", "placebo"))+
  
  geom_errorbar(mapping = aes(y= trt, xmin= CI_L, xmax= CI_U), data = result, width= 0.1, size=1.2, position= position_nudge(x = 0, y = -0.1))+
  geom_point(mapping = aes(y=trt, x=median, fill= trt),data= result, size=5, position= position_nudge(x = 0, y = -0.1), pch=21, stroke=2.2)+
  
  
  labs(y= "", x = "Departure date in May")+
  
  theme_classic()+
  theme(text = element_text(size = 20),
        legend.position = "none")
A

a %>%
  group_by(date_capt) %>%
  summarize(n=n(),
            deploy= unique(Deploy.On.Timestamp)) %>%
  print(n=100)

hist(a$date_capt,breaks=15)#### Deciding which day to present on FIG 2. http://127.0.0.1:46845/graphics/plot_zoom_png?width=1920&height=1009

bh_pred_fig <- bh_pred_fig %>%
  dplyr::rename(date_capt = capt)

data_dep <- merge(bh_pred_fig, instance_df, by= c("day_since_trt","trt","date_capt"))

library(ggridges)
B <-ggplot(bh_pred_fig%>%filter(date_capt%in%c(7),day_since_trt%in%c(7:19)) , aes(day_since_trt-0.5, est, color= trt)) +
  geom_step()+
  geom_point(aes(day_since_trt, est, color= trt))+
  pammtools::geom_stepribbon(aes(ymin = sd_plus, ymax = sd_moins, fill = trt), alpha = .3)+
  labs(y= "Migratory departure probability", x = "Day after pellet implantation", fill="", color="")+
  scale_x_continuous(n.breaks=12)+
  scale_color_manual(values=c("#E73000","#529CB5"),labels=c('CORT','PLACEBO'))+
  scale_fill_manual(values=c("#E73000","#529CB5"),labels=c('CORT','PLACEBO'),
                    guide = guide_legend(override.aes = list(color = NA)))+
  theme_classic()+
  theme(text = element_text(size = 20),
        legend.position='top')
B

bh_pred_fig%>%filter(date_capt%in%c(7),day_since_trt%in%c(14))

library(gridExtra)
grid.arrange(A, B, nrow = 1)

# TABLE 2. 
library(coxme)
instance_df <- instance_df %>% mutate(trt=factor(trt, levels = c("placebo", "cort")))## Relevel factors to obtain  CORT estimates
cme <- coxme(Surv(day_since_trt, is_departed) ~ trt*date_capt + (1|paire), data= instance_df, x=T)
summary(cme)
  
p <- sjPlot::tab_model(cme, dv.labels = "Probability of departure",transform=NULL)
p


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



