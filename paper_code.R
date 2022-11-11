############################################################################
############################################################################
##########################      SET UP       ###############################
############################################################################
############################################################################
##### Set up ####
# The data used in these analyses are publicly available at DIETFITS' 
# Open Science Framework repository.

# This code downloads automatically (no need to save the data file).

# Address correspondence about this code to adrian.sotom@tec.mx

# Note: Some graphs were edited in our reanalysis paper for readability.

# If you are only interested in corroborating our results and figures, just 
# click "Ctrl + Shift + Enter" and get a coffee (results will take several
# minutes).

# If you want to run the code line by line, make sure to use
# the same package versions. Run lines 35- 46 in this code to ensure
# the correct versions are running.

#Working directory setup
setwd("C:/Users/adria/Dropbox/UIEM/LEAD/Proyectos/DIETFITS")

# Since "easystats" is currently not available in CRAN, install manually if you
# don't already have it.

# install.packages("easystats", repos = "https://easystats.r-universe.dev")

# Now, confirm you have "pacman" installed. If you don't have it, 
# remove the "#" in the line below and press "Enter".
# install.packages("pacman") 

#Packages setup
pacman::p_load(dplyr,tidyr,ggstatsplot,readxl,tableone,easystats,dagitty,lme4,
               patchwork,MASS,see,qqplotr,bootStepAIC,performance,ggdag,multcomp,
               rpart,rpart.plot,gtools,broom,lmtest,visdat,report,beepr,Rmisc,
               parameters,ggcharts,conflicted,car,rattle,cvms,lavaan,
               missForest,mlogit,MLmetrics,beepr,readr,haven,dagitty,mediation)
library(mice)

#Solving duplicate functions conflicts
conflict_prefer("select","dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("arrange", "dplyr")

##### Data Upload ####
data <- read.csv(url("https://osf.io/ztysq/download"))

#Participant ID=2001 is excluded due to an implausible lipid value (TG).
data <- data %>% filter(study_id != "2001") 

############################################################################
############################################################################
######################      DATA FILTERING       ###########################
############################################################################
############################################################################
##### Filtering by time points ####
# Note, in the original dataset "Blue" = LCD and "Purple" = LFD.
data <- data %>% mutate(diet=ifelse(diet=="Blue","LCD","LFD"))

data <- data %>% mutate(timepoint=
                          case_when(redcap_event_name=="12_months_arm_1"~"12m",
                                    redcap_event_name=="6_months_arm_1"~"6m",
                                    redcap_event_name=="3_months_arm_1"~"3m",
                                    redcap_event_name=="baseline_arm_1"~"Baseline"))

data_bl <- data %>%  filter(redcap_event_name == "baseline_arm_1") %>% 
  select(study_id,baseline_ins,gluc,diet,timepoint,scr_gender,
         dxa_percentfat,weight_gcrc,lipid_trig_v2,lipid_ldl_v2,
         lipid_hdl_v2,MetSyn,cal,fat..,carb..,protein..,bmi,
         saturated_fat..,fiber.g,carb.g,fat.g,protein.g,
         GI_glucose,GL_glucose,added_sugars,sugar.cal,
         thirty_ins) %>% 
  mutate(bl_carb_cals = carb.g*4) %>% 
  mutate(bl_prot_cals = protein.g*4) %>% 
  mutate(bl_fat_cals = fat.g*9) %>% 
  mutate(bl_tgtohdl = lipid_trig_v2/lipid_hdl_v2) %>% 
  mutate(bl_ldlplushdl=lipid_hdl_v2+lipid_ldl_v2) %>% 
  mutate(weight_change_bl=ifelse(diet=="LCD",0,0)) %>% 
  mutate(across(where(is.numeric), round, 2)) %>%
  distinct()

data_3m <- data %>% filter(redcap_event_name == "3_months_arm_1") %>% 
  select(study_id,baseline_ins,gluc,diet,scr_gender,dxa_percentfat,
         weight_gcrc,lipid_trig_v2,lipid_ldl_v2,timepoint,
         lipid_hdl_v2,MetSyn,cal,fat..,carb..,protein..,bmi,
         saturated_fat..,fiber.g,carb.g,fat.g,protein.g,
         GI_glucose,GL_glucose,added_sugars,sugar.cal,thirty_ins) %>% 
  mutate(delta_weight_3m = (data_bl$weight_gcrc - weight_gcrc)) %>% 
  mutate(weight_change_3m = (weight_gcrc-data_bl$weight_gcrc)) %>% 
  mutate(delta_weight_3m_2 = ((data_bl$weight_gcrc - weight_gcrc)/data_bl$weight_gcrc)*100) %>% 
  mutate(qdelta_weight_3m = quantcut(delta_weight_3m,q=5)) %>% 
  mutate(bl_bmi = data_bl$bmi) %>% 
  mutate(bl_carb_cals = data_bl$bl_carb_cals) %>% 
  mutate(carb_cals_3m = carb.g*4) %>% 
  mutate(delta_carbcals_3m = carb_cals_3m-bl_carb_cals) %>% 
  mutate(z_delta_carbcals_3m = (mean(delta_carbcals_3m,na.rm=T)-delta_carbcals_3m)/sd(delta_carbcals_3m,na.rm=T)) %>% 
  mutate(qdelta_carbcals_3m = quantcut(delta_carbcals_3m,q=5)) %>% 
  mutate(addedsugarcals = added_sugars*4) %>% 
  mutate(bl_prot_cals = data_bl$bl_prot_cals) %>% 
  mutate(prot_cals_3m = protein.g*4) %>% 
  mutate(delta_protcals_3m = bl_prot_cals - prot_cals_3m) %>% 
  mutate(qdelta_protcals_3m = quantcut(delta_protcals_3m,q=5)) %>% 
  mutate(bl_fat_cals = data_bl$bl_fat_cals) %>% 
  mutate(fat_cals_3m = fat.g*4) %>% 
  mutate(delta_fatcals_3m = fat_cals_3m-bl_fat_cals) %>% 
  mutate(z_delta_fatcals_3m = (mean(delta_fatcals_3m,na.rm=T)-delta_fatcals_3m)/sd(delta_fatcals_3m,na.rm=T)) %>% 
  mutate(qdelta_fatcals_3m = quantcut(delta_fatcals_3m,q=5)) %>% 
  mutate(bl_GI = data_bl$GI_glucose) %>% 
  mutate(delta_GI_3m = bl_GI - GI_glucose) %>% 
  mutate(qdelta_GI_3m = quantcut(delta_GI_3m,q=5)) %>% 
  mutate(bl_GL = data_bl$GL_glucose) %>% 
  mutate(delta_GL_3m = GL_glucose - bl_GL) %>%
  mutate(z_delta_GL_3m = (mean(delta_GL_3m,na.rm=T)-delta_GL_3m)/sd(delta_GL_3m,na.rm=T)) %>% 
  mutate(qdelta_GL_3m  = quantcut(delta_GL_3m,q=5)) %>%  
  mutate(deltaGLq5 = delta_GL_3m <= -107) %>% 
  mutate(q3delta_GL_3m  = quantcut(delta_GL_3m,q=3)) %>%  
  mutate(deltaGLq3 = delta_GL_3m <= -86.4) %>% 
  mutate(q4delta_GL_3m  = quantcut(delta_GL_3m,q=4)) %>%  
  mutate(deltaGLq4 = delta_GL_3m <= -100) %>% 
  mutate(bl_cal = data_bl$cal) %>% 
  mutate(delta_cal_3m = cal - bl_cal) %>%
  mutate(z_delta_cals_3m = (mean(delta_cal_3m,na.rm=T)-delta_cal_3m)/sd(delta_cal_3m,na.rm=T)) %>%
  mutate(qdelta_cal_3m  = quantcut(delta_cal_3m,q=5)) %>% 
  mutate(tgtohdl_bl = data_bl$bl_tgtohdl) %>% 
  mutate(tgtohdl_3m = lipid_trig_v2/lipid_hdl_v2) %>% 
  mutate(delta_tgtohdl_3m = tgtohdl_3m-tgtohdl_bl) %>% 
  mutate(z_delta_tgtohdl_3m= (mean(delta_tgtohdl_3m,na.rm=T)-delta_tgtohdl_3m)/sd(delta_tgtohdl_3m,na.rm=T)) %>% 
  mutate(qdelta_tgtohdl = quantcut(delta_tgtohdl_3m,q=5)) %>% 
  mutate(bl_30ins = data_bl$thirty_ins) %>%
  mutate(qbl_30ins = quantcut(bl_30ins,q=5)) %>%
  mutate(bl30q5 = bl_30ins >= 127) %>% 
  mutate(q4bl_30ins = quantcut(bl_30ins,q=4)) %>%
  mutate(bl30q4 = bl_30ins >= 118) %>% 
  mutate(q3bl_30ins = quantcut(bl_30ins,q=3)) %>%
  mutate(bl30q3 = bl_30ins >= 103) %>% 
  mutate(bl_satfat = data_bl$saturated_fat..) %>% 
  mutate(delta_satfat = bl_satfat - saturated_fat..) %>%
  mutate(qdelta_satfat = quantcut(delta_satfat,q=5)) %>% 
  mutate(bl_added_sugars = data_bl$added_sugars) %>% 
  mutate(delta_addsugars = bl_added_sugars - added_sugars) %>%
  mutate(qdelta_addsugars = quantcut(delta_addsugars,q=5)) %>% 
  mutate(bl_sugarcal = data_bl$sugar.cal) %>% 
  mutate(delta_sugarcal = bl_sugarcal - sugar.cal) %>%
  mutate(bl_ldl = data_bl$lipid_ldl_v2) %>% 
  mutate(delta_ldl = bl_ldl - lipid_ldl_v2) %>%
  mutate(qdelta_ldl = quantcut(delta_ldl,q=5)) %>%
  mutate(bl_ldlplushdl = data_bl$bl_ldlplushdl) %>% 
  mutate(ldlplushdl_3m = lipid_ldl_v2 + lipid_hdl_v2) %>%
  mutate(delta_ldlplushdl = ldlplushdl_3m - bl_ldlplushdl) %>% 
  mutate(z_delta_ldlplushdl= (mean(delta_ldlplushdl,na.rm=T)-delta_ldlplushdl)/sd(delta_ldlplushdl,na.rm=T)) %>%
  mutate(qdelta_ldlplushdl_3m = quantcut(delta_ldlplushdl,q=5)) %>%
  mutate(q5insGLr=deltaGLq5==T & bl30q5==T) %>%
  mutate(q4insGLr=deltaGLq4==T & bl30q4==T) %>%
  mutate(q3insGLr=deltaGLq3==T & bl30q4==T) %>%
  mutate(across(where(is.numeric), round, 2)) %>%
  distinct() 

data_6m <- data %>% filter(redcap_event_name == "6_months_arm_1") %>% 
  select(study_id,baseline_ins,gluc,diet,scr_gender,dxa_percentfat,
         weight_gcrc,lipid_trig_v2,lipid_ldl_v2,timepoint,
         lipid_hdl_v2,MetSyn,cal,fat..,carb..,protein..,bmi,
         saturated_fat..,fiber.g,carb.g,fat.g,protein.g,
         GI_glucose,GL_glucose,added_sugars,sugar.cal,thirty_ins) %>% 
  mutate(delta_weight_6m = (data_bl$weight_gcrc - weight_gcrc)) %>% 
  mutate(weight_change_6m = (weight_gcrc-data_bl$weight_gcrc)) %>%
  mutate(delta_weight_6m_2 = ((data_bl$weight_gcrc - weight_gcrc)/data_bl$weight_gcrc)*100) %>% 
  mutate(qdelta_weight_6m = quantcut(delta_weight_6m,q=5)) %>% 
  mutate(bl_bmi = data_bl$bmi) %>% 
  mutate(bl_carb_cals = data_bl$bl_carb_cals) %>% 
  mutate(carb_cals_6m = carb.g*4) %>% 
  mutate(delta_carbcals_6m = carb_cals_6m-bl_carb_cals) %>% 
  mutate(z_delta_carbcals_6m = (mean(delta_carbcals_6m,na.rm=T)-delta_carbcals_6m)/sd(delta_carbcals_6m,na.rm=T)) %>% 
  mutate(qdelta_carbcals_6m = quantcut(delta_carbcals_6m,q=5)) %>%
  mutate(qdelta_carbcals_6m = quantcut(delta_carbcals_6m,q=5)) %>%
  mutate(addedsugarcals = added_sugars*4) %>%
  mutate(bl_prot_cals = data_bl$bl_prot_cals) %>% 
  mutate(prot_cals_6m = protein.g*4) %>% 
  mutate(delta_protcals_6m = bl_prot_cals - prot_cals_6m) %>% 
  mutate(qdelta_protcals_6m = quantcut(delta_protcals_6m,q=5)) %>% 
  mutate(bl_fat_cals = data_bl$bl_fat_cals) %>% 
  mutate(fat_cals_6m = fat.g*4) %>% 
  mutate(delta_fatcals_6m = fat_cals_6m-bl_fat_cals) %>% 
  mutate(z_delta_fatcals_6m = (mean(delta_fatcals_6m,na.rm=T)-delta_fatcals_6m)/sd(delta_fatcals_6m,na.rm=T)) %>% 
  mutate(qdelta_fatcals_6m = quantcut(delta_fatcals_6m,q=5)) %>% 
  mutate(bl_GI = data_bl$GI_glucose) %>% 
  mutate(delta_GI_6m = bl_GI - GI_glucose) %>% 
  mutate(qdelta_GI_6m = quantcut(delta_GI_6m,q=5)) %>% 
  mutate(bl_GL = data_bl$GL_glucose) %>% 
  mutate(delta_GL_6m = GL_glucose - bl_GL) %>%
  mutate(z_delta_GL_6m = (mean(delta_GL_6m,na.rm=T)-delta_GL_6m)/sd(delta_GL_6m,na.rm=T)) %>% 
  mutate(qdelta_GL_6m  = quantcut(delta_GL_6m,q=5)) %>%  
  mutate(deltaGLq5 = delta_GL_6m <= -99.6) %>% 
  mutate(q3delta_GL_6m  = quantcut(delta_GL_6m,q=3)) %>%  
  mutate(deltaGLq3 = delta_GL_6m <= -80.7) %>% 
  mutate(q4delta_GL_6m  = quantcut(delta_GL_6m,q=4)) %>%  
  mutate(deltaGLq4 = delta_GL_6m <= -89.9) %>% 
  mutate(bl_cal = data_bl$cal) %>% 
  mutate(delta_cal_6m = cal - bl_cal) %>%
  mutate(z_delta_cals_6m = (mean(delta_cal_6m,na.rm=T)-delta_cal_6m)/sd(delta_cal_6m,na.rm=T)) %>%
  mutate(qdelta_cal_6m  = quantcut(delta_cal_6m,q=5)) %>% 
  mutate(tgtohdl_bl = data_bl$bl_tgtohdl) %>% 
  mutate(tgtohdl_6m = lipid_trig_v2/lipid_hdl_v2) %>% 
  mutate(delta_tgtohdl_6m = tgtohdl_6m-tgtohdl_bl) %>% 
  mutate(z_delta_tgtohdl_6m= (mean(delta_tgtohdl_6m,na.rm=T)-delta_tgtohdl_6m)/sd(delta_tgtohdl_6m,na.rm=T)) %>% 
  mutate(qdelta_tgtohdl = quantcut(delta_tgtohdl_6m,q=5)) %>% 
  mutate(bl_30ins = data_bl$thirty_ins) %>%
  mutate(qbl_30ins = quantcut(bl_30ins,q=5)) %>%
  mutate(bl30q5 = bl_30ins >= 127) %>% 
  mutate(q4bl_30ins = quantcut(bl_30ins,q=4)) %>%
  mutate(bl30q4 = bl_30ins >= 118) %>% 
  mutate(q3bl_30ins = quantcut(bl_30ins,q=3)) %>%
  mutate(bl30q3 = bl_30ins >= 103) %>% 
  mutate(bl_satfat = data_bl$saturated_fat..) %>% 
  mutate(delta_satfat = bl_satfat - saturated_fat..) %>%
  mutate(qdelta_satfat = quantcut(delta_satfat,q=5)) %>% 
  mutate(bl_added_sugars = data_bl$added_sugars) %>% 
  mutate(delta_addsugars = bl_added_sugars - added_sugars) %>%
  mutate(qdelta_addsugars = quantcut(delta_addsugars,q=5)) %>% 
  mutate(bl_sugarcal = data_bl$sugar.cal) %>% 
  mutate(delta_sugarcal = bl_sugarcal - sugar.cal) %>%
  mutate(bl_ldl = data_bl$lipid_ldl_v2) %>% 
  mutate(delta_ldl = bl_ldl - lipid_ldl_v2) %>%
  mutate(qdelta_ldl = quantcut(delta_ldl,q=5)) %>%
  mutate(bl_ldlplushdl = data_bl$bl_ldlplushdl) %>% 
  mutate(ldlplushdl_6m = lipid_ldl_v2 + lipid_hdl_v2) %>%
  mutate(delta_ldlplushdl = ldlplushdl_6m - bl_ldlplushdl) %>% 
  mutate(z_delta_ldlplushdl= (mean(delta_ldlplushdl,na.rm=T)-delta_ldlplushdl)/sd(delta_ldlplushdl,na.rm=T)) %>%
  mutate(qdelta_ldlplushdl_6m = quantcut(delta_ldlplushdl,q=5)) %>%
  mutate(q5insGLr=deltaGLq5==T & bl30q5==T) %>%
  mutate(q4insGLr=deltaGLq4==T & bl30q4==T) %>%
  mutate(q3insGLr=deltaGLq3==T & bl30q4==T) %>%
  mutate(across(where(is.numeric), round, 2)) %>%
  distinct()

data_12m <- data %>% filter(redcap_event_name == "12_months_arm_1") %>% 
  select(study_id,baseline_ins,gluc,diet,scr_gender,dxa_percentfat,
         weight_gcrc,lipid_trig_v2,lipid_ldl_v2,timepoint,
         lipid_hdl_v2,MetSyn,cal,fat..,carb..,protein..,bmi,
         saturated_fat..,fiber.g,carb.g,fat.g,protein.g,
         GI_glucose,GL_glucose,added_sugars,sugar.cal,thirty_ins) %>% 
  mutate(delta_weight_12m = (data_bl$weight_gcrc - weight_gcrc)) %>% 
  mutate(weight_change_12m = (weight_gcrc-data_bl$weight_gcrc)) %>%
  mutate(delta_weight_12m_2 = ((data_bl$weight_gcrc - weight_gcrc)/data_bl$weight_gcrc)*100) %>% 
  mutate(qdelta_weight_12m = quantcut(delta_weight_12m,q=5)) %>% 
  mutate(bl_bmi = data_bl$bmi) %>% 
  mutate(bl_carb_cals = data_bl$bl_carb_cals) %>% 
  mutate(carb_cals_12m = carb.g*4) %>% 
  mutate(delta_carbcals_12m = carb_cals_12m-bl_carb_cals) %>% 
  mutate(z_delta_carbcals_12m = (mean(delta_carbcals_12m,na.rm=T)-delta_carbcals_12m)/sd(delta_carbcals_12m,na.rm=T)) %>% 
  mutate(qdelta_carbcals_12m = quantcut(delta_carbcals_12m,q=5)) %>% 
  mutate(addedsugarcals = added_sugars*4) %>%
  mutate(bl_prot_cals = data_bl$bl_prot_cals) %>% 
  mutate(prot_cals_12m = protein.g*4) %>% 
  mutate(delta_protcals_12m = bl_prot_cals - prot_cals_12m) %>% 
  mutate(qdelta_protcals_12m = quantcut(delta_protcals_12m,q=5)) %>% 
  mutate(bl_fat_cals = data_bl$bl_fat_cals) %>% 
  mutate(fat_cals_12m = fat.g*4) %>% 
  mutate(delta_fatcals_12m = fat_cals_12m-bl_fat_cals) %>% 
  mutate(z_delta_fatcals_12m = (mean(delta_fatcals_12m,na.rm=T)-delta_fatcals_12m)/sd(delta_fatcals_12m,na.rm=T)) %>% 
  mutate(qdelta_fatcals_12m = quantcut(delta_fatcals_12m,q=5)) %>% 
  mutate(bl_GI = data_bl$GI_glucose) %>% 
  mutate(delta_GI_12m = bl_GI - GI_glucose) %>% 
  mutate(qdelta_GI_12m = quantcut(delta_GI_12m,q=5)) %>% 
  mutate(bl_GL = data_bl$GL_glucose) %>% 
  mutate(delta_GL_12m = GL_glucose - bl_GL) %>%
  mutate(z_delta_GL_12m = (mean(delta_GL_12m,na.rm=T)-delta_GL_12m)/sd(delta_GL_12m,na.rm=T)) %>% 
  mutate(qdelta_GL_12m  = quantcut(delta_GL_12m,q=5)) %>%  
  mutate(deltaGLq5 = delta_GL_12m <= -89.4) %>% 
  mutate(q3delta_GL_12m  = quantcut(delta_GL_12m,q=3)) %>%  
  mutate(deltaGLq3 = delta_GL_12m <= -66.4) %>% 
  mutate(q4delta_GL_12m  = quantcut(delta_GL_12m,q=4)) %>%  
  mutate(deltaGLq4 = delta_GL_12m <= -81.7) %>% 
  mutate(bl_cal = data_bl$cal) %>% 
  mutate(delta_cal_12m = cal - bl_cal) %>%
  mutate(z_delta_cals_12m = (mean(delta_cal_12m,na.rm=T)-delta_cal_12m)/sd(delta_cal_12m,na.rm=T)) %>%
  mutate(qdelta_cal_12m  = quantcut(delta_cal_12m,q=5)) %>% 
  mutate(tgtohdl_bl = data_bl$bl_tgtohdl) %>% 
  mutate(tgtohdl_12m = lipid_trig_v2/lipid_hdl_v2) %>% 
  mutate(delta_tgtohdl_12m = tgtohdl_12m-tgtohdl_bl) %>% 
  mutate(z_delta_tgtohdl_12m= (mean(delta_tgtohdl_12m,na.rm=T)-delta_tgtohdl_12m)/sd(delta_tgtohdl_12m,na.rm=T)) %>% 
  mutate(qdelta_tgtohdl = quantcut(delta_tgtohdl_12m,q=5)) %>% 
  mutate(bl_30ins = data_bl$thirty_ins) %>%
  mutate(qbl_30ins = quantcut(bl_30ins,q=5)) %>%
  mutate(bl30q5 = bl_30ins >= 127) %>% 
  mutate(q4bl_30ins = quantcut(bl_30ins,q=4)) %>%
  mutate(bl30q4 = bl_30ins >= 118) %>% 
  mutate(q3bl_30ins = quantcut(bl_30ins,q=3)) %>%
  mutate(bl30q3 = bl_30ins >= 103) %>% 
  mutate(bl_satfat = data_bl$saturated_fat..) %>% 
  mutate(delta_satfat = bl_satfat - saturated_fat..) %>%
  mutate(qdelta_satfat = quantcut(delta_satfat,q=5)) %>% 
  mutate(bl_added_sugars = data_bl$added_sugars) %>% 
  mutate(delta_addsugars = bl_added_sugars - added_sugars) %>%
  mutate(qdelta_addsugars = quantcut(delta_addsugars,q=5)) %>% 
  mutate(bl_sugarcal = data_bl$sugar.cal) %>% 
  mutate(delta_sugarcal = bl_sugarcal - sugar.cal) %>%
  mutate(bl_ldl = data_bl$lipid_ldl_v2) %>% 
  mutate(delta_ldl = bl_ldl - lipid_ldl_v2) %>%
  mutate(qdelta_ldl = quantcut(delta_ldl,q=5)) %>%
  mutate(bl_ldlplushdl = data_bl$bl_ldlplushdl) %>% 
  mutate(ldlplushdl_12m = lipid_ldl_v2 + lipid_hdl_v2) %>%
  mutate(delta_ldlplushdl = ldlplushdl_12m - bl_ldlplushdl) %>% 
  mutate(z_delta_ldlplushdl= (mean(delta_ldlplushdl,na.rm=T)-delta_ldlplushdl)/sd(delta_ldlplushdl,na.rm=T)) %>%
  mutate(qdelta_ldlplushdl_12m = quantcut(delta_ldlplushdl,q=5)) %>%
  mutate(q5insGLr=deltaGLq5==T & bl30q5==T) %>%
  mutate(q4insGLr=deltaGLq4==T & bl30q4==T) %>%
  mutate(q3insGLr=deltaGLq3==T & bl30q4==T) %>%
  mutate(across(where(is.numeric), round, 2)) %>%
  distinct() 

##### Filtering by diet groups ####
data_lfd_3m <-data_3m %>% filter (diet == "LFD")
data_lcd_3m <-data_3m %>% filter (diet == "LCD")
data_lfd_6m <-data_6m %>% filter (diet == "LFD")
data_lcd_6m <-data_6m %>% filter (diet == "LCD")
data_lfd_12m <-data_12m %>% filter (diet == "LFD")
data_lcd_12m <-data_12m %>% filter (diet == "LCD")

lcd <- data %>% filter(diet=="LCD")
lfd <- data %>% filter(diet=="LFD")
############################################################################
############################################################################
#####################      DATA IMPUTATION       ###########################
############################################################################
############################################################################
##### Preparing Data #####
# To examine effect modification, we first created a Boolean variable (q5insGLr).
# Then, to impute data, Boolean and character variables need to defined as factors.
data_3m$diet <- as.factor(data_12m$diet)
data_3m$timepoint <- as.factor(data_3m$timepoint)
data_3m$q5insGLr <- as.factor(data_12m$q5insGLr)
data_6m$diet <- as.factor(data_6m$diet)
data_6m$timepoint <- as.factor(data_6m$timepoint)
data_6m$q5insGLr <- as.factor(data_6m$q5insGLr)
data_12m$diet <- as.factor(data_12m$diet)
data_12m$timepoint <- as.factor(data_12m$timepoint)
data_12m$q5insGLr <- as.factor(data_12m$q5insGLr)


# Selecting informative variables for informing the algorithm.
data_3m_imp <- data_3m %>% select(diet,cal,carb..,bmi,scr_gender,q5insGLr,
                                  fat..,fiber.g,weight_change_3m,baseline_ins,
                                  delta_GL_3m,delta_weight_3m,bl_bmi,bl_GI,timepoint)

data_6m_imp <- data_6m %>% select(diet,cal,carb..,bmi,scr_gender,q5insGLr,
                                  fat..,fiber.g,weight_change_6m,baseline_ins,
                                  delta_GL_6m,delta_weight_6m,bl_bmi,bl_GI,timepoint)

data_12m_imp <- data_12m %>% select(diet,cal,carb..,bmi,scr_gender,q5insGLr,
                                    fat..,fiber.g,weight_change_12m,baseline_ins,
                                    delta_GL_12m,delta_weight_12m,bl_bmi,bl_GI,timepoint)

##### Data imputation with random forests #####
set.seed(4321)
imp3m <- missForest(data_3m_imp)
set.seed(4321)
imp6m <- missForest(data_6m_imp)
set.seed(4321)
imp12m <- missForest(data_12m_imp)

data_3m_rf <- as.data.frame(imp3m$ximp)
data_6m_rf <- as.data.frame(imp6m$ximp)
data_12m_rf <- as.data.frame(imp12m$ximp)

##### Data imputation MICE (Multiple Imputation with Chained Equations) #####
data_3m_mice_m <- mice(data_3m_imp,m=5,seed = 4321)
data_6m_mice_m <- mice(data_6m_imp,m=5, seed = 4321)
data_12m_mice_m <- mice(data_12m_imp,m=5, seed = 4321)

data_3m_mice <- complete(data_3m_mice_m,1)
data_6m_mice <- complete(data_6m_mice_m,1)
data_12m_mice <- complete(data_12m_mice_m,1)


##### Weight loss between groups after data imputation #####
# Evaluating weight differences between diet groups with and w/o Data Imputation.
# This calculation produces p-values. Mean weigth loss is obtained below.
kruskal.test(g=data_3m$diet,x=data_3m$weight_change_3m)
kruskal.test(g=data_6m$diet,x=data_6m$weight_change_6m)
kruskal.test(g=data_12m$diet,x=data_12m$weight_change_12m)
kruskal.test(g=data_3m_rf$diet,x=data_3m_rf$weight_change_3m)
kruskal.test(g=data_6m_rf$diet,x=data_6m_rf$weight_change_6m)
kruskal.test(g=data_12m_rf$diet,x=data_12m_rf$weight_change_12m)
kruskal.test(g=data_3m_mice$diet,x=data_3m_mice$weight_change_3m)
kruskal.test(g=data_6m_mice$diet,x=data_6m_mice$weight_change_6m)
kruskal.test(g=data_12m_mice$diet,x=data_12m_mice$weight_change_12m)

# Evaluating effect modification.
emm5q_12m_rf <-  lm(weight_change_12m~q5insGLr,data=data_12m_rf)
summary(emm5q_12m_rf)
emm5q_12m_mice <-  lm(weight_change_12m~q5insGLr,data=data_12m_mice)
summary(emm5q_12m_mice)


############################################################################
##### FIGURE 1. Weight Differences with and without imputation #####
##### Panel A. Without imputation ######
data_or_bl <- data_bl %>% select(diet,weight_change_bl,timepoint)
colnames(data_or_bl) <- c("Diet","WeightChange","Timepoint")
data_or_3m <- data_3m %>% select(diet,weight_change_3m,timepoint)
colnames(data_or_3m) <- c("Diet","WeightChange","Timepoint")
data_or_6m <- data_6m %>% select(diet,weight_change_6m,timepoint)
colnames(data_or_6m) <- c("Diet","WeightChange","Timepoint")
data_or_12m <- data_12m %>% select(diet,weight_change_12m,timepoint)
colnames(data_or_12m) <- c("Diet","WeightChange","Timepoint")
data_or <- bind_rows(data_or_bl,data_or_3m,data_or_6m,data_or_12m)


dfiga <- summarySE(data_or, measurevar="WeightChange", groupvars=c("Timepoint","Diet"),na.rm = T)
# Mean weight loss can be obtained here.
view(dfiga)

dfiga$Timepoint <- factor(dfiga$Timepoint, levels=c("Baseline","3m", "6m", "12m"))
pd <- position_dodge(0.1) 
pfiga <- ggplot(dfiga, aes(x=Timepoint, y=WeightChange, colour=Diet, group=Diet)) + 
  geom_errorbar(aes(ymin=WeightChange-se, ymax=WeightChange+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd,aes(linetype = Diet), colour="black",size=2) +
  geom_point(position=pd, colour="black", size=4)
pfiga <- pfiga+ggtitle("A) Without Imputation")+
  xlab("Time (months)")+ylab("Weight change (kg)")+labs(colour="Diet")+theme_classic()
pfiga <- pfiga+ggplot2::theme(plot.title = element_text(size=22,face = "bold", margin = margin(b=1, unit = "cm")),
                              plot.title.position = "plot",
                              plot.subtitle = element_text(size = 20, face="bold.italic" ,margin = margin(b=1,unit = "cm")),
                              plot.margin = margin(t=1,b=1,l=1,r=1,unit = "cm"),
                              axis.title.x = element_text(size=20, colour = "black",face="bold",margin = margin(t=1,unit = "cm")),
                              axis.title.y = element_text(size=20, colour = "black",face="bold",margin = margin(r=1,unit = "cm")),
                              axis.text = element_text(size=20, colour = "black", face="bold"),
                              legend.box.margin = margin(l=1,unit = "cm"),
                              legend.title = element_text(size=20, colour = "black",face = "bold",margin = margin(b=0.3,unit = "cm")),
                              legend.text = element_text(size=20, colour = "black")
)
pfiga
ggsave("Fig1_A.Weight change without data imputation.bmp", units="cm", width=30, height=20, dpi=1000)

##### Panel B. Random Forest Imputation ######
data_3m_rf <- data_3m_rf %>% select(diet,weight_change_3m,timepoint)
colnames(data_3m_rf) <- c("Diet","WeightChange","Timepoint")
data_6m_rf <- data_6m_rf %>% select(diet,weight_change_6m,timepoint)
colnames(data_6m_rf) <- c("Diet","WeightChange","Timepoint")
data_12m_rf <- data_12m_rf %>% select(diet,weight_change_12m,timepoint)
colnames(data_12m_rf) <- c("Diet","WeightChange","Timepoint")
data_rf <- bind_rows(data_or_bl,data_3m_rf,data_6m_rf,data_12m_rf)


dfigb <- summarySE(data_rf, measurevar="WeightChange", groupvars=c("Timepoint","Diet"),na.rm = T)
# Mean weight loss can be obtained here.
view(dfigb)

dfigb$Timepoint <- factor(dfigb$Timepoint, levels=c("Baseline","3m", "6m", "12m"))
pd <- position_dodge(0.1) 
pfigb <- ggplot(dfigb, aes(x=Timepoint, y=WeightChange, colour=Diet, group=Diet)) + 
  geom_errorbar(aes(ymin=WeightChange-se, ymax=WeightChange+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd,aes(linetype = Diet), colour="black",size=2) +
  geom_point(position=pd, colour="black", size=4)
pfigb <- pfigb+ggtitle("B) With Imputation Using Random Forest Method")+
  xlab("Time (months)")+ylab("Weight change (kg)")+labs(colour="Diet")+theme_classic()
pfigb <- pfigb+ggplot2::theme(plot.title = element_text(size=22,face = "bold", margin = margin(b=1, unit = "cm")),
                              plot.title.position = "plot",
                              plot.subtitle = element_text(size = 20, face="bold.italic" ,margin = margin(b=1,unit = "cm")),
                              plot.margin = margin(t=1,b=1,l=1,r=1,unit = "cm"),
                              axis.title.x = element_text(size=20, colour = "black",face="bold",margin = margin(t=1,unit = "cm")),
                              axis.title.y = element_text(size=20, colour = "black",face="bold",margin = margin(r=1,unit = "cm")),
                              axis.text = element_text(size=20, colour = "black", face="bold"),
                              legend.box.margin = margin(l=1,unit = "cm"),
                              legend.title = element_text(size=20, colour = "black",face = "bold",margin = margin(b=0.3,unit = "cm")),
                              legend.text = element_text(size=20, colour = "black")
)
pfigb
ggsave("Fig1_B.Weight change with random forest.bmp", units="cm", width=30, height=20, dpi=1000)

##### Panel C. MICE Imputation ######
data_3m_mice <- data_3m_mice %>% select(diet,weight_change_3m,timepoint)
colnames(data_3m_mice) <- c("Diet","WeightChange","Timepoint")
data_6m_mice <- data_6m_mice %>% select(diet,weight_change_6m,timepoint)
colnames(data_6m_mice) <- c("Diet","WeightChange","Timepoint")
data_12m_mice <- data_12m_mice %>% select(diet,weight_change_12m,timepoint)
colnames(data_12m_mice) <- c("Diet","WeightChange","Timepoint")
data_mice <- bind_rows(data_or_bl,data_3m_mice,data_6m_mice,data_12m_mice)


dfigc <- summarySE(data_mice, measurevar="WeightChange", groupvars=c("Timepoint","Diet"),na.rm = T)
# Mean weight loss can be obtained here.
view(dfigc)

dfigc$Timepoint <- factor(dfigc$Timepoint, levels=c("Baseline","3m", "6m", "12m"))
pd <- position_dodge(0.1) 
pfigc <- ggplot(dfigc, aes(x=Timepoint, y=WeightChange, colour=Diet, group=Diet)) + 
  geom_errorbar(aes(ymin=WeightChange-se, ymax=WeightChange+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd,aes(linetype = Diet), colour="black",size=2) +
  geom_point(position=pd, colour="black", size=4)
pfigc <- pfigc+ggtitle("C) With Imputation Using MICE Method")+
  xlab("Time (months)")+ylab("Weight change (kg)")+labs(colour="Diet")+theme_classic()
pfigc <- pfigc+ggplot2::theme(plot.title = element_text(size=22,face = "bold", margin = margin(b=1, unit = "cm")),
                              plot.title.position = "plot",
                              plot.subtitle = element_text(size = 20, face="bold.italic" ,margin = margin(b=1,unit = "cm")),
                              plot.margin = margin(t=1,b=1,l=1,r=1,unit = "cm"),
                              axis.title.x = element_text(size=20, colour = "black",face="bold",margin = margin(t=1,unit = "cm")),
                              axis.title.y = element_text(size=20, colour = "black",face="bold",margin = margin(r=1,unit = "cm")),
                              axis.text = element_text(size=20, colour = "black", face="bold"),
                              legend.box.margin = margin(l=1,unit = "cm"),
                              legend.title = element_text(size=20, colour = "black",face = "bold",margin = margin(b=0.3,unit = "cm")),
                              legend.text = element_text(size=20, colour = "black")
)
pfigc
ggsave("Fig1_C.Weight change with MICE.bmp", units="cm", width=30, height=20, dpi=1000)


############################################################################
############################################################################
###############   MAIN RESULTS (WITHOUT IMPUTATION)   ######################
############################################################################
############################################################################
##### FIGURE 2. Dietary Mediators (3-months) #####
##### Weight change adjusted models BL to 3m #####
w_3m_adj_totcal <- lm(delta_weight_3m~cal+scr_gender+bl_bmi+bl_cal+diet,
                      data=data_3m)

w_3m_adj_carbcal <- lm(delta_weight_3m~carb_cals_3m+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_3m)

w_3m_adj_protcal <- lm(delta_weight_3m~prot_cals_3m+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_3m)

w_3m_adj_fatcal <- lm(delta_weight_3m~fat_cals_3m+scr_gender+bl_bmi+bl_cal
                      +diet,data=data_3m)

w_3m_adj_sugarcal <- lm(delta_weight_3m~addedsugarcals+scr_gender+bl_bmi+bl_cal
                        +diet,data=data_3m)

w_3m_adj_fiber <- lm(delta_weight_3m~fiber.g+scr_gender+bl_bmi+bl_cal
                     +diet,data=data_3m)

w_3m_adj_satfat <- lm(delta_weight_3m~saturated_fat..+scr_gender+bl_bmi+bl_cal
                      +diet,data=data_3m)

w_3m_adj_gi <- lm(delta_weight_3m~GI_glucose+scr_gender+bl_bmi+bl_cal
                  +diet,data=data_3m)

w_3m_cal_adj_comp <- compare_performance(w_3m_adj_totcal,w_3m_adj_carbcal,
                                         w_3m_adj_protcal,w_3m_adj_fatcal,
                                         w_3m_adj_sugarcal,w_3m_adj_fiber,
                                         w_3m_adj_satfat,w_3m_adj_gi)
plot(w_3m_cal_adj_comp)
ggsave("Fig2_Adjusted_models_3m.Baseline to 3m.tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

############################################################################
##### FIGURE 3. Dietary Mediators #####
##### Panel A. Glycemic load reduction vs weight loss #####
gl_3m_trend <- ggscatterstats(
  data  = data_3m,
  x     = z_delta_GL_3m,
  y     = delta_weight_3m,
  xlab  = "Glycemic load reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "A) Glycemic Load Reduction",
  subtitle = "Beta = 1.13, SE = 1.19, p = 3.6 x"~10^-9,
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,
                          na.rm = TRUE)
  ) + 
  ggplot2::scale_y_continuous(
  limits = c(-20,20),
  breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
                 ) + 
  coord_cartesian(expand = T)
gl_3m_trend
ggsave("Fig3_A.Weight reduction vs Glycemic Load reduction.tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

fig_3a_beta <- lm(delta_weight_3m~z_delta_GL_3m,data=data_3m)
summary(fig_3a_beta)
length(fig_3a_beta$residuals)

##### Panel B. Fat intake reduction vs weight loss #####
fatcal_3m_trend <- ggscatterstats(
  data  = data_3m,
  x     = z_delta_fatcals_3m,
  y     = delta_weight_3m,
  xlab  = "Fat intake reduction, calories (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "B) Fat Intake Reduction",
  subtitle = "Beta = 0.33, SE = 0.19, p = 0.08",
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE)
  ) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
fatcal_3m_trend
ggsave("Fig3_B.Weight reduction vs fat intake reduction.tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

fig_3b_beta <- lm(delta_weight_3m~z_delta_fatcals_3m,data=data_3m)
summary(fig_3b_beta)
length(fig_3b_beta$residuals)

############################################################################
##### FIGURE 4. Adherence measures #####
##### Panel A. Carbohydrate intake reduction vs weight loss in the LCD #####
lcd_adherence <- ggscatterstats(
  data  = data_lcd_3m,
  x     = z_delta_carbcals_3m,
  y     = delta_weight_3m,
  xlab  = "Carbohydrate intake reduction, calories (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "A) Carbohydrate Reduction (LCD Group)",
  subtitle = "Beta = 2.08, SE = 0.33, p = 2.0 x"~10^-9,
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE) 
  ) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
lcd_adherence
ggsave("Fig4_A.Weight reduction vs carbohydrate intake reduction in the LCD group.tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

fig_4a_beta <- lm(delta_weight_3m~z_delta_carbcals_3m,
                  data=data_lcd_3m)
summary(fig_4a_beta)
length(fig_4a_beta$residuals)

##### Panel B. Fat intake reduction vs weight loss in the LFD #####
lfd_adherence <- ggscatterstats(
  data  = data_lfd_3m,
  x     = z_delta_fatcals_3m,
  y     = delta_weight_3m,
  xlab  = "Fat intake reduction, calories (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "B) Fat Reduction (LFD Group)",
  subtitle = "Beta = 1.09, SE = 0.25, p = 1.3 x"~10^-5,
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE) 
  ) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
lfd_adherence
ggsave("Fig4_B.Weight reduction vs fat intake reduction in the LFD group.tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

fig_4b_beta <- lm(delta_weight_3m~z_delta_fatcals_3m,
                  data=data_lfd_3m)
summary(fig_4b_beta)
length(fig_4b_beta$residuals)

############################################################################
##### FIGURE 5. Biomeasures #####
##### Panel A. Reduction of Tg/HDL vs weight loss #####
tgtohdl_3m_trend <- ggscatterstats(
  data  = data_3m,
  x     = z_delta_tgtohdl_3m,
  y     = delta_weight_3m,
  xlab  = "TG/HDLc reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "A) TG/HDLc Reduction",
  subtitle = "Beta = 1.10, SE = 0.18, p = 3.5 x"~10^-9,
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE) 
  ) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
tgtohdl_3m_trend
ggsave("Fig5_A.Weight reduction vs TG to HDL reduction.tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

fig_5a_beta <- lm(delta_weight_3m~z_delta_tgtohdl_3m,data=data_3m)
summary(fig_5a_beta)
length(fig_5a_beta$residuals)

##### Panel B. Reduction of LDL+HDL vs weight loss #####
ldlplushdl_3m_trend <- ggscatterstats(
  data  = data_3m,
  x     = z_delta_ldlplushdl,
  y     = delta_weight_3m,
  xlab  = "LDLc + HDLc reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "B) LDLc + HDLc Reduction",
  subtitle = "Beta = 0.30, SE = 1.19, p = 0.12",
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE) 
  ) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
ldlplushdl_3m_trend
ggsave("Fig5_B.Weight reduction vs LDL + HDL reduction.tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

fig_5b_beta <- lm(delta_weight_3m~z_delta_ldlplushdl,data=data_3m)
summary(fig_5b_beta)
length(fig_5b_beta$residuals)

############################################################################
##### FIGURE 6. Effect Modification #####
##### Weight loss vs insulin at 30 min vs GL reduction #####
w_3m_df_5q_gl <- aggregate(x=data_3m$delta_weight_3m,
                           by= list(data_3m$qdelta_GL,data_3m$qbl_30ins),
                           FUN=mean,na.rm=TRUE)
colnames(w_3m_df_5q_gl) <- c("GL","qbl30min","kg")

w_3m_hmap_plot_5q <- ggplot(data=w_3m_df_5q_gl,aes(x=GL,y=qbl30min))+
  geom_tile(aes(fill=kg))+
  ggtitle("Weight loss vs Glycemic Load Reduction and Baseline OGTT 30' insulin quintile",
          subtitle = "3 month time point")+
  scale_y_discrete(name="Baseline OGTT insulin 30' quintile",breaks=waiver())+
  scale_x_discrete(name="Change in glycemic goad quintile",breaks=waiver())+
  scale_fill_gradient(name="kg lost",low="#FC6B64",high="#6B77F8",n.breaks=5)+
  geom_text(aes(GL,qbl30min,label=round(kg,1)),color="black",size=6)+
  ggplot2::theme(plot.title = element_text(size=22,face = "bold"),
                 plot.title.position = "plot",
                 plot.subtitle = element_text(size = 18, face="bold.italic" ,margin = margin(b=1,unit = "cm")),
                 plot.margin = margin(t=1,b=1,l=1,r=1,unit = "cm"),
                 axis.title.x = element_text(size=18, colour = "black",face="bold",margin = margin(t=1,unit = "cm")),
                 axis.title.y = element_text(size=18, colour = "black",face="bold",margin = margin(r=1,unit = "cm")),
                 axis.text = element_text(size=14, colour = "black"),
                 legend.box.margin = margin(l=1,unit = "cm"),
                 legend.title = element_text(size=18, colour = "black",face = "bold",margin = margin(b=0.3,unit = "cm")),
                 legend.text = element_text(size=14, colour = "black")
  )
w_3m_hmap_plot_5q
ggsave("Fig6.Effect Modification 3-months.tiff", units="cm", width=35, height=35, dpi=1000, compression = 'lzw')

##### Weight loss vs insulin at 30 min vs GL reduction #####
# To examine effect modification, we transformed "q5insGLr" back to Boolean.
data_3m <- data_3m %>% mutate(q5insGLr=deltaGLq5==T & bl30q5==T)

# Effect Modification using quintiles.
emm5q_3m <-  lm(delta_weight_3m~q5insGLr,data=data_3m)
summary(emm5q_3m)
sum(data_3m$q5insGLr==T, na.rm = TRUE)
length(emm5q_3m$residuals)

# Effect Modification using quartiles.
emm4q_3m <-  lm(delta_weight_3m~q4insGLr,data=data_3m)
summary(emm4q_3m)
sum(data_3m$q4insGLr==T, na.rm = TRUE)
length(emm5q_3m$residuals)

############################################################################
##### TABLE 2 (Table 1 comes from Gardner et al. DIETFITS paper) ############
##### SEM mediation analyses #####
mediation <- data_3m %>% select(delta_weight_3m,delta_cal_3m,delta_GL_3m,diet) %>% 
  na.omit()

m0 <- lm(delta_cal_3m~delta_GL_3m+diet,
         data=mediation)
summary(m0)

m1 <- lm(delta_weight_3m~delta_cal_3m+diet,
         data=mediation)
summary(m1)

m2 <- lm(delta_weight_3m~delta_GL_3m+delta_cal_3m+diet,
         data=mediation)
summary(m2)

set.seed(4321)
med_model <- mediate(m0,m2,treat="diet",
                     mediator="delta_GL_3m",
                     robustSE=T,
                     sims=100,)
summary(med_model)

############################################################################
############################################################################
######################   SUPPLEMENTAL RESULTS  #############################
############################################################################
############################################################################
##### SUPPLEMENTAL TABLE 1 (6-months)  #####
##### SEM mediation analyses #####
mediation_st1 <- data_6m %>% select(delta_weight_6m,delta_cal_6m,delta_GL_6m,diet) %>% 
  na.omit()

m0_st1 <- lm(delta_cal_6m~delta_GL_6m+diet,
         data=mediation_st1)
summary(m0_st1)

m1_st1 <- lm(delta_weight_6m~delta_cal_6m+diet,
         data=mediation_st1)
summary(m1_st1)

m2_st1 <- lm(delta_weight_6m~delta_GL_6m+delta_cal_6m+diet,
         data=mediation_st1)
summary(m2_st1)

set.seed(4321)
med_model_st1 <- mediate(m0_st1,m2_st1,treat="diet",
                     mediator="delta_GL_6m",
                     robustSE=T,
                     sims=100)
summary(med_model_st1)

############################################################################
##### SUPPLEMENTAL FIGURE 1. Dietary Mediators (6-months) #####
##### Weight change adjusted models BL to 6m #####
w_6m_adj_totcal <- lm(delta_weight_6m~cal+scr_gender+bl_bmi+bl_cal+diet,
                      data=data_6m)

w_6m_adj_carbcal <- lm(delta_weight_6m~carb_cals_6m+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_6m)

w_6m_adj_protcal <- lm(delta_weight_6m~prot_cals_6m+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_6m)

w_6m_adj_fatcal <- lm(delta_weight_6m~fat_cals_6m+scr_gender+bl_bmi+bl_cal
                      +diet,data=data_6m)

w_6m_adj_sugarcal <- lm(delta_weight_6m~addedsugarcals+scr_gender+bl_bmi+bl_cal
                        +diet,data=data_6m)

w_6m_adj_fiber <- lm(delta_weight_6m~fiber.g+scr_gender+bl_bmi+bl_cal
                     +diet,data=data_6m)

w_6m_adj_satfat <- lm(delta_weight_6m~saturated_fat..+scr_gender+bl_bmi+bl_cal
                      +diet,data=data_6m)

w_6m_adj_gi <- lm(delta_weight_6m~GI_glucose+scr_gender+bl_bmi+bl_cal
                  +diet,data=data_6m)

w_6m_cal_adj_comp <- compare_performance(w_6m_adj_totcal,w_6m_adj_carbcal,
                                         w_6m_adj_protcal,w_6m_adj_fatcal,
                                         w_6m_adj_sugarcal,w_6m_adj_fiber,
                                         w_6m_adj_satfat,w_6m_adj_gi)
plot(w_6m_cal_adj_comp)
ggsave("SuppFig1_Adjusted_models_6m.Baseline to 6m.tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

############################################################################
##### SUPPLEMENTAL FIGURE 2. Dietary Mediators (6-months)  #####
##### Panel A. Glycemic load reduction vs weight loss #####
gl_6m_trend <- ggscatterstats(
  data  = data_6m,
  x     = z_delta_GL_6m,
  y     = delta_weight_6m,
  xlab  = "Glycemic load reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "A) Glycemic Load Reduction",
  subtitle = "Beta = 1.00, SE = 0.28, p = 0.0003",
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE) 
  ) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
gl_6m_trend
ggsave("SuppFig2_A.Weight reduction vs Glycemic Load reduction.tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

fig_2a_beta_6m <- lm(delta_weight_6m~z_delta_GL_6m,data=data_6m)
summary(fig_2a_beta_6m)
length(fig_2a_beta_6m$residuals)


##### Panel B. Fat intake reduction vs weight loss #####
fatcal_6m_trend <- ggscatterstats(
  data  = data_6m,
  x     = z_delta_fatcals_6m,
  y     = delta_weight_6m,
  xlab  = "Fat intake reduction, calories (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "B) Fat Intake Reduction",
  subtitle = "Beta = 0.16, SE = 0.28, p = 0.57",
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE) 
  ) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
fatcal_6m_trend
ggsave("SuppFig2_B.Weight reduction vs fat intake reduction.tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

fig_2b_beta_6m <- lm(delta_weight_6m~z_delta_fatcals_6m,data=data_6m)
summary(fig_2b_beta_6m)
length(fig_2b_beta_6m$residuals)

############################################################################
##### SUPPLEMENTAL FIGURE 3. Adherence measures (6-months)  #####
##### Panel A. Carbohydrate intake reduction vs weight loss in the LCD #####
lcd_adherence_6m <- ggscatterstats(
  data  = data_lcd_6m,
  x     = z_delta_carbcals_6m,
  y     = delta_weight_6m,
  xlab  = "Carbohydrate intake reduction, calories (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "A) Fat Reduction (LCD Group)",
  subtitle = "Beta = 2.35, SE = 0.45, p = 5.5 x"~10^-7,
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE) 
  ) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
lcd_adherence_6m
ggsave("SuppFig3_A.Weight reduction vs carbohydrate intake reduction in the LCD group.tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

fig_3a_beta_6m <- lm(delta_weight_6m~z_delta_carbcals_6m,
                  data=data_lcd_6m)
summary(fig_3a_beta_6m)
length(fig_3a_beta_6m$residuals)

##### Panel B. Fat intake reduction vs weight loss in the LFD #####
lfd_adherence_6m <- ggscatterstats(
  data  = data_lfd_6m,
  x     = z_delta_fatcals_6m,
  y     = delta_weight_6m,
  xlab  = "Fat intake reduction, calories (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "B) Fat Reduction (LFD Group)",
  subtitle = "Beta = 1.20, SE = 0.40, p = 0.003",
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE) 
  ) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
lfd_adherence_6m
ggsave("SuppFig3_B.Weight reduction vs fat intake reduction in the LFD group.tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

fig_3b_beta_6m <- lm(delta_weight_6m~z_delta_fatcals_6m,
                  data=data_lfd_6m)
summary(fig_3b_beta_6m)
length(fig_3b_beta_6m$residuals)

############################################################################
##### SUPPLEMENTAL FIGURE 4. Biomeasures (6-months) #####
##### Panel A. Reduction of Tg/HDL vs weight loss #####
tgtohdl_6m_trend <- ggscatterstats(
  data  = data_6m,
  x     = z_delta_tgtohdl_6m,
  y     = delta_weight_6m,
  xlab  = "TG/HDLc reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "A) TG/HDLc Reduction",
  subtitle = "Beta = 1.66, SE = 0.27, p = 1.1 x"~10^-9,
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE) 
  ) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
tgtohdl_6m_trend
ggsave("SuppFig4_A.Weight reduction vs TG to HDL reduction.tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

fig_4a_beta_6m <- lm(delta_weight_6m~z_delta_tgtohdl_6m,data=data_6m)
summary(fig_4a_beta_6m)
length(fig_4a_beta_6m$residuals)

##### Panel B. Reduction of LDL+HDL vs weight loss #####
ldlplushdl_6m_trend <- ggscatterstats(
  data  = data_6m,
  x     = z_delta_ldlplushdl,
  y     = delta_weight_6m,
  xlab  = "LDLc + HDLc reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "B) LDLc + HDLc Reduction",
  subtitle = "Beta = -0.15, SE = 0.28, p = 0.58",
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE) 
  ) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
ldlplushdl_6m_trend
ggsave("SuppFig4_B.Weight reduction vs LDL + HDL reduction.tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

fig_4b_beta_6m <- lm(delta_weight_6m~z_delta_ldlplushdl,data=data_6m)
summary(fig_4b_beta_6m)
length(fig_4b_beta_6m$residuals)

############################################################################
##### SUPPLEMENTAL FIGURE 5. Effect Modification (6-months) #####
##### Weight loss vs insulin at 30 min vs GL reduction #####
w_6m_df_5q_gl <- aggregate(x=data_6m$delta_weight_6m,
                           by= list(data_6m$qdelta_GL,data_6m$qbl_30ins),
                           FUN=mean,na.rm=TRUE)
colnames(w_6m_df_5q_gl) <- c("GL","qbl30min","kg")

w_6m_hmap_plot_5q <- ggplot(data=w_6m_df_5q_gl,aes(x=GL,y=qbl30min))+
  geom_tile(aes(fill=kg))+
  ggtitle("Weight loss vs Glycemic Load Reduction and Baseline OGTT 30' insulin quintile",
          subtitle = "6 month time point")+
  scale_y_discrete(name="Baseline OGTT insulin 30' quintile",breaks=waiver())+
  scale_x_discrete(name="Change in glycemic goad quintile",breaks=waiver())+
  scale_fill_gradient(name="kg lost",low="#FC6B64",high="#6B77F8",n.breaks=5)+
  geom_text(aes(GL,qbl30min,label=round(kg,1)),color="black",size=6)+
  ggplot2::theme(plot.title = element_text(size=22,face = "bold"),
                 plot.title.position = "plot",
                 plot.subtitle = element_text(size = 18, face="bold.italic" ,margin = margin(b=1,unit = "cm")),
                 plot.margin = margin(t=1,b=1,l=1,r=1,unit = "cm"),
                 axis.title.x = element_text(size=18, colour = "black",face="bold",margin = margin(t=1,unit = "cm")),
                 axis.title.y = element_text(size=18, colour = "black",face="bold",margin = margin(r=1,unit = "cm")),
                 axis.text = element_text(size=14, colour = "black"),
                 legend.box.margin = margin(l=1,unit = "cm"),
                 legend.title = element_text(size=18, colour = "black",face = "bold",margin = margin(b=0.3,unit = "cm")),
                 legend.text = element_text(size=14, colour = "black")
  )
w_6m_hmap_plot_5q
ggsave("SuppFig5.Effect Modification 6-months.tiff", units="cm", width=35, height=35, dpi=1000, compression = 'lzw')
##### Weight loss vs insulin at 30 min vs GL reduction #####
emm5q_6m <-  lm(delta_weight_6m~q5insGLr,data=data_6m)
summary(emm5q_6m)
sum(data_6m$q5insGLr==T, na.rm = TRUE)
length(emm5q_6m$residuals)

emm4q_6m <-  lm(delta_weight_6m~q4insGLr,data=data_6m)
summary(emm4q_6m)
sum(data_6m$q4insGLr==T, na.rm = TRUE)
############################################################################
##### SUPPLEMENTAL TABLE 2 (12-months) #####
##### SEM MEDIATION ANALYSIS #####
mediation_st2 <- data_12m %>% select(delta_weight_12m,delta_cal_12m,delta_GL_12m,diet) %>% 
  na.omit()

m0_st2 <- lm(delta_cal_12m~delta_GL_12m+diet,
             data=mediation_st2)
summary(m0_st2)

m1_st2 <- lm(delta_weight_12m~delta_cal_12m+diet,
             data=mediation_st2)
summary(m1_st2)

m2_st2 <- lm(delta_weight_12m~delta_GL_12m+delta_cal_12m+diet,
             data=mediation_st2)
summary(m2_st2)

set.seed(4321)
med_model_st2 <- mediate(m0_st2,m2_st2,treat="diet",
                         mediator="delta_GL_12m",
                         robustSE=T,
                         sims=100)
summary(med_model_st2)

############################################################################
##### SUPPLEMENTAL FIGURE 6 (12-months) #####
##### Weight change adjusted models BL to 12m #####
w_12m_adj_totcal <- lm(delta_weight_12m~cal+scr_gender+bl_bmi+bl_cal+diet,
                       data=data_12m)

w_12m_adj_carbcal <- lm(delta_weight_12m~carb_cals_12m+scr_gender+bl_bmi+bl_cal
                        +diet,data=data_12m)

w_12m_adj_protcal <- lm(delta_weight_12m~prot_cals_12m+scr_gender+bl_bmi+bl_cal
                        +diet,data=data_12m)

w_12m_adj_fatcal <- lm(delta_weight_12m~fat_cals_12m+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_12m)

w_12m_adj_sugarcal <- lm(delta_weight_12m~addedsugarcals+scr_gender+bl_bmi+bl_cal
                         +diet,data=data_12m)

w_12m_adj_fiber <- lm(delta_weight_12m~fiber.g+scr_gender+bl_bmi+bl_cal
                      +diet,data=data_12m)

w_12m_adj_satfat <- lm(delta_weight_12m~saturated_fat..+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_12m)

w_12m_adj_gi <- lm(delta_weight_12m~GI_glucose+scr_gender+bl_bmi+bl_cal
                   +diet,data=data_12m)

w_12m_cal_adj_comp <- compare_performance(w_12m_adj_totcal,w_12m_adj_carbcal,
                                          w_12m_adj_protcal,w_12m_adj_fatcal,
                                          w_12m_adj_sugarcal,w_12m_adj_fiber,
                                          w_12m_adj_satfat,w_12m_adj_gi)
plot(w_12m_cal_adj_comp)
ggsave("SuppFig6_Adjusted_models_12m.Baseline to 12m.tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

############################################################################
##### SUPPLEMENTAL FIGURE 7. Dietary Mediators (12-months) #####
##### Panel A. Glycemic load reduction vs weight loss #####
gl_12m_trend <- ggscatterstats(
  data  = data_12m,
  x     = z_delta_GL_12m,
  y     = delta_weight_12m,
  xlab  = "Glycemic load reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "A) Glycemic Load Reduction",
  subtitle = "Beta = 0.94, SE = 0.35, p = 0.008",
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE) 
  ) +  
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
gl_12m_trend
ggsave("SupFig7_A.Weight reduction vs Glycemic Load reduction(12m).tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

supfig_7a_beta_12 <- lm(delta_weight_12m~z_delta_GL_12m,
                  data=data_12m)
summary(supfig_7a_beta_12)
length(supfig_7a_beta_12$residuals)

##### Panel B. Fat intake reduction vs weight loss #####
cal_12m_trend <- ggscatterstats(
  data  = data_12m,
  x     = z_delta_fatcals_12m,
  y     = delta_weight_12m,
  xlab  = "Fat intake reduction, calories (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "B) Fat Intake Reduction",
  subtitle = "Beta = 0.05, SE = 0.34, p = 0.89",
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE) 
  ) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
cal_12m_trend
ggsave("SupFig7_B.Weight reduction vs fat intake reduction (12m).tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

supfig_7b_beta <- lm(delta_weight_12m~z_delta_fatcals_12m,
                  data=data_12m)
summary(supfig_7b_beta)
length(supfig_7b_beta$residuals)

############################################################################
##### SUPPLEMENTAL FIGURE 8. Adherence measures (12-months) #####
##### Panel A. Carbohydrate intake reduction vs weight loss in the LCD #####
lcd_adherence_12 <- ggscatterstats(
  data  = data_lcd_12m,
  x     = z_delta_carbcals_12m,
  y     = delta_weight_12m,
  xlab  = "Carbohydrate intake reduction, calories (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "A) Carbohydrate Reduction (LCD Group)",
  subtitle = "Beta = 1.25, SE = 0.55, p = 0.02",
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE) 
  ) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
lcd_adherence_12
ggsave("SupFig8_A.Weight reduction vs carbohydrate intake reduction in the LCD group (12m).tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

supfig_8a_beta <- lm(delta_weight_12m~z_delta_carbcals_12m,
                  data=data_lcd_12m)
summary(supfig_8a_beta)
length(supfig_8a_beta$residuals)

##### Panel B. Fat intake reduction vs weight loss in the LFD #####
lfd_adherence_12 <- ggscatterstats(
  data  = data_lfd_12m,
  x     = z_delta_fatcals_12m,
  y     = delta_weight_12m,
  xlab  = "Fat intake reduction, calories (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "B) Fat Reduction (LFD Group)",
  subtitle = "Beta = 1.34, SE = 0.51, p = 0.009",
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE) 
  ) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
lfd_adherence_12
ggsave("SupFig8_B.Weight reduction vs fat intake reduction in the LFD group (12m).tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

supfig_8b_beta <- lm(delta_weight_12m~z_delta_fatcals_12m,
                  data=data_lfd_12m)
summary(supfig_8b_beta)
length(supfig_8b_beta$residuals)

############################################################################
##### SUPPLEMENTAL FIGURE 9. Biomeasures (12-months) #####
##### Panel A. Reduction of Tg/HDL vs weight loss #####
tgtohdl_12m_trend <- ggscatterstats(
  data  = data_12m,
  x     = z_delta_tgtohdl_12m,
  y     = delta_weight_12m,
  xlab  = "Change in TG/HDLc (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "A) TG/HDLc Reduction",
  subtitle = "Beta = 2.59, SE = 0.31, p = 1.5 x"~10^-15,
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE) 
  ) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
tgtohdl_12m_trend
ggsave("SupFig9_A.Weight reduction vs TG to HDL reduction (12m).tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

supfig_4a_beta <- lm(delta_weight_12m~z_delta_tgtohdl_12m,
                  data=data_12m)
summary(supfig_4a_beta)
length(supfig_4a_beta$residuals)

##### Panel B. Reduction of LDL+HDL vs weight loss #####
ldlplushdl_12m_trend <- ggscatterstats(
  data  = data_12m,
  x     = z_delta_ldlplushdl,
  y     = delta_weight_12m,
  xlab  = "Change in LDLc + HDLc (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "B) LDLc + HDLc Reduction",
  subtitle = "Beta = -0.31, SE = 0.34, p = 0.35",
  results.subtitle = F, marginal = F,
  smooth.line.args = list(size = 1.5, color = "black", method = "lm", formula = y ~ x,                           na.rm = TRUE) 
  ) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=26,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.line = element_line(color = "black",size = 1, linetype = "solid"),
                 axis.title.y = element_text(size=24,colour = "black",face="bold",margin = margin(r=0.5,unit="cm")),
                 axis.title.x = element_text(size=24,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=24,colour = "black",face="bold", margin = margin(t=1,r=1,unit="cm")),
                 plot.subtitle = element_text(size=22,colour = "black",hjust = 0.5,margin = margin(b=1,t=1,unit="cm")),
                 plot.caption = element_text(size=22,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) + 
  coord_cartesian(expand = T)
ldlplushdl_12m_trend
ggsave("SupFig9_B.Weight reduction vs LDL + HDL reduction(12m).tiff", units="cm", width=30, height=30, dpi=1000, compression = 'lzw')

supfig_4b_beta <- lm(delta_weight_12m~z_delta_ldlplushdl,
                  data=data_12m)
summary(supfig_4b_beta)
length(supfig_4b_beta$residuals)

############################################################################
##### SUPPLEMENTAL FIGURE 10. Effect Modification (12-months) #####
##### Weight loss vs insulin at 30 min vs GL reduction #####
w_12m_df_5q_gl <- aggregate(x=data_12m$delta_weight_12m,
                            by= list(data_12m$qdelta_GL,data_12m$qbl_30ins),
                            FUN=mean,na.rm=TRUE)
colnames(w_12m_df_5q_gl) <- c("GL","qbl30min","kg")

w_12m_hmap_plot_5q <- ggplot(data=w_12m_df_5q_gl,aes(x=GL,y=qbl30min))+
  geom_tile(aes(fill=kg))+
  ggtitle("Weight loss vs Glycemic Load Reduction and Baseline OGTT 30' insulin quintile",
          subtitle = "12 month time point")+
  scale_y_discrete(name="Baseline OGTT insulin 30' quintile",breaks=waiver())+
  scale_x_discrete(name="Change in glycemic goad quintile",breaks=waiver())+
  scale_fill_gradient(name="kg lost",low="#FC6B64",high="#6B77F8",n.breaks=5)+
  geom_text(aes(GL,qbl30min,label=round(kg,1)),color="black",size=6)+
  ggplot2::theme(plot.title = element_text(size=22,face = "bold"),
                 plot.title.position = "plot",
                 plot.subtitle = element_text(size = 18, face="bold.italic" ,margin = margin(b=1,unit = "cm")),
                 plot.margin = margin(t=1,b=1,l=1,r=1,unit = "cm"),
                 axis.title.x = element_text(size=18, colour = "black",face="bold",margin = margin(t=1,unit = "cm")),
                 axis.title.y = element_text(size=18, colour = "black",face="bold",margin = margin(r=1,unit = "cm")),
                 axis.text = element_text(size=14, colour = "black"),
                 legend.box.margin = margin(l=1,unit = "cm"),
                 legend.title = element_text(size=18, colour = "black",face = "bold",margin = margin(b=0.3,unit = "cm")),
                 legend.text = element_text(size=14, colour = "black")
  )
w_12m_hmap_plot_5q
ggsave("SupFig10.Effect Modification 12-months.tiff", units="cm", width=35, height=35, dpi=1000, compression = 'lzw')

##### Weight loss vs insulin at 30 min vs GL reduction #####
emm5q12 <-  lm(delta_weight_12m~q5insGLr,data=data_12m)
summary(emm5q12)
sum(data_12m$q5insGLr==T, na.rm = TRUE)
length(emm5q12$residuals)

emm4q_12m <-  lm(delta_weight_12m~q4insGLr,data=data_12m)
summary(emm4q_12m)
sum(data_12m$q4insGLr==T, na.rm = TRUE)
############################################################################
############################################################################
