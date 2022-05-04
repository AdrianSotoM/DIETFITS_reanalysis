#### SET UP ####
#Data available at https://github.com/AdrianSotoM
#Please, address all correspondence about this code to adrian.sotom@incmnsz_mx
#Running it takes a few minutes, I've set up a sound to let you know 
#when it's done in case you want to do something else while it runs.

#Working directory setup
setwd("C:/Users/adria/Dropbox/UIEM/LEAD/Proyectos/DIETFITS")

# Since "easystats" is currently not available in CRAN, if you don't have it,
# you'll need to install manually.
# install.packages("easystats", repos = "https://easystats.r-universe.dev")

# Now, confirm you have "pacman" installed. If you don't have "pacman" but want
# to install it, remove the # in the line below and press "Enter".
# install.packages("pacman") 

#Packages setup
pacman::p_load(dplyr,tidyr,ggstatsplot,readxl,tableone,easystats,dagitty,
               patchwork,MASS,see,qqplotr,bootStepAIC,performance,ggdag,
               rpart,rpart.plot,gtools,broom,lmtest,visdat,report,
               parameters,ggcharts,conflicted,car,rattle,cvms,lavaan,
               mlogit,MLmetrics,beepr,readr,haven,dagitty,mediation)

#Solving duplicate functions conflicts
conflict_prefer("select","dplyr")
conflict_prefer("filter", "dplyr")

#### DATA UPLOAD ####
data <- read.csv(url("https://osf.io/ztysq/download"))

#Participant 2001 will be excluded due to a likely erroneous lipid measurement.
data <- data %>% filter(study_id != "2001")
#############################################################################
##### DATA FILTERING BY DIFFERENT TIMEPOINTS####
data_bl <- data %>%  filter(redcap_event_name == "baseline_arm_1") %>% 
  select(id,baseline_ins,gluc,diet,redcap_event_name,scr_gender,
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
  mutate(across(where(is.numeric), round, 2)) %>%
  distinct()

data_3m <- data %>% filter(redcap_event_name == "3_months_arm_1") %>% 
  select(baseline_ins,gluc,diet,scr_gender,dxa_percentfat,
         weight_gcrc,lipid_trig_v2,lipid_ldl_v2,
         lipid_hdl_v2,MetSyn,cal,fat..,carb..,protein..,bmi,
         saturated_fat..,fiber.g,carb.g,fat.g,protein.g,
         GI_glucose,GL_glucose,added_sugars,sugar.cal,thirty_ins) %>% 
  mutate(delta_weight_3m = data_bl$weight_gcrc - weight_gcrc) %>% 
  mutate(qdelta_weight_3m = quantcut(delta_weight_3m,q=5)) %>% 
  mutate(bl_bmi = data_bl$bmi) %>% 
  mutate(bl_carb_cals = data_bl$bl_carb_cals) %>% 
  mutate(carb_cals_3m = carb.g*4) %>% 
  mutate(delta_carbcals_3m = carb_cals_3m-bl_carb_cals) %>% 
  mutate(z_delta_carbcals_3m = (mean(delta_carbcals_3m,na.rm=T)-delta_carbcals_3m)/sd(delta_carbcals_3m,na.rm=T)) %>% 
  mutate(qdelta_carbcals_3m = quantcut(delta_carbcals_3m,q=5)) %>% 
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
  mutate(deltaGLq5 = delta_GL_3m < -107) %>% 
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
  mutate(bl30q5 = bl_30ins > 127) %>% 
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
  mutate(across(where(is.numeric), round, 2)) %>%
  distinct() 

data_6m <- data %>% filter(redcap_event_name == "6_months_arm_1") %>% 
  select(baseline_ins,gluc,diet,scr_gender,dxa_percentfat,
         weight_gcrc,lipid_trig_v2,lipid_ldl_v2,
         lipid_hdl_v2,MetSyn,cal,fat..,carb..,protein..,bmi,
         saturated_fat..,fiber.g,carb.g,fat.g,protein.g,
         GI_glucose,GL_glucose,added_sugars,sugar.cal,thirty_ins) %>% 
  mutate(delta_weight_6m = data_bl$weight_gcrc - weight_gcrc) %>% 
  mutate(qdelta_weight_6m = quantcut(delta_weight_6m,q=5)) %>% 
  mutate(bl_bmi = data_bl$bmi) %>% 
  mutate(bl_carb_cals = data_bl$bl_carb_cals) %>% 
  mutate(carb_cals_6m = carb.g*4) %>% 
  mutate(delta_carbcals_6m = carb_cals_6m-bl_carb_cals) %>% 
  mutate(z_delta_carbcals_6m = (mean(delta_carbcals_6m,na.rm=T)-delta_carbcals_6m)/sd(delta_carbcals_6m,na.rm=T)) %>% 
  mutate(qdelta_carbcals_6m = quantcut(delta_carbcals_6m,q=5)) %>% 
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
  mutate(deltaGLq5 = delta_GL_6m < -107) %>% 
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
  mutate(bl30q5 = bl_30ins > 127) %>% 
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
  mutate(across(where(is.numeric), round, 2)) %>%
  distinct() 

data_12m <- data %>% filter(redcap_event_name == "12_months_arm_1") %>% 
  select(baseline_ins,gluc,diet,scr_gender,dxa_percentfat,
         weight_gcrc,lipid_trig_v2,lipid_ldl_v2,
         lipid_hdl_v2,MetSyn,cal,fat..,carb..,protein..,bmi,
         saturated_fat..,fiber.g,carb.g,fat.g,protein.g,
         GI_glucose,GL_glucose,added_sugars,sugar.cal,thirty_ins) %>% 
  mutate(delta_weight_12m = data_bl$weight_gcrc - weight_gcrc) %>% 
  mutate(qdelta_weight_12m = quantcut(delta_weight_12m,q=5)) %>% 
  mutate(bl_bmi = data_bl$bmi) %>% 
  mutate(bl_carb_cals = data_bl$bl_carb_cals) %>% 
  mutate(carb_cals_12m = carb.g*4) %>% 
  mutate(delta_carbcals_12m = carb_cals_12m-bl_carb_cals) %>% 
  mutate(z_delta_carbcals_12m = (mean(delta_carbcals_12m,na.rm=T)-delta_carbcals_12m)/sd(delta_carbcals_12m,na.rm=T)) %>% 
  mutate(qdelta_carbcals_12m = quantcut(delta_carbcals_12m,q=5)) %>% 
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
  mutate(deltaGLq5 = delta_GL_12m < -107) %>% 
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
  mutate(bl30q5 = bl_30ins > 127) %>% 
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
  mutate(delta_weight_6m12m = weight_gcrc - data_6m$weight_gcrc) %>%
  mutate(across(where(is.numeric), round, 2)) %>%
  distinct()

##### DATA FILTERING BY DIFFERENT DIET GROUPS####
data_lfd_3m <-data_3m %>% filter (diet == "Purple")
data_lcd_3m <-data_3m %>% filter (diet == "Blue")
data_lfd_6m <-data_6m %>% filter (diet == "Purple")
data_lcd_6m <-data_6m %>% filter (diet == "Blue")
data_lfd_12m <-data_12m %>% filter (diet == "Purple")
data_lcd_12m <-data_12m %>% filter (diet == "Blue")



##### TABLES. WEIGHT LOSS AND LIPID MEASURES TIMEPOINT COMPARISONS #####
##### Weightloss #####
weight_loss <- data.frame(data_3m$delta_weight_3m,
                          data_6m$delta_weight_6m,
                          data_12m$delta_weight_12m,data$diet)
names(weight_loss) <- c('3m', '6m', '12m','diet')

table1 <- CreateTableOne(data=weight_loss,strata = "diet",
                         vars = c("3m","6m","12m"))
summary(table1)

##### Lipid measures #####
biomeasures <- data.frame(data_3m$tgtohdl_bl,
                          data_3m$bl_ldl,
                          data_3m$tgtohdl_3m,
                          data_6m$tgtohdl_6m,
                          data_12m$tgtohdl_12m,
                          data_3m$ldlplushdl,
                          data_6m$ldlplushdl,
                          data_12m$ldlplushdl,
                          data$diet)
names(biomeasures) <- c('bltgtohdl','blldlplushdl','tgtohdl3m', 
                        'tgtohdl6m', 'tgtohdl12m',
                        'ldlplushdl3m','ldlplushdl6m',
                        'ldlplushdl12m','diet')

table2 <- CreateContTable(data=biomeasures,strata = "diet",
                         vars = c('tgtohdl3m', 'tgtohdl6m', 'tgtohdl12m',
                                  'ldlplushdl3m','ldlplushdl6m','ldlplushdl12m'))
summary(table2)

##### Lipid within group comparison #####
lcdbiom <- biomeasures %>% filter(diet=="Blue")
lfdbiom <- biomeasures %>% filter(diet=="Purple")


lcdtghdlbl3 <- wilcox.test(lcdbiom$bltgtohdl,lcdbiom$tgtohdl3m,paired = T)
lcdtghdl36 <- wilcox.test(lcdbiom$bltgtohdl,lcdbiom$tgtohdl6m,paired = T)
lcdtghdl12 <- wilcox.test(lcdbiom$bltgtohdl,lcdbiom$tgtohdl12m,paired = T)
lcdldlplushdlbl3 <- wilcox.test(lcdbiom$blldlplushdl,lcdbiom$ldlplushdl3m,paired = T)
lcdldlplushdl36 <- wilcox.test(lcdbiom$blldlplushdl,lcdbiom$ldlplushdl6m,paired = T)
lcdldlplushdl12 <- wilcox.test(lcdbiom$blldlplushdl,lcdbiom$ldlplushdl12m,paired = T)

lfdtghdlbl3 <- wilcox.test(lfdbiom$bltgtohdl,lfdbiom$tgtohdl3m,paired = T)
lfdtghdl36 <- wilcox.test(lfdbiom$bltgtohdl,lfdbiom$tgtohdl6m,paired = T)
lfdtghdl12 <- wilcox.test(lfdbiom$bltgtohdl,lfdbiom$tgtohdl12m,paired = T)
lfdldlplushdlbl3 <- wilcox.test(lfdbiom$blldlplushdl,lfdbiom$ldlplushdl3m,paired = T)
lfdldlplushdl36 <- wilcox.test(lfdbiom$blldlplushdl,lfdbiom$ldlplushdl6m,paired = T)
lfdldlplushdl12 <- wilcox.test(lfdbiom$blldlplushdl,lfdbiom$ldlplushdl12m,paired = T)

lfdtghdlbl3
lfdtghdl36
lfdtghdl12
lfdldlplushdlbl3
lfdldlplushdl36
lfdldlplushdl12
lcdtghdlbl3
lcdtghdl36
lcdtghdl12
lcdldlplushdlbl3
lcdldlplushdl36
lcdldlplushdl12

############################################################################
##### FIGURE 1 #####
##### Panel A. Weight change adjusted models BL to 3m #####
w_3m_adj_totcal <- lm(delta_weight_3m~cal+scr_gender+bl_bmi+bl_cal+diet,
                      data=data_3m)
summary(w_3m_adj_totcal)
w_3m_adj_totcal_betas <- ggcoefstats(w_3m_adj_totcal)
w_3m_adj_totcal_betas

w_3m_adj_carbcal <- lm(delta_weight_3m~carb_cals_3m+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_3m)
summary(w_3m_adj_carbcal)
w_3m_adj_carbcal_betas <- ggcoefstats(w_3m_adj_carbcal)
w_3m_adj_carbcal_betas

w_3m_adj_protcal <- lm(delta_weight_3m~prot_cals_3m+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_3m)
summary(w_3m_adj_protcal)
w_3m_adj_protcal_betas <- ggcoefstats(w_3m_adj_protcal)
w_3m_adj_protcal_betas

w_3m_adj_fatcal <- lm(delta_weight_3m~fat_cals_3m+scr_gender+bl_bmi+bl_cal
                      +diet,data=data_3m)
summary(w_3m_adj_fatcal)
w_3m_adj_fatcal_betas <- ggcoefstats(w_3m_adj_fatcal)
w_3m_adj_fatcal_betas

w_3m_adj_sugarcal <- lm(delta_weight_3m~sugar.cal+scr_gender+bl_bmi+bl_cal
                        +diet,data=data_3m)
summary(w_3m_adj_sugarcal)
w_3m_adj_sugarcal_betas <- ggcoefstats(w_3m_adj_sugarcal)
w_3m_adj_sugarcal_betas

w_3m_adj_fiber <- lm(delta_weight_3m~fiber.g+scr_gender+bl_bmi+bl_cal
                     +diet,data=data_3m)
summary(w_3m_adj_fiber)
w_3m_adj_fiber_betas <- ggcoefstats(w_3m_adj_fiber)
w_3m_adj_fiber_betas

w_3m_adj_satfat <- lm(delta_weight_3m~saturated_fat..+scr_gender+bl_bmi+bl_cal
                      +diet,data=data_3m)
summary(w_3m_adj_satfat)
w_3m_adj_satfat_betas <- ggcoefstats(w_3m_adj_satfat)
w_3m_adj_satfat_betas

w_3m_cal_adj_comp <- compare_performance(w_3m_adj_totcal,w_3m_adj_carbcal,
                                         w_3m_adj_protcal,w_3m_adj_fatcal,
                                         w_3m_adj_sugarcal,w_3m_adj_fiber,
                                         w_3m_adj_satfat)
w_3m_cal_adj_comp
plot(w_3m_cal_adj_comp)
##### Panel B. Weight change adjusted models 6 to 12m #####
w_12m_adj_totcal <- lm(delta_weight_12m~cal+scr_gender+bl_bmi+bl_cal+diet,
                       data=data_12m)
summary(w_12m_adj_totcal)
w_12m_adj_totcal_betas <- ggcoefstats(w_12m_adj_totcal)
w_12m_adj_totcal_betas

w_12m_adj_carbcal <- lm(delta_weight_12m~carb_cals_12m+scr_gender+bl_bmi+bl_cal
                        +diet,data=data_12m)
summary(w_12m_adj_carbcal)
w_12m_adj_carbcal_betas <- ggcoefstats(w_12m_adj_carbcal)
w_12m_adj_carbcal_betas

w_12m_adj_protcal <- lm(delta_weight_12m~prot_cals_12m+scr_gender+bl_bmi+bl_cal
                        +diet,data=data_12m)
summary(w_12m_adj_protcal)
w_12m_adj_protcal_betas <- ggcoefstats(w_12m_adj_protcal)
w_12m_adj_protcal_betas

w_12m_adj_fatcal <- lm(delta_weight_12m~fat_cals_12m+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_12m)
summary(w_12m_adj_fatcal)
w_12m_adj_fatcal_betas <- ggcoefstats(w_12m_adj_fatcal)
w_12m_adj_fatcal_betas

w_12m_adj_sugarcal <- lm(delta_weight_12m~sugar.cal+scr_gender+bl_bmi+bl_cal
                         +diet,data=data_12m)
summary(w_12m_adj_sugarcal)
w_12m_adj_sugarcal_betas <- ggcoefstats(w_12m_adj_sugarcal)
w_12m_adj_sugarcal_betas

w_12m_adj_fiber <- lm(delta_weight_12m~fiber.g+scr_gender+bl_bmi+bl_cal
                      +diet,data=data_12m)
summary(w_12m_adj_fiber)
w_12m_adj_fiber_betas <- ggcoefstats(w_12m_adj_fiber)
w_12m_adj_fiber_betas

w_12m_adj_satfat <- lm(delta_weight_12m~saturated_fat..+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_12m)
summary(w_12m_adj_satfat)
w_12m_adj_satfat_betas <- ggcoefstats(w_12m_adj_satfat)
w_12m_adj_satfat_betas

w_12m_cal_adj_comp <- compare_performance(w_12m_adj_totcal,w_12m_adj_carbcal,
                                          w_12m_adj_protcal,w_12m_adj_fatcal,
                                          w_12m_adj_sugarcal,w_12m_adj_fiber,
                                          w_12m_adj_satfat)
w_12m_cal_adj_comp
plot(w_12m_cal_adj_comp)
############################################################################
##### FIGURE 2. Dietary Mediators #####
##### Panel A. Glycemic load reduction vs weight loss #####
gl_3m_trend <- ggscatterstats(
  data  = data_3m,
  x     = z_delta_GL_3m,
  y     = delta_weight_3m,
  xlab  = "GL reduction (z-score)",
  ylab  = "Weight reduction(kg)",
  title = "Weight reduction vs Glycemic Load reduction"
)
gl_3m_trend

fig_2a_beta <- lm(delta_weight_3m~delta_GL_3m+scr_gender+bl_bmi+bl_cal+diet,
                       data=data_3m)
summary(fig_2a_beta)
plot_fig_2a_beta <- ggcoefstats(fig_2a_beta)
plot_fig_2a_beta

##### Panel B. Fat intake reduction vs weight loss #####
cal_3m_trend <- ggscatterstats(
  data  = data_3m,
  x     = z_delta_fatcals_3m,
  y     = delta_weight_3m,
  xlab  = "Fat intake reduction (total calories z-score)",
  ylab  = "Weight reduction (kg)",
  title = "Weight reduction vs fat intake reduction"
)
cal_3m_trend

fig_2b_beta <- lm(delta_weight_3m~delta_cal_3m+scr_gender+bl_bmi+bl_cal+diet,
                  data=data_3m)

summary(fig_2b_beta)
beta_fig_2b_beta <- ggcoefstats(fig_2b_beta)
beta_fig_2b_beta

############################################################################
##### FIGURE 3. Adherence measures #####
##### Panel A. Carbohydrate intake reduction vs weight loss in the LCD #####
lcd_adherence <- ggscatterstats(
  data  = data_lcd_3m,
  x     = z_delta_carbcals_3m,
  y     = delta_weight_3m,
  xlab  = "Carbohydrate intake reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "Weight reduction vs carbohydrate intake reduction in the LCD group"
)
lcd_adherence

fig_3a_beta <- lm(delta_weight_3m~delta_carbcals_3m+scr_gender+bl_bmi+bl_cal,
                       data=data_lcd_3m)

summary(fig_3a_beta)

plot_fig_3a_beta <- ggcoefstats(fig_3a_beta)
plot_fig_3a_beta

##### Panel B. Fat intake reduction vs weight loss in the LFD #####
lfd_adherence <- ggscatterstats(
  data  = data_lfd_3m,
  x     = z_delta_fatcals_3m,
  y     = delta_weight_3m,
  xlab  = "Fat intake reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "Weight reduction vs fat intake reduction in the LFD group"
)
lfd_adherence

fig_3b_beta <- lm(delta_weight_3m~delta_fatcals_3m+scr_gender+bl_bmi+bl_cal,
                  data=data_lfd_3m)

summary(fig_3b_beta)

plot_fig_3b_beta <- ggcoefstats(fig_3b_beta)
plot_fig_3b_beta

############################################################################
##### FIGURE 4. Lipid Biomeasures #####
##### Panel A. Reduction of Tg/HDL vs weight loss #####
tgtohdl_3m_trend <- ggscatterstats(
  data  = data_3m,
  x     = z_delta_tgtohdl_3m,
  y     = delta_weight_3m,
  xlab  = "TG/HDL reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "Weight reduction vs TG/HDL reduction"
)
tgtohdl_3m_trend

fig_4a_beta <- lm(delta_weight_3m~delta_tgtohdl_3m+scr_gender+bl_bmi+bl_cal+diet,
                  data=data_3m)

summary(fig_4a_beta)

plot_fig_4a_beta <- ggcoefstats(fig_4a_beta)
plot_fig_4a_beta

##### Panel B. Reduction of LDL+HDL vs weight loss #####
ldlplushdl_3m_trend <- ggscatterstats(
  data  = data_3m,
  x     = z_delta_ldlplushdl,
  y     = delta_weight_3m,
  xlab  = "LDL + HDL reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "Weight reduction vs LDL + HDL reduction"
)
ldlplushdl_3m_trend

fig_4b_beta <- lm(delta_weight_3m~delta_ldlplushdl+scr_gender+bl_bmi+bl_cal+diet,
                  data=data_3m)

summary(fig_4b_beta)

plot_fig_4b_beta <- ggcoefstats(fig_4b_beta)
plot_fig_4b_beta

############################################################################
##### FIGURE 5. Effect Modification #####
##### Panel A. Weight loss vs insulin at 30 min vs GL reduction #####
w_3m_df_5q_gl <- aggregate(x=data_3m$delta_weight_3m,
                           by= list(data_3m$qdelta_GL,data_3m$qbl_30ins),
                           FUN=mean,na.rm=TRUE)
colnames(w_3m_df_5q_gl) <- c("GL","qbl30min","weightloss")

w_3m_hmap_plot_5q <- ggplot(data=w_3m_df_5q_gl,aes(x=GL,y=qbl30min))+
  geom_tile(aes(fill=weightloss))+
  ggtitle("Weightloss by delta GL and qbl30min LFD")+
  scale_y_discrete(name="Baseline insulin 30 quintile",breaks=waiver())+
  scale_x_discrete(name="Change in Glycemic Load quintile",breaks=waiver())+
  scale_fill_gradient(name="weightloss",low="#FC6B64",high="#6B77F8",n.breaks=10)+
  geom_text(aes(GL,qbl30min,label=round(weightloss,1)),color="black",size=4)+
  ggplot2::theme(text = element_text(size = 12),
                 plot.title = element_text(size=20,face = "bold"),
                 axis.title.x = element_text(size=16, colour = "black"),
                 axis.title.y = element_text(size=16, colour = "black")
  )
w_3m_hmap_plot_5q

##### Panel B. Weight loss vs insulin at 30 min vs GL reduction #####
em <- lm(delta_weight_3m~bl30q5*deltaGLq5+diet+bl_bmi+bl_cal,
         data=data_3m)
summary(em)
em_betas <- ggcoefstats(em)
em_betas

############################################################################
##### SUPPLEMENTAL MATERIAL 
##### SUPPLEMENTAL TABLE 1 #####
##### SEM MEDIATION ANALYSIS #####
mediation <- data_3m %>% select(delta_weight_3m,delta_cal_3m,delta_GL_3m,diet) %>% 
  na.omit()

m1 <- lm(delta_cal_3m~delta_GL_3m+diet,
                    data=mediation)
summary(m1)

m2 <- lm(delta_weight_3m~delta_GL_3m+delta_cal_3m+diet,
                 data=mediation)
summary(m2)

m3 <- lm(delta_weight_3m~delta_cal_3m+diet,
         data=mediation)
summary(m3)

med_model <- mediate(m1,m2,treat="diet",
                     mediator="delta_GL_3m",
                     robustSE=T,
                     sims=100)
summary(med_model)


#############################################################################
##### SUPPLEMENTAL FIGURE 1 #####
##### DIETFITS' DAG #####
dag <- dagify (
  weightloss ~ glr,
  weightloss ~ enired,
  enired ~ lowcarb,
  enired ~ lowfat,
  glr ~ nutedu,
  glr ~ lowcarb,
  nutedu ~ lowcarb,
  nutedu ~ lowfat,
  lowcarb ~ enrollment,
  lowfat ~ enrollment,
  outcome = "weightloss",
  exposure = "enrollment",
  labels = c(
              weightloss = "Weightloss",
              enrollment = "DIETFITS enrollment",
              nutedu = "Nutritional education",
              enired="Energy intake reduction",
              glr="Glycemic load reduction",
              lowcarb="Carbohydrate intake reduction",
              lowfat="Fat intake reduction"
             )
                )
set.seed(555)
ggdag_status(dag, use_labels = "label",text=F) + theme_dag()

##### DIETFITS' DAG OPEN PATHS#####
set.seed(555)
ggdag_paths(dag,shadow = T,use_labels = "label",text=F) + theme_dag()

