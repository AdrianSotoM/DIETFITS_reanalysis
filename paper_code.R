#### SET UP ####
# Data is openly available at DIETFITS' Open Science Framework repository.

# Actually, this code downloads it automatically, you don't need to save the
# data file in your computer.

# Please, address all correspondence about this code to adrian.sotom@tec.mx
# Also please, note that some graphs were edited in the paper to improve their
#readability.

# If you are only interested in corroborating our results and figures, just 
# click "Ctrl + Shift + Enter" and go get a coffee. I programmed a funny sound 
# to let you it finished running the whole code.

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
               rpart,rpart.plot,gtools,broom,lmtest,visdat,report,beepr,
               parameters,ggcharts,conflicted,car,rattle,cvms,lavaan,
               mlogit,MLmetrics,beepr,readr,haven,dagitty,mediation)

#Solving duplicate functions conflicts
conflict_prefer("select","dplyr")
conflict_prefer("filter", "dplyr")

#### DATA UPLOAD ####
data <- read.csv(url("https://osf.io/ztysq/download"))
#Participant 2001 will be excluded due to a likely erroneous lipid measurement.
data <- data %>% filter(study_id != "2001")
############################################################################
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
  mutate(delta_weight_3m = (data_bl$weight_gcrc - weight_gcrc)) %>% 
  mutate(delta_weight_3m_2 = ((data_bl$weight_gcrc - weight_gcrc)/data_bl$weight_gcrc)*100) %>% 
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
  mutate(delta_weight_6m = (data_bl$weight_gcrc - weight_gcrc)) %>% 
  mutate(delta_weight_6m_2 = ((data_bl$weight_gcrc - weight_gcrc)/data_bl$weight_gcrc)*100) %>%
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
  mutate(delta_weight_12m = (data_bl$weight_gcrc - weight_gcrc)) %>% 
  mutate(delta_weight_12m_2 = ((data_bl$weight_gcrc - weight_gcrc)/data_bl$weight_gcrc)*100) %>%
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



############################################################################
##### FIGURE 1 #####
##### Weight change adjusted models BL to 3m #####
w_3m_adj_totcal <- lm(delta_weight_3m~cal+scr_gender+bl_bmi+bl_cal+diet,
                      data=data_3m)

w_3m_adj_carbcal <- lm(delta_weight_3m~carb_cals_3m+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_3m)

w_3m_adj_protcal <- lm(delta_weight_3m~prot_cals_3m+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_3m)

w_3m_adj_fatcal <- lm(delta_weight_3m~fat_cals_3m+scr_gender+bl_bmi+bl_cal
                      +diet,data=data_3m)

w_3m_adj_sugarcal <- lm(delta_weight_3m~sugar.cal+scr_gender+bl_bmi+bl_cal
                        +diet,data=data_3m)

w_3m_adj_fiber <- lm(delta_weight_3m~fiber.g+scr_gender+bl_bmi+bl_cal
                     +diet,data=data_3m)

w_3m_adj_satfat <- lm(delta_weight_3m~saturated_fat..+scr_gender+bl_bmi+bl_cal
                      +diet,data=data_3m)

w_3m_cal_adj_comp <- compare_performance(w_3m_adj_totcal,w_3m_adj_carbcal,
                                         w_3m_adj_protcal,w_3m_adj_fatcal,
                                         w_3m_adj_sugarcal,w_3m_adj_fiber,
                                         w_3m_adj_satfat)
w_3m_cal_adj_comp
plot(w_3m_cal_adj_comp)
ggsave("Fig1_Adjusted models. Baseline to 3m.tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

############################################################################
##### FIGURE 2. Dietary Mediators #####
##### Panel A. Glycemic load reduction vs weight loss #####
gl_3m_trend <- ggscatterstats(
  data  = data_3m,
  x     = z_delta_GL_3m,
  y     = delta_weight_3m,
  xlab  = "Glycemic load reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "A) Weight Reduction vs Glycemic Load Reduction",
  subtitle = "ß = 6.0, p = 3.64 x 10^9",
  caption = "3 month timepoint",
  results.subtitle = F,
  marginal = F
  ) + 
  ggplot2::scale_y_continuous(
  limits = c(-20,20),
  breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
                 ) +
coord_cartesian(expand = F)
gl_3m_trend
ggsave("Fig2_A.Weight reduction vs Glycemic Load reduction.tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

fig_2a_beta <- lm(delta_weight_3m~z_delta_GL_3m,data=data_3m)
summary(fig_2a_beta)
plot_fig_2a_beta <- ggcoefstats(fig_2a_beta)
plot_fig_2a_beta

##### Panel B. Fat intake reduction vs weight loss #####
fatcal_3m_trend <- ggscatterstats(
  data  = data_3m,
  x     = z_delta_fatcals_3m,
  y     = delta_weight_3m,
  xlab  = "Calories from fat intake reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "B) Weight Reduction vs Fat Intake Reduction",
  subtitle = "ß = 1.75, p = 0.08",
  caption = "3 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
fatcal_3m_trend
ggsave("Fig2_B.Weight reduction vs fat intake reduction.tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

fig_2b_beta <- lm(delta_weight_3m~z_delta_fatcals_3m,data=data_3m)
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
  xlab  = "Calories form carbohydrate intake reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "A) Weight Reduction vs Carbohydrate Reduction (LCD Group)",
  subtitle = "ß = 6.22, p = 1.96 x 10^9",
  caption = "12 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
lcd_adherence
ggsave("Fig3_A.Weight reduction vs carbohydrate intake reduction in the LCD group.tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

fig_3a_beta <- lm(delta_weight_3m~z_delta_carbcals_3m,
                  data=data_lcd_3m)

summary(fig_3a_beta)

plot_fig_3a_beta <- ggcoefstats(fig_3a_beta)
plot_fig_3a_beta

##### Panel B. Fat intake reduction vs weight loss in the LFD #####
lfd_adherence <- ggscatterstats(
  data  = data_lfd_3m,
  x     = z_delta_fatcals_3m,
  y     = delta_weight_3m,
  xlab  = "Calories from fat intake reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "B) Weight Reduction vs Fat Reduction (LFD Group)",
  subtitle = "ß = 4.45, p = 1.26 x 10^5",
  caption = "3 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
lfd_adherence
ggsave("Fig3_B.Weight reduction vs fat intake reduction in the LFD group.tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

fig_3b_beta <- lm(delta_weight_3m~z_delta_fatcals_3m,
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
  title = "A) Weight Reduction vs TG/HDL Reduction",
  subtitle = "ß = 6.01, p = 3.51 x 10^9",
  caption = "3 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
tgtohdl_3m_trend
ggsave("Fig4_A.Weight reduction vs TG to HDL reduction.tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

fig_4a_beta <- lm(delta_weight_3m~z_delta_tgtohdl_3m,data=data_3m)

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
  title = "B) Weight Reduction vs LDL+HDL Reduction",
  subtitle = "ß = 1.57, p = 0.12",
  caption = "3 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
ldlplushdl_3m_trend
ggsave("Fig4_B.Weight reduction vs LDL + HDL reduction.tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

fig_4b_beta <- lm(delta_weight_3m~z_delta_ldlplushdl,data=data_3m)

summary(fig_4b_beta)

plot_fig_4b_beta <- ggcoefstats(fig_4b_beta)
plot_fig_4b_beta

############################################################################
##### FIGURE 5. Effect Modification #####
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
ggsave("Fig5.Effect Modification 3-months.tiff", units="cm", width=35, height=35, dpi=600, compression = 'lzw')
##### Weight loss vs insulin at 30 min vs GL reduction #####
em <- lm(delta_weight_3m~bl30q5*deltaGLq5+diet,
         data=data_3m)
summary(em)
em_betas <- ggcoefstats(em)
em_betas
############################################################################
##### TABLE 1 #####
##### SEM MEDIATION ANALYSIS #####
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

med_model <- mediate(m0,m2,treat="diet",
                     mediator="delta_GL_3m",
                     robustSE=T,
                     sims=100)
summary(med_model)

#########################SUPPLEMENTAL MATERIAL##############################
############################################################################
##### SUPPLEMENTAL TABLE 1 (12-months) #####
##### SEM MEDIATION ANALYSIS #####
mediation_12 <- data_12m %>% select(delta_weight_12m,delta_cal_12m,delta_GL_12m,
                                    diet) %>% 
  na.omit()

m0_12 <- lm(delta_cal_12m~delta_GL_12m+diet,
            data=mediation_12)
summary(m0_12)

m1_12 <- lm(delta_weight_12m~delta_cal_12m+diet,
            data=mediation_12)
summary(m1_12)

m2_12 <- lm(delta_weight_12m~delta_GL_12m+delta_cal_12m+diet,
            data=mediation_12)
summary(m2_12)

med_model_12 <- mediate(m0_12,m2_12,treat="diet",
                        mediator="delta_GL_12m",
                        robustSE=T,
                        sims=100)
summary(med_model_12)

############################################################################
##### SUPPLEMENTAL FIGURE 1 #####
##### Weight change adjusted models BL to 12m #####
w_12m_adj_totcal <- lm(delta_weight_12m~cal+scr_gender+bl_bmi+bl_cal+diet,
                       data=data_12m)

w_12m_adj_carbcal <- lm(delta_weight_12m~carb_cals_12m+scr_gender+bl_bmi+bl_cal
                        +diet,data=data_12m)

w_12m_adj_protcal <- lm(delta_weight_12m~prot_cals_12m+scr_gender+bl_bmi+bl_cal
                        +diet,data=data_12m)

w_12m_adj_fatcal <- lm(delta_weight_12m~fat_cals_12m+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_12m)

w_12m_adj_sugarcal <- lm(delta_weight_12m~sugar.cal+scr_gender+bl_bmi+bl_cal
                         +diet,data=data_12m)

w_12m_adj_fiber <- lm(delta_weight_12m~fiber.g+scr_gender+bl_bmi+bl_cal
                      +diet,data=data_12m)

w_12m_adj_satfat <- lm(delta_weight_12m~saturated_fat..+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_12m)

w_12m_cal_adj_comp <- compare_performance(w_12m_adj_totcal,w_12m_adj_carbcal,
                                          w_12m_adj_protcal,w_12m_adj_fatcal,
                                          w_12m_adj_sugarcal,w_12m_adj_fiber,
                                          w_12m_adj_satfat)
w_12m_cal_adj_comp
plot(w_12m_cal_adj_comp)
ggsave("SupFig1_Adjusted models. Baseline to 12m.tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')
############################################################################
##### SUPPLEMENTAL FIGURE 2. Dietary Mediators (12-months) #####
##### Panel A. Glycemic load reduction vs weight loss #####
gl_12m_trend <- ggscatterstats(
  data  = data_12m,
  x     = z_delta_GL_12m,
  y     = delta_weight_12m,
  xlab  = "Glycemic load reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "A) Weight Reduction vs Glycemic Load Reduction",
  subtitle = "ß = -2.65, p = 8.25 x 10^3",
  caption = "12 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
gl_12m_trend
ggsave("SupFig2_A.Weight reduction vs Glycemic Load reduction(12m).tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')


supfig_2a_beta_12 <- lm(delta_weight_12m~delta_GL_12m,
                  data=data_12m)
summary(supfig_2a_beta_12)
plot_supfig_2a_beta_12m <- ggcoefstats(supfig_2a_beta_12)
plot_supfig_2a_beta_12m

##### Panel B. Fat intake reduction vs weight loss #####
cal_12m_trend <- ggscatterstats(
  data  = data_12m,
  x     = z_delta_fatcals_12m,
  y     = delta_weight_12m,
  xlab  = "Calories from fat reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "B) Weight Reduction vs Fat Intake Reduction",
  subtitle = "ß = -0.14, p = 0.89",
  caption = "12 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
cal_12m_trend
ggsave("SupFig2_B.Weight reduction vs fat intake reduction (12m).tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')


supfig_2b_beta <- lm(delta_weight_12m~delta_fatcals_12m,
                  data=data_12m)

summary(supfig_2b_beta)
beta_supfig_2b_beta <- ggcoefstats(supfig_2b_beta)
beta_supfig_2b_beta

############################################################################
##### SUPPLEMENTAL FIGURE 3. Adherence measures (12-months) #####
##### Panel A. Carbohydrate intake reduction vs weight loss in the LCD #####
lcd_adherence_12 <- ggscatterstats(
  data  = data_lcd_12m,
  x     = z_delta_carbcals_12m,
  y     = delta_weight_12m,
  xlab  = "Calories from carbohydrate intake reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "A) Weight Reduction vs Carbohydrate Reduction (LCD Group)",
  subtitle = "ß = -2.27, p = 0.02",
  caption = "12 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
lcd_adherence_12
ggsave("SupFig3_A.Weight reduction vs carbohydrate intake reduction in the LCD group (12m).tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

supfig_3a_beta <- lm(delta_weight_12m~delta_carbcals_12m,
                  data=data_lcd_12m)

summary(supfig_3a_beta)

plot_supfig_3a_beta <- ggcoefstats(supfig_3a_beta)
plot_supfig_3a_beta

##### Panel B. Fat intake reduction vs weight loss in the LFD #####
lfd_adherence_12 <- ggscatterstats(
  data  = data_lfd_12m,
  x     = z_delta_fatcals_12m,
  y     = delta_weight_12m,
  xlab  = "Calories from fat intake reduction (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "B) Weight Reduction vs Fat Reduction (LFD Group)",
  subtitle = "ß = -2.65, p = 8.06 x 10^3",
  caption = "12 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
lfd_adherence_12
ggsave("SupFig3_B.Weight reduction vs fat intake reduction in the LFD group (12m).tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')


supfig_3b_beta <- lm(delta_weight_12m~delta_fatcals_12m,
                  data=data_lfd_12m)

summary(supfig_3b_beta)

plot_supfig_3b_beta <- ggcoefstats(supfig_3b_beta)
plot_supfig_3b_beta

############################################################################
##### SUPPLEMENTAL FIGURE 4. Lipid Biomeasures (12-months) #####
##### Panel A. Reduction of Tg/HDL vs weight loss #####
tgtohdl_12m_trend <- ggscatterstats(
  data  = data_12m,
  x     = z_delta_tgtohdl_12m,
  y     = delta_weight_12m,
  xlab  = "Change in TG/HDL (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "A) Weight Reduction vs TG/HDL Reduction",
  subtitle = "ß = -8.29, p = 1.44 x 10^15",
  caption = "12 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
tgtohdl_12m_trend
ggsave("SupFig4_A.Weight reduction vs TG to HDL reduction (12m).tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')


supfig_4a_beta <- lm(delta_weight_12m~delta_tgtohdl_12m,
                  data=data_12m)

summary(supfig_4a_beta)

plot_supfig_4a_beta <- ggcoefstats(supfig_4a_beta)
plot_supfig_4a_beta

##### Panel B. Reduction of LDL+HDL vs weight loss #####
ldlplushdl_12m_trend <- ggscatterstats(
  data  = data_12m,
  x     = z_delta_ldlplushdl,
  y     = delta_weight_12m,
  xlab  = "Change in LDL + HDL (z-score)",
  ylab  = "Weight reduction (kg)",
  title = "B) Weight Reduction vs LDL+HDL Reduction",
  subtitle = "ß = 0.93, p = 0.35",
  caption = "12 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
ldlplushdl_12m_trend
ggsave("SupFig4_B.Weight reduction vs LDL + HDL reduction(12m).tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

supfig_4b_beta <- lm(delta_weight_12m~delta_ldlplushdl,
                  data=data_12m)

summary(supfig_4b_beta)

plot_supfig_4b_beta <- ggcoefstats(supfig_4b_beta)
plot_supfig_4b_beta

############################################################################
##### SUPPLEMENTAL FIGURE 5. Effect Modification (12-months) #####
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
ggsave("SupFig5.Effect Modification 12-months.tiff", units="cm", width=35, height=35, dpi=600, compression = 'lzw')

##### Weight loss vs insulin at 30 min vs GL reduction #####
em12 <- lm(delta_weight_12m~bl30q5*deltaGLq5+diet+bl_bmi+bl_cal,
         data=data_12m)
summary(em12)
em_betas12 <- ggcoefstats(em12)
em_betas12

############################################################################
############################################################################
######### WEIGHT LOSS AS THE % DIFFERENCE FROM BASELINE ###################
############################################################################
############################################################################
##### FIGURE 1 Weight loss (Def2) #####
##### Weight change adjusted models BL to 3m #####
w_3m_adj_totcal_2 <- lm(delta_weight_3m_2~cal+scr_gender+bl_bmi+bl_cal+diet,
                      data=data_3m)

w_3m_adj_carbcal_2 <- lm(delta_weight_3m_2~carb_cals_3m+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_3m)

w_3m_adj_protcal_2 <- lm(delta_weight_3m_2~prot_cals_3m+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_3m)

w_3m_adj_fatcal_2 <- lm(delta_weight_3m_2~fat_cals_3m+scr_gender+bl_bmi+bl_cal
                      +diet,data=data_3m)

w_3m_adj_sugarcal_2 <- lm(delta_weight_3m_2~sugar.cal+scr_gender+bl_bmi+bl_cal
                        +diet,data=data_3m)

w_3m_adj_fiber_2 <- lm(delta_weight_3m_2~fiber.g+scr_gender+bl_bmi+bl_cal
                     +diet,data=data_3m)

w_3m_adj_satfat_2 <- lm(delta_weight_3m_2~saturated_fat..+scr_gender+bl_bmi+bl_cal
                      +diet,data=data_3m)

w_3m_cal_adj_comp_2 <- compare_performance(w_3m_adj_totcal_2,w_3m_adj_carbcal_2,
                                         w_3m_adj_protcal_2,w_3m_adj_fatcal_2,
                                         w_3m_adj_sugarcal_2,w_3m_adj_fiber_2,
                                         w_3m_adj_satfat_2)
w_3m_cal_adj_comp_2
plot(w_3m_cal_adj_comp_2)
ggsave("Def2_Fig1_Adjusted models. Baseline to 3m.tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')
############################################################################
##### FIGURE 2. Dietary Mediators (Def2) #####
##### Panel A. Glycemic load reduction vs weight loss #####
gl_3m_trend_2 <- ggscatterstats(
  data  = data_3m,
  x     = z_delta_GL_3m,
  y     = delta_weight_3m_2,
  xlab  = "GL reduction (z-score)",
  ylab  = "Weight reduction(%)",
  title = "Weight reduction vs Glycemic Load reduction",
  subtitle = "ß = 6.08, p = 2.35 x 10^9",
  caption = "3 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
gl_3m_trend_2
ggsave("Def2_Fig2_A.Weight reduction vs Glycemic Load reduction.tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')


d2fig_2a_beta_2 <- lm(delta_weight_3m_2~z_delta_GL_3m,
                  data=data_3m)
summary(d2fig_2a_beta_2)
d2plot_fig_2a_beta_2 <- ggcoefstats(d2fig_2a_beta_2)
d2plot_fig_2a_beta_2

##### Panel B. Fat intake reduction vs weight loss #####
fatcal_3m_trend_2 <- ggscatterstats(
  data  = data_3m,
  x     = z_delta_fatcals_3m,
  y     = delta_weight_3m_2,
  xlab  = "Fat intake reduction (total calories z-score)",
  ylab  = "Weight reduction (%)",
  title = "Weight reduction vs fat intake reduction",
  subtitle = "ß = 0.54, p = 0.59",
  caption = "3 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
fatcal_3m_trend_2
ggsave("Def2_Fig2_B.Weight reduction vs fat intake reduction.tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

d2fig_2b_beta_2 <- lm(delta_weight_3m_2~z_delta_fatcals_3m,
                  data=data_3m)

summary(d2fig_2b_beta_2)
d2beta_fig_2b_beta_2 <- ggcoefstats(d2fig_2b_beta_2)
d2beta_fig_2b_beta_2

############################################################################
##### FIGURE 3. Adherence measures (Def2) #####
##### Panel A. Carbohydrate intake reduction vs weight loss in the LCD #####
lcd_adherence_2 <- ggscatterstats(
  data  = data_lcd_3m,
  x     = z_delta_carbcals_3m,
  y     = delta_weight_3m_2,
  xlab  = "Carbohydrate intake reduction (z-score)",
  ylab  = "Weight reduction (%)",
  title = "Weight reduction vs carbohydrate intake reduction in the LCD group",
  subtitle = "ß = 5.79, p = 2 x 10^8",
  caption = "3 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
lcd_adherence_2
ggsave("Def2_Fig3_A.Weight reduction vs carbohydrate intake reduction in the LCD group.tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')


d2fig_3a_beta_2 <- lm(delta_weight_3m_2~z_delta_carbcals_3m,
                  data=data_lcd_3m)

summary(d2fig_3a_beta_2)

d2plot_fig_3a_beta_2 <- ggcoefstats(d2fig_3a_beta_2)
d2plot_fig_3a_beta_2

##### Panel B. Fat intake reduction vs weight loss in the LFD #####
lfd_adherence_2 <- ggscatterstats(
  data  = data_lfd_3m,
  x     = z_delta_fatcals_3m,
  y     = delta_weight_3m_2,
  xlab  = "Fat intake reduction (z-score)",
  ylab  = "Weight reduction (%)",
  title = "Weight reduction vs fat intake reduction in the LFD group",
  subtitle = "ß = 3.46, p = 6.35 x 10^4",
  caption = "3 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
lfd_adherence_2
ggsave("Def2_Fig3_B.Weight reduction vs fat intake reduction in the LFD group.tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

d2fig_3b_beta_2 <- lm(delta_weight_3m_2~z_delta_fatcals_3m,
                  data=data_lfd_3m)

summary(d2fig_3b_beta_2)

d2plot_fig_3b_beta_2 <- ggcoefstats(d2fig_3b_beta_2)
d2plot_fig_3b_beta_2

############################################################################
##### FIGURE 4. Lipid Biomeasures (Def2) #####
##### Panel A. Reduction of Tg/HDL vs weight loss #####
tgtohdl_3m_trend_2 <- ggscatterstats(
  data  = data_3m,
  x     = z_delta_tgtohdl_3m,
  y     = delta_weight_3m_2,
  xlab  = "TG/HDL reduction (z-score)",
  ylab  = "Weight reduction (%)",
  title = "Weight reduction vs TG to HDL reduction",
  subtitle = "ß = 6.13, p = 1.77 x 10^9",
  caption = "3 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
tgtohdl_3m_trend_2
ggsave("Def2_Fig4_A.Weight reduction vs TG to HDL reduction.tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

d2fig_4a_beta_2 <- lm(delta_weight_3m_2~z_delta_tgtohdl_3m,
                  data=data_3m)

summary(d2fig_4a_beta_2)

d2plot_fig_4a_beta_2 <- ggcoefstats(d2fig_4a_beta_2)
d2plot_fig_4a_beta_2

##### Panel B. Reduction of LDL+HDL vs weight loss #####
ldlplushdl_3m_trend_2 <- ggscatterstats(
  data  = data_3m,
  x     = z_delta_ldlplushdl,
  y     = delta_weight_3m_2,
  xlab  = "LDL + HDL reduction (z-score)",
  ylab  = "Weight reduction (%)",
  title = "Weight reduction vs LDL + HDL reduction",
  subtitle = "ß = 1.26, p = 0.21",
  caption = "3 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
ldlplushdl_3m_trend_2
ggsave("Def2_Fig4_B.Weight reduction vs LDL + HDL reduction.tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')


d2fig_4b_beta_2 <- lm(delta_weight_3m_2~z_delta_ldlplushdl,
                  data=data_3m)

summary(d2fig_4b_beta_2)

d2plot_fig_4b_beta_2 <- ggcoefstats(d2fig_4b_beta_2)
d2plot_fig_4b_beta_2

############################################################################
##### FIGURE 5. Effect Modification (Def2) #####
##### Weight loss vs insulin at 30 min vs GL reduction #####
w_3m_df_5q_gl_2 <- aggregate(x=data_3m$delta_weight_3m_2,
                           by= list(data_3m$qdelta_GL,data_3m$qbl_30ins),
                           FUN=mean,na.rm=TRUE)
colnames(w_3m_df_5q_gl_2) <- c("GL","qbl30min","weightloss")

w_3m_hmap_plot_5q_2 <- ggplot(data=w_3m_df_5q_gl_2,aes(x=GL,y=qbl30min))+
  geom_tile(aes(fill=weightloss))+
  ggtitle("Weight loss vs Glycemic Load Reduction and Baseline OGTT 30' insulin quintile",
          subtitle = "3 month time point")+
  scale_y_discrete(name="Baseline OGTT insulin 30' quintile",breaks=waiver())+
  scale_x_discrete(name="Change in glycemic goad quintile",breaks=waiver())+
  scale_fill_gradient(name="% change",low="#FC6B64",high="#6B77F8",n.breaks=5)+
  geom_text(aes(GL,qbl30min,label=round(weightloss,1)),color="black",size=6)+
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
w_3m_hmap_plot_5q_2
ggsave("Def2_Fig5.Effect Modification 3-months.tiff", units="cm", width=35, height=35, dpi=600, compression = 'lzw')

##### Weight loss vs insulin at 30 min vs GL reduction #####
em_2 <- lm(delta_weight_3m_2~bl30q5*deltaGLq5+diet,
         data=data_3m)
summary(em_2)
em_betas_2 <- ggcoefstats(em_2)
em_betas_2

############################################################################
##### TABLE 1 (Def2) #####
##### SEM MEDIATION ANALYSIS (Def2) #####
mediation_2 <- data_3m %>% select(delta_weight_3m_2,delta_cal_3m,delta_GL_3m,diet) %>% 
  na.omit()

m0_2 <- lm(delta_cal_3m~delta_GL_3m+diet,
         data=mediation_2)
summary(m0_2)

m1_2 <- lm(delta_weight_3m_2~delta_cal_3m+diet,
         data=mediation_2)
summary(m1_2)

m2_2 <- lm(delta_weight_3m_2~delta_GL_3m+delta_cal_3m+diet,
         data=mediation_2)
summary(m2_2)

med_model_2 <- mediate(m0_2,m2_2,treat="diet",
                     mediator="delta_GL_3m",
                     robustSE=T,
                     sims=100)
summary(med_model_2)

####################SUPPLEMENTAL MATERIAL (Def2)############################
############################################################################
##### SUPPLEMENTAL TABLE 1 (12-months) (Def2) #####
##### SEM MEDIATION ANALYSIS (Def2) #####
mediation_12_2 <- data_12m %>% select(delta_weight_12m_2,delta_cal_12m,delta_GL_12m,
                                    diet) %>% 
  na.omit()

m0_12_2 <- lm(delta_cal_12m~delta_GL_12m+diet,
            data=mediation_12_2)
summary(m0_12_2)

m1_12_2 <- lm(delta_weight_12m_2~delta_cal_12m+diet,
            data=mediation_12_2)
summary(m1_12_2)

m2_12_2 <- lm(delta_weight_12m_2~delta_GL_12m+delta_cal_12m+diet,
            data=mediation_12_2)
summary(m2_12_2)

med_model_12_2 <- mediate(m0_12_2,m2_12_2,treat="diet",
                        mediator="delta_GL_12m",
                        robustSE=T,
                        sims=100)
summary(med_model_12_2)

############################################################################
##### SUPPLEMENTAL FIGURE 1 (Def2) #####
##### Weight change adjusted models 12m #####
w_12m_adj_totcal_2 <- lm(delta_weight_12m_2~cal+scr_gender+bl_bmi+bl_cal+diet,
                       data=data_12m)

w_12m_adj_carbcal_2 <- lm(delta_weight_12m_2~carb_cals_12m+scr_gender+bl_bmi+bl_cal
                        +diet,data=data_12m)

w_12m_adj_protcal_2 <- lm(delta_weight_12m_2~prot_cals_12m+scr_gender+bl_bmi+bl_cal
                        +diet,data=data_12m)

w_12m_adj_fatcal_2 <- lm(delta_weight_12m_2~fat_cals_12m+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_12m)

w_12m_adj_sugarcal_2 <- lm(delta_weight_12m_2~sugar.cal+scr_gender+bl_bmi+bl_cal
                         +diet,data=data_12m)

w_12m_adj_fiber_2 <- lm(delta_weight_12m_2~fiber.g+scr_gender+bl_bmi+bl_cal
                      +diet,data=data_12m)

w_12m_adj_satfat_2 <- lm(delta_weight_12m_2~saturated_fat..+scr_gender+bl_bmi+bl_cal
                       +diet,data=data_12m)

w_12m_cal_adj_comp_2 <- compare_performance(w_12m_adj_totcal_2,w_12m_adj_carbcal_2,
                                          w_12m_adj_protcal_2,w_12m_adj_fatcal_2,
                                          w_12m_adj_sugarcal_2,w_12m_adj_fiber_2,
                                          w_12m_adj_satfat_2)
w_12m_cal_adj_comp_2
plot(w_12m_cal_adj_comp_2)
ggsave("Def2_SupFig1_Adjusted models. Baseline to 12m.tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')
############################################################################
##### SUPPLEMENTAL FIGURE 2. Dietary Mediators (12-months) (Def2) #####
##### Panel A. Glycemic load reduction vs weight loss (Def2) #####
gl_12m_trend_2 <- ggscatterstats(
  data  = data_12m,
  x     = z_delta_GL_12m,
  y     = delta_weight_12m_2,
  xlab  = "GL reduction (z-score)",
  ylab  = "Weight change(kg)",
  title = "Weight change vs change in Glycemic Load",
  subtitle = "ß = -2.73, p = 6.53 x 10^3",
  caption = "12 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
gl_12m_trend_2
ggsave("Def2_SupFig2_A.Weight reduction vs Glycemic Load reduction(12m).tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

d2supfig_2a_beta_2 <- lm(delta_weight_12m_2~delta_GL_12m,
                    data=data_12m)
summary(d2supfig_2a_beta_2)

d2plot_supfig_2a_beta_2 <- ggcoefstats(d2supfig_2a_beta_2)
d2plot_supfig_2a_beta_2

##### Panel B. Fat intake reduction vs weight loss (Def2) #####
cal_12m_trend_2 <- ggscatterstats(
  data  = data_12m,
  x     = z_delta_fatcals_12m,
  y     = delta_weight_12m_2,
  xlab  = "Calories from fat reduction (z-score)",
  ylab  = "Weight change (kg)",
  title = "Weight change vs fat intake reduction",
  subtitle = "ß = 0.71, p = 0.48",
  caption = "12 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
cal_12m_trend_2
ggsave("Def2_SupFig2_B.Weight reduction vs fat intake reduction (12m).tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

d2supfig_2b_beta_2 <- lm(delta_weight_12m_2~delta_fatcals_12m,
                    data=data_12m)

summary(d2supfig_2b_beta_2)

d2beta_supfig_2b_beta_2 <- ggcoefstats(d2supfig_2b_beta_2)
d2beta_supfig_2b_beta_2

############################################################################
##### SUPPLEMENTAL FIGURE 3. Adherence measures (12-months) (Def2) #####
##### Panel A. Carb intake reduction vs weight loss in the LCD (Def2) #####
lcd_adherence_12_2 <- ggscatterstats(
  data  = data_lcd_12m,
  x     = z_delta_carbcals_12m,
  y     = delta_weight_12m_2,
  xlab  = "Change in carbohydrate intake (z-score)",
  ylab  = "Weight change (kg)",
  title = "Weight change vs change in carbohydrate intake in the LCD group",
  subtitle = "ß = -2.25, p = 0.03",
  caption = "12 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
lcd_adherence_12_2
ggsave("Def2_SupFig3_A.Weight reduction vs carbohydrate intake reduction in the LCD group (12m).tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

d2supfig_3a_beta_2 <- lm(delta_weight_12m_2~delta_carbcals_12m,
                    data=data_lcd_12m)

summary(d2supfig_3a_beta_2)

d2plot_supfig_3a_beta_2 <- ggcoefstats(d2supfig_3a_beta_2)
d2plot_supfig_3a_beta_2

##### Panel B. Fat intake reduction vs weight loss in the LFD (Def2) #####
lfd_adherence_12_2 <- ggscatterstats(
  data  = data_lfd_12m,
  x     = z_delta_fatcals_12m,
  y     = delta_weight_12m_2,
  xlab  = "Change in fat intake (z-score)",
  ylab  = "Weight change (kg)",
  title = "Weight change vs change in fat intake in the LFD group",
  subtitle = "ß = -2.0, p = 0.05",
  caption = "12 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
lfd_adherence_12_2
ggsave("Def2_SupFig3_B.Weight reduction vs fat intake reduction in the LFD group (12m).tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

d2supfig_3b_beta_2 <- lm(delta_weight_12m_2~delta_fatcals_12m,
                    data=data_lfd_12m)

summary(d2supfig_3b_beta_2)

d2plot_supfig_3b_beta_2 <- ggcoefstats(d2supfig_3b_beta_2)
d2plot_supfig_3b_beta_2

############################################################################
##### SUPPLEMENTAL FIGURE 4. Lipid Biomeasures (12-months) (Def2) #####
##### Panel A. Reduction of Tg/HDL vs weight loss (Def2) #####
tgtohdl_12m_trend_2 <- ggscatterstats(
  data  = data_12m,
  x     = z_delta_tgtohdl_12m,
  y     = delta_weight_12m_2,
  xlab  = "Change in TG/HDL (z-score)",
  ylab  = "Weight change (kg)",
  title = "Weight change vs Change in TG/HDL",
  subtitle = "ß = -8.41, p = 6.0 x 10^16",
  caption = "12 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
tgtohdl_12m_trend_2
ggsave("Def2_SupFig4_A.Weight reduction vs TG to HDL reduction (12m).tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

d2supfig_4a_beta_2 <- lm(delta_weight_12m_2~delta_tgtohdl_12m,
                    data=data_12m)

summary(d2supfig_4a_beta_2)

d2plot_fig_4a_beta_2 <- ggcoefstats(d2supfig_4a_beta_2)
d2plot_fig_4a_beta_2

##### Panel B. Reduction of LDL+HDL vs weight loss (Def2) #####
ldlplushdl_12m_trend_2 <- ggscatterstats(
  data  = data_12m,
  x     = z_delta_ldlplushdl,
  y     = delta_weight_12m_2,
  xlab  = "Change in LDL + HDL (z-score)",
  ylab  = "Weight change (kg)",
  title = "Weight change vs change in LDL + HDL",
  subtitle = "ß = 0.7, p = 0.49",
  caption = "12 month timepoint",
  results.subtitle = F,
  marginal = F
) + 
  ggplot2::scale_y_continuous(
    limits = c(-20,20),
    breaks = seq(-20, 20, by= 10))+
  ggplot2::scale_x_continuous(
    limits = c(-3,3),
    breaks = seq(-3, 3, by= 1))+
  ggplot2::theme(plot.title = element_text(size=24,face = "bold"),
                 plot.title.position = "plot",
                 plot.margin = margin(1,1,1,1,unit="cm"),
                 axis.title.y = element_text(size=20,colour = "black",face="bold",margin = margin(r=1,unit="cm")),
                 axis.title.x = element_text(size=20,colour = "black",face="bold",margin = margin(t=1,unit="cm")),
                 axis.text = element_text(size=16,colour = "black"),
                 plot.subtitle = element_text(size=20,colour = "blue",face="bold",hjust = 0.5,margin = margin(b=1,unit="cm")),
                 plot.caption = element_text(size=18,colour="gray40",face="bold.italic",hjust = 1,margin = margin(t=0.5,unit="cm"))
  ) +
  coord_cartesian(expand = F)
ldlplushdl_12m_trend_2
ggsave("Def2_SupFig4_B.Weight reduction vs LDL + HDL reduction(12m).tiff", units="cm", width=30, height=30, dpi=600, compression = 'lzw')

d2supfig_4b_beta_2 <- lm(delta_weight_12m_2~delta_ldlplushdl,
                  data=data_12m)

summary(d2supfig_4b_beta_2)

d2plot_fig_4b_beta_2 <- ggcoefstats(d2supfig_4b_beta_2)
d2plot_fig_4b_beta_2

############################################################################
##### SUPPLEMENTAL FIGURE 5. Effect Modification (12-months) (Def2) #####
##### Weight loss vs insulin at 30 min vs GL reduction (Def2) #####
w_12m_df_5q_gl_2 <- aggregate(x=data_12m$delta_weight_12m_2,
                            by= list(data_12m$qdelta_GL,data_12m$qbl_30ins),
                            FUN=mean,na.rm=TRUE)
colnames(w_12m_df_5q_gl_2) <- c("GL","qbl30min","weightloss")

w_12m_hmap_plot_5q_2 <- ggplot(data=w_12m_df_5q_gl_2,aes(x=GL,y=qbl30min))+
  geom_tile(aes(fill=weightloss))+
  ggtitle("Weight loss vs Glycemic Load Reduction and Baseline OGTT 30' insulin quintile",
          subtitle = "12 month time point")+
  scale_y_discrete(name="Baseline OGTT insulin 30' quintile",breaks=waiver())+
  scale_x_discrete(name="Change in glycemic goad quintile",breaks=waiver())+
  scale_fill_gradient(name="% change",low="#FC6B64",high="#6B77F8",n.breaks=5)+
  geom_text(aes(GL,qbl30min,label=round(weightloss,1)),color="black",size=6)+
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
w_12m_hmap_plot_5q_2
ggsave("Def2_SupFig5.Effect Modification 12-months.tiff", units="cm", width=35, height=35, dpi=600, compression = 'lzw')

##### Weight loss vs insulin at 30 min vs GL reduction (Def2) #####
em_12m_2 <- lm(delta_weight_12m_2~bl30q5*deltaGLq5+diet+bl_bmi+bl_cal,
         data=data_12m)
summary(em_12m_2)

em_12m_betas_2 <- ggcoefstats(em_12m_2)
em_12m_betas_2
#########################LETTING YOU KNOW IT FINISHED#######################
beep(8)
