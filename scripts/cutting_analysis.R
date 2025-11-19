
## Functions ====
source("scripts/custom_theme.R")
source("scripts/functions.R")

# homing and mosaic data read-in====

homing_data <- read_csv("data/homing/homing_data.csv")
  

# ================================
 
# Analytical models====
 
## Homing model====
 


   homing_model <- glmmTMB(cbind(win,loss)~ gRNA_type*pre_cut+(1|group_letter/id), family=binomial, data=homing_data)
   summary(homing_model)
   sim <- simulateResiduals(homing_model)
   plot(sim, asFactor = T)
 
   

# Cutting model ====

   # Summary across observed data
homing_data %>% 
  group_by(gRNA_type) %>% 
  summarise(mean = mean((cleavage/(win+loss)), na.rm = T))

   # Precise model
cutting_model <- glmmTMB(cbind(cleavage,((win+loss)-cleavage))~ gRNA_type+(1|group_letter/id), family=binomial, data=homing_data)
summary(cutting_model)
sim <- simulateResiduals(cutting_model)
plot(sim, asFactor = T)

##  Robust model ====
# Wide confidence intervals for 338-384 due to complete homing rate
robust_cutting_model <- brglm(cbind(cleavage,((win+loss)-cleavage))~ gRNA_type, family=binomial, data=homing_data)

emmeans::emmeans(robust_cutting_model, specs = ~ gRNA_type, type = "response")

means <- emmeans::emmeans(robust_cutting_model, specs = ~ gRNA_type)
pairs(means)
   



                  