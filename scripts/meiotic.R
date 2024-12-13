source("scripts/functions.R")
source("scripts/custom_theme.R")
library(DHARMa)
library(lmerTest)
library(glmmTMB)
font_add_google("Open Sans", "Sans")

### meiotic driver

file_paths <- list.files(path="./data/meiotic", pattern = "*.xlsx", full.names = T) 
sheet_nums <- c(2,3)

meiotic_df <- map2_df(file_paths, sheet_nums, ~read_plus(.x, sheet = .y, range = ("C4:L32")))


####

df_2301 <- meiotic_df %>% 
  filter(filename == "./data/meiotic/Crossing data summary A2301 B_1590 M_1928B12 x QA383P - CLEAN.xlsx") %>% 
  mutate(win = (ABM+AM),
       loss = (AB+A)) 

meiotic_model <- glmmTMB(cbind(win,loss)~ 1+(1|`Female no.`), family=binomial, data=df_2301)

sim <- simulateResiduals(meiotic_model)
plot(sim, asFactor = T)

emmeans::emmeans(meiotic_model, specs = ~1, type = "response") 



####

df_1759 <- meiotic_df %>% 
  filter(filename == "./data/meiotic/Crossing data summary A_1759 B_1590 Marker_1928B12 - CLEAN.xlsx") %>% 
  drop_na(`A`) %>% 
  mutate(win = (ABM+AM),
         loss = (AB+A)) 
  
  meiotic_model2 <- glm(cbind(win,loss)~ 1, family=binomial, data=df_1759)

  emmeans::emmeans(meiotic_model2, specs = ~1, type = "response") 
  