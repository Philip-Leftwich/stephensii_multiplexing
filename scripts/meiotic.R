source("scripts/functions.R")
source("scripts/custom_theme.R")


# Read in data to check for meiotic drive====

file_paths <- list.files(path="./data/meiotic", pattern = "*.xlsx", full.names = T) 
sheet_nums <- c(2,3)

meiotic_df <- map2_df(file_paths, sheet_nums, ~read_plus(.x, sheet = .y, range = ("C4:L32")))



# cdg338-384====
df_2301 <- meiotic_df %>% 
  filter(filename == "./data/meiotic/Crossing data summary A2301 B_1590 M_1928B12 x QA383P - CLEAN.xlsx") %>% 
  mutate(win = (ABM+AM),
       loss = (AB+A)) 

meiotic_model <- glmmTMB(cbind(win,loss)~ 1+(1|`Female no.`), family=binomial, data=df_2301)

sim <- simulateResiduals(meiotic_model)
plot(sim, asFactor = T)

emmeans::emmeans(meiotic_model, specs = ~1, type = "response") 



# cdg384====
df_1759 <- meiotic_df %>% 
  filter(filename == "./data/meiotic/Crossing data summary A_1759 B_1590 Marker_1928B12 - CLEAN.xlsx") %>% 
  drop_na(`A`) %>% 
  mutate(win = (ABM+AM),
         loss = (AB+A)) 
  
  meiotic_model2 <- glm(cbind(win,loss)~ 1, family=binomial, data=df_1759)

  emmeans::emmeans(meiotic_model2, specs = ~1, type = "response") 
  
  
# Overall model====
  
  meiotic_df <- meiotic_df %>% drop_na(`A`) %>% 
    mutate(win = (ABM+AM+BM+M),
           loss = (AB+A+B+WT)) %>% 
    mutate(
      `gRNA_type` = if_else(
        str_detect(filename, "2301"),
        "cd<sup><i>g338-384</i></sup>",
        "cd<sup><i>g384</i></sup>"
      )
    ) %>% 
    mutate(filename = factor(`gRNA_type`)) %>% 
    mutate(filename = fct_rev(`gRNA_type`))
  
  meiotic_model3 <- glm(cbind(win,loss)~ `gRNA_type`, family=binomial, data=meiotic_df)
  
  homing_model <- glm(cbind(A+AB+ABM+AM,BM+B+M+WT)~ `gRNA_type`, family=binomial, data=meiotic_df)
  
  emmeans::emmeans(meiotic_model3, specs = ~ `gRNA_type`, type = "response") 
  
  emmeans::emmeans(homing_model, specs = ~ `gRNA_type`, type = "response") 
  