

source("scripts/custom_theme.R")
source("scripts/functions.R")

# homing and mosaic data read-in====

file_paths <- list.files(path="./data/", pattern = "*.xlsx", full.names = T) 
sheet_nums <- c(3,2,2,3,3,3,2,2)

result_df <- map2_df(file_paths, sheet_nums, ~read_plus(.x, sheet = .y, range = ("C4:I300")))%>% 
    fill(`F1 cross`) %>% 
  drop_na(`Female no.`)

cut_result <- map2_df(file_paths, sheet_nums, ~read_plus(.x, sheet = .y, range = ("C4:AA300")))%>% 
  fill(`F1 cross`) %>% 
  drop_na(`Female no.`) %>% 
  select(`...25`) %>%
  rename("cleavage"=`...25`) 

result_df <- cbind(result_df, cut_result)
  
# Create a data frame of unique strings and assign group letters
 unique_strings <- result_df %>%
      distinct(`F1 cross`) %>%
       mutate(group_letter = letters[1:n()])

# Merge the original data frame with the unique strings data frame
   new <- result_df %>%
       left_join(unique_strings, by = "F1 cross")
   
   homing_data <- new %>% 
     mutate(cas9_parent = if_else(group_letter %in% c("a","b","d","e","g","i","k","l","m","p","q","r","s"), "Male", "Female")) %>% 
     mutate(gRNA_type = case_when(group_letter %in% c("a","b","c","d","e","f","g","h","i","j","k") ~ "1759",
                                  group_letter %in% c("l", "m", "n", "o") ~ "2072",
                                  group_letter %in% c("p", "q") ~ "2273",
                                  group_letter %in% c("r", "s") ~ "2301")) %>% 
     mutate(target_locus = if_else(gRNA_type %in% c("1759", "2072", "2301"), "locus A", "locus B")) %>% 
     mutate(multiplex = if_else(gRNA_type == "2301", "yes", "no")) %>% 
     mutate(pre_cut = case_when(group_letter %in% c("g","h","i","j") ~ "D251",
                                group_letter %in% c("k","q","s") ~ "QA383P",
                                .default = "wild-type"))%>% 
     mutate(pre_cut = factor(pre_cut, levels = c("wild-type", "D251", "QA383P"))) %>% 
     mutate(mosaic = case_when(group_letter %in% c("a", "c", "e") ~ "yes", 
                               group_letter %in% c("b", "d","f") ~ "no",
                               .default = "mixed"))%>% 
     mutate(id=row_number())%>% 
     mutate(win = (A+AB),
            loss = (B+WT)) %>% 
     mutate(pre_cut_original = pre_cut) %>% 
     mutate(pre_cut = factor(pre_cut, labels = c(
                               "wild-type" = "wild-type",                          # Plain text
                               "D251"      = "cd<sup><i>225</i></sup>",           # cd with superscript italic "225"
                               "QA383P"    = "cd<sup><i>384</i></sup>"            # cd with superscript italic "384"
                             ))) %>% 
     mutate(gRNA_type_original = gRNA_type) %>% 
     mutate(gRNA_type = case_when(
       gRNA_type=="2072" ~ "cd<sup><i>g384_del</i></sup>",
       gRNA_type=="1759" ~ "cd<sup><i>g384</i></sup>",
       gRNA_type=="2301" ~ "cd<sup><i>g338-384</i></sup>",
       gRNA_type=="2273" ~ "cd<sup><i>g225</i></sup>")) %>% 
     mutate(cleavage = if_else(gRNA_type_original == "2273", cleavage + win, cleavage)) # 2273 has cd recoding
  


#=======================
   
  # Mosaicism vs Dark eyes
   
   mosaic_data <- homing_data %>% 
     filter(mosaic %in% c("yes","no"))
   

#================================
   
   # Different sgRNAs
   
   complex_data <- homing_data %>% 
     filter(cas9_parent == "Male") 

   complex_model <- glmmTMB(cbind(win,loss)~ gRNA_type*pre_cut+(1|group_letter/id), family=binomial, data=complex_data)
   summary(complex_model)
   sim <- simulateResiduals(complex_model)
   plot(sim, asFactor = T)
 
   
   ## Cutting model
   
   complex_data %>% group_by(gRNA_type) %>% summarise(mean = mean((cleavage/(win+loss)), na.rm = T))
   
   cutting_model <- glmmTMB(cbind(cleavage,((win+loss)-cleavage))~ gRNA_type+(1|group_letter/id), family=binomial, data=complex_data)
   summary(cutting_model)
   sim <- simulateResiduals(cutting_model)
   plot(sim, asFactor = T)
   
   
   ## Robust glm
   
   complex_data2 <- complex_data %>% 
     mutate(gRNA_type = factor(gRNA_type, levels = c("cd<sup><i>g338-384</i></sup>",
                                                     "cd<sup><i>g384</i></sup>",
                                                     "cd<sup><i>g384_del</i></sup>",
                                                     "cd<sup><i>g225</i></sup>")))
   
   cutting_model <- brglm(cbind(cleavage,((win+loss)-cleavage))~ gRNA_type, family=binomial, data=complex_data2)
   
   ## Means
   emmeans::emmeans(cutting_model, specs = ~gRNA_type, type  ="response")
   
   emmeans::emmeans(cutting_model, specs = pairwise~gRNA_type, adjust = "none")
   
   complex_data3 <- complex_data %>% 
     mutate(gRNA_type = factor(gRNA_type, levels = c("cd<sup><i>g384</i></sup>",
                                                     "cd<sup><i>g338-384</i></sup>",
                                                     "cd<sup><i>g384_del</i></sup>",
                                                     "cd<sup><i>g225</i></sup>")))
   
   cutting_model <- brglm(cbind(cleavage,((win+loss)-cleavage))~ gRNA_type, family=binomial, data=complex_data3)
#===================================
   
  ## Figures
   
   homing_summary <- emmeans::emmeans(complex_model, specs= ~ gRNA_type*pre_cut, type="response") %>% as_tibble() 
   
   count_data <- complex_data %>% 
     group_by(gRNA_type, pre_cut) %>% 
     summarise(n = n())
   
   homing_summary <- full_join(homing_summary, count_data, by = c("gRNA_type", "pre_cut"))

   plot_title <- "Homing rates"
   
   custom_labels <- c(
     "2072" = expression(cd^italic("g338_del")),
     "1759" = expression(cd^italic("g384")),
     "2301" = expression(cd^italic("g338-384")),
     "2273" = expression(cd^italic("g225"))
     )
   
   `g384_del` <- "#fcb001"
   `g338-384` <- "#A5AB81"
   `g384` <- "#94B6D2"
   `g225` <- "#BF9000"
   
   
 colours_fig_1b <- c(`g338-384`, `g384_del`)
 colours_fig_2c <- c(`g338-384`, `g384`)
 colours_fig_3b <- c(`g384`)
 colours_fig_3d <- c(`g225`)

   fig_1b <- complex_data %>% 
     filter(gRNA_type == "cd<sup><i>g384_del</i></sup>" | gRNA_type == "cd<sup><i>g338-384</i></sup>") %>% 
     filter(pre_cut == "wild-type")
   
   fig_1b_mean <- homing_summary %>% 
     filter(gRNA_type == "cd<sup><i>g384_del</i></sup>" | gRNA_type == "cd<sup><i>g338-384</i></sup>") %>% 
     filter(pre_cut == "wild-type")
   
   fig_2c <- complex_data %>% 
     filter(gRNA_type == "cd<sup><i>g384</i></sup>" | gRNA_type == "cd<sup><i>g338-384</i></sup>") %>% 
     filter(pre_cut == "cd<sup><i>384</i></sup>")
   
   fig_2c_mean <- homing_summary %>% 
     filter(gRNA_type == "cd<sup><i>g384</i></sup>" | gRNA_type == "cd<sup><i>g338-384</i></sup>") %>% 
     filter(pre_cut == "cd<sup><i>384</i></sup>")
   
   fig_3b <- complex_data %>% 
     filter(gRNA_type == "cd<sup><i>g384</i></sup>") %>% 
     filter(pre_cut == "cd<sup><i>225</i></sup>" | pre_cut == "wild-type")
   
   fig_3b_mean <- homing_summary %>% 
     filter(gRNA_type == "cd<sup><i>g384</i></sup>") %>% 
     filter(pre_cut == "cd<sup><i>225</i></sup>" | pre_cut == "wild-type")
   
   fig_3d <- complex_data %>% 
     filter(gRNA_type == "cd<sup><i>g225</i></sup>") %>% 
     filter(pre_cut == "cd<sup><i>384</i></sup>" | pre_cut == "wild-type")
   
   fig_3d_mean <- homing_summary %>% 
     filter(gRNA_type == "cd<sup><i>g225</i></sup>") %>% 
     filter(pre_cut == "cd<sup><i>384</i></sup>" | pre_cut == "wild-type")
   

ggsave("figures/homing_plot.png", dpi = 900, width = 9, height = 6, units = "in")
   

   #================================
   
   # 1759 effects  
   
   data_1759 <- homing_data %>% 
     filter(cas9_parent == "Male") %>% 
     filter(gRNA_type == "1759") 
   



                  