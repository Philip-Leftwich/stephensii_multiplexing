
source("scripts/functions.R")
source("scripts/custom_theme.R")
library(DHARMa)
library(lmerTest)
library(glmmTMB)
font_add_google("Open Sans", "Sans")

### homing and mosaic data

file_paths <- list.files(path="./data/", pattern = "*.xlsx", full.names = T) 
sheet_nums <- c(3,2,2,3,3,3,2,2)

result_df <- map2_df(file_paths, sheet_nums, ~read_plus(.x, sheet = .y, range = ("C4:I300")))%>% 
    fill(`F1 cross`) %>% 
  drop_na(`Female no.`)




  
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
       gRNA_type=="2273" ~ "cd<sup><i>g225</i></sup>")) 
  



   
#=======================
   
  # Mosaicism vs Dark eyes
   
   mosaic_data <- homing_data %>% 
     filter(mosaic %in% c("yes","no"))
   
   # No significant effect of mosaicism on homing rates
   mosaic_model <- glmmTMB(cbind(win,loss)~ cas9_parent*mosaic+(1|group_letter/id), family=binomial, data=mosaic_data)
   summary(mosaic_model)
   sim <- simulateResiduals(mosaic_model)
   plot(sim, asFactor = T)
     
#================================
   
   # Different sgRNAs
   
   complex_data <- homing_data %>% 
     filter(cas9_parent == "Male") 
 #    filter(!(gRNA_type == "1759" & pre_cut == "QA383P")) %>% 
#    filter(pre_cut != "D251")
   
   complex_model <- glmmTMB(cbind(win,loss)~ gRNA_type*pre_cut+(1|group_letter/id), family=binomial, data=complex_data)
   summary(complex_model)
   sim <- simulateResiduals(complex_model)
   plot(sim, asFactor = T)
 
   
#===================================
   
  ## Figures
   
   homing_summary <- emmeans::emmeans(complex_model, specs= ~ gRNA_type*pre_cut, type="response") %>% as_tibble() 
   
   colour <- "grey10"
   plot_title <- "Homing rates"
   
   custom_labels <- c(
     "2072" = expression(cd^italic("g338_del")),
     "1759" = expression(cd^italic("g384")),
     "2301" = expression(cd^italic("g338-384")),
     "2273" = expression(cd^italic("g225"))
     )
   

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
   
   homing_plots <- function(.df, .means, facet = TRUE){
     
     homing_plot <- .df %>% 
     mutate(homing=((win/(win+loss)*100))) %>% 
     ggplot(aes(x=gRNA_type, y=homing, group=id))+
     geom_quasirandom(aes(size=win+loss), 
                      fill="white", 
                      shape = 21,
                      alpha=0.6,
                      colour = colour)+
     geom_errorbar(data = .means, colour = colour, aes(min=(asymp.LCL*100), max=(asymp.UCL*100), y=prob, group=NA),width=0,  linewidth=1.2, position=position_nudge(x=0.4))+
     geom_point(data= .means, colour = colour, aes(y=prob*100, group=NA, fill=after_scale(desaturate(lighten(colour, .6), .6))), size=3, position=position_nudge(x=0.4), stroke=0.6, shape = 21)+
     scale_size(range=c(0,3),
                breaks=c(50,100,150))+
     guides(fill=FALSE, colour=FALSE)+
     labs(x="", 
          y="Percentage of individuals scored",
          size="Number \nof \noffspring",
          shape="")+
     scale_x_discrete()+
     scale_y_continuous(limits=c(10,100),
                        labels=scales::percent_format(scale=1) # automatic percentages
     )+
     theme_custom()+
     theme(strip.text = element_markdown(),
           axis.text.x = element_markdown())
     
     
     if(facet == TRUE){
       homing_plot <- homing_plot+
        facet_wrap(~ pre_cut, scales = "free")
     }
     print(homing_plot)
   }
   
ggsave("figures/homing_plot.png", dpi = 900, width = 9, height = 6, units = "in")
   
   #================================
   
   # 1759 effects  
   
   data_1759 <- homing_data %>% 
     filter(cas9_parent == "Male") %>% 
     filter(gRNA_type == "1759") 
   



                  