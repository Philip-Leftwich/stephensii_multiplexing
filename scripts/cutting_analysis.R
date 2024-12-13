
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

### 03.10.23====


### 04.10.23====

#new <- result_df %>% 
#  fill(`F1 cross`) %>% 
#  mutate(phenotype = case_when(
#    str_detect(`F1 cross`, "\\bME\\b") ~ "Mosaic",
#    str_detect(`F1 cross`, "\\bDE\\b") ~ "Dark",
#    .default = "other")) %>% 
  
  
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
     mutate(pre_cut = factor(pre_cut, labels = c(
                               "wild-type" = "wild-type",                          # Plain text
                               "D251"      = "cd<sup><i>225</i></sup>",           # cd with superscript italic "225"
                               "QA383P"    = "cd<sup><i>384</i></sup>"            # cd with superscript italic "384"
                             ))) %>% 
     mutate(gRNA_type = case_when(
       gRNA_type=="2072" ~ "cd<sup><i>g338_del</i></sup>",
       gRNA_type=="1759" ~ "cd<sup><i>g384</i></sup>",
       gRNA_type=="2301" ~ "cd<sup><i>g338-384</i></sup>",
       gRNA_type=="2273" ~ "cd<sup><i>g225</i></sup>")) 
  



#======================
   
# All crosses
   
#   simple_model <- glmer(cbind(win,loss)~ cas9_parent*gRNA_type*pre_cut+(1|group_letter/id), family=binomial, data=homing_data)
#   
#   second_model <- glmmTMB(cbind(win,loss)~ group_letter+(1|group_letter/id), family=binomial, data=homing_data)
#   sim <- simulateResiduals(second_model)
#   plot(sim, asFactor = T)
   
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
   

   

   
   
   homing_plot <-   complex_data %>% 
     mutate(homing=((win/(win+loss)*100))) %>% 
     ggplot(aes(x=gRNA_type, y=homing, group=id))+
     geom_quasirandom(aes(size=win+loss), 
                      fill="white", 
                      shape = 21,
                      alpha=0.6,
                      colour = colour)+
     geom_errorbar(data = homing_summary, colour = colour, aes(min=(asymp.LCL*100), max=(asymp.UCL*100), y=prob, group=NA),width=0,  linewidth=1.2, position=position_nudge(x=0.4))+
     geom_point(data=homing_summary, colour = colour, aes(y=prob*100, group=NA, fill=after_scale(desaturate(lighten(colour, .6), .6))), size=3, position=position_nudge(x=0.4), stroke=0.6, shape = 21)+
     scale_size(range=c(0,3),
                breaks=c(50,100,150))+
     
     guides(fill=FALSE, colour=FALSE)+
     labs(x="", 
          y="Percentage of individuals scored",
          size="Number \nof \noffspring",
          shape="")+
     scale_x_discrete(guide = guide_axis(n.dodge = 2))+
     scale_y_continuous(limits=c(10,100),
                        labels=scales::percent_format(scale=1) # automatic percentages
     )+
     
     theme_custom()+
     theme(strip.text = element_markdown(),
           axis.text.x = element_markdown())+
     facet_wrap(~ pre_cut, scales = "free")
     
   
 
ggsave("homing_plot.png", dpi = 900, width = 9, height = 6, units = "in")
   
   #================================
   
   # 1759 effects  
   
   data_1759 <- homing_data %>% 
     filter(cas9_parent == "Male") %>% 
     filter(gRNA_type == "1759") 
   
 #  model_1759 <- glmmTMB(cbind(win,loss)~ pre_cut+(1|group_letter/id), family=binomial, data=data_1759)
 #  summary(model_1759)
 #  sim <- simulateResiduals(model_1759)
 #  plot(sim, asFactor = T)  
   
   
   
#=====================
   
#cutting_data <- cutting_tbl %>% 
#  rename(`F1 cross`=1,
 #        genotype_A=4,
  #       genotype_AB=5,
   #      genotype_B=6,
    #     genotype_WT=7,
     #    mosaic_A=15,
       #  mosaic_AB=16,
      #   mosaic_B=17,
        # mosaic_WT=18) %>% 
#  select(c(1:7, 15:18, 19)) %>% 
#  fill(`F1 cross`) %>% 
#  drop_na(genotype_A:mosaic_WT) %>% 
#  mutate(cross= str_remove_all(`F1 cross`, "[()Xx]")) %>% ### keep F1 cross variable
#  mutate(cross = str_remove_all(cross, "Cd ")) %>% 
#  separate(cross, into =c("male_parent", "female_parent"), 
#           sep=" ", extra = "merge", fill = "left") %>% 
#  mutate(Cas9_parent=if_else(male_parent=="KO", "Female", "Male")) %>% 
#  separate(male_parent, into=c("m_male_grandparent", "m_female_grandparent"), sep=":", remove=FALSE) %>% 
#  separate(female_parent, into=c("f_male_grandparent", "f_female_grandparent"), sep=":", remove=FALSE) %>% 
 # mutate(m_male_grandparent=str_trim(m_male_grandparent, side = "both")) %>% # extra white spaces causing issues with ifelse statement
#  mutate(f_male_grandparent=str_trim(f_male_grandparent, side = "both")) %>% 
  #mutate(Cas9_grandparent=ifelse(m_male_grandparent=="1590B", "Male",
   #                              ifelse(f_male_grandparent=="1590B", "Male", "Female"))) %>% 
 # rename("Line"=filename)

### Models


#cutting_data_model <- cutting_data %>% 
#  mutate(win = (genotype_A+genotype_AB),
#         loss = (genotype_B+genotype_WT),
#         win1=(mosaic_A+mosaic_AB),
#         loss1=((win+loss)- win1)) %>% 
#  mutate(id=row_number())

# use the same parameters as the models for crosses to WT (homing_mosaic_analysis.R)

#cut.model.1 <- glmmTMB((cbind(win,loss))~ Cas9_parent+Cas9_grandparent+Line+Cas9_parent:Cas9_grandparent+Cas9_parent:Line+Cas9_grandparent:Line+(1|Line/Cas9_grandparent/Cas9_parent/id), family=binomial, data=cutting_data_model, REML=FALSE)

#sim <- DHARMa::simulateResiduals(cut.model.1) 
#plot(sim, asFactor=T)

#mosaic.model.1 <- glmmTMB((cbind(win1,loss1))~ Cas9_parent+Cas9_grandparent+Line+Cas9_parent:Cas9_grandparent+Cas9_parent:Line+Cas9_grandparent:Line+(1|Line/Cas9_grandparent/Cas9_parent/id), family=binomial, data=cutting_data_model, REML=FALSE)#

#sim <- DHARMa::simulateResiduals(mosaic.model.1) 
#plot(sim, asFactor=T)




### Summary tables

#gtsummary::tbl_regression(model2.fit, exponentiate=TRUE)
#sjPlot::tab_model(model2.fit)



### Figure

# cutting_summary <- emmeans::emmeans(cut.model.1, specs=pairwise~Cas9_parent+Cas9_grandparent+Line+Cas9_parent:Cas9_grandparent+Cas9_parent:Line+Cas9_grandparent:Line, type="response") %>% .$emmeans %>% as_tibble() %>% bind_cols(homing_or_cutting="cutting")

# homing_summary <- emmeans::emmeans(mosaic.model.1, specs=pairwise~Cas9_parent+Cas9_grandparent+Line+Cas9_parent:Cas9_grandparent+Cas9_parent:Line+Cas9_grandparent:Line, type="response") %>% .$emmeans %>% as_tibble() %>% bind_cols(homing_or_cutting="homing")

# total_CI <- rbind(cutting_summary, homing_summary) %>% 
#  mutate(Cas9_parent = factor(Cas9_parent, levels=c("Female", "Male"), labels=c("\u2640", "\u2642"))) %>%
#  mutate(Cas9_grandparent = factor(Cas9_grandparent, levels=c("Female", "Male"), labels=c("\u2640", "\u2642")))


# cutting_data_model %>% 
#   mutate(cutting=((win1/(win+loss))*100)) %>% 
#   mutate(homing=((win/(win+loss)*100))) %>% 
#   mutate(Cas9_parent = factor(Cas9_parent, levels=c("Female", "Male"), labels=c("\u2640", "\u2642"))) %>%
#   mutate(Cas9_grandparent = factor(Cas9_grandparent, levels=c("Female", "Male"), labels=c("\u2640", "\u2642"))) %>%
#   pivot_longer(cols=homing:cutting, names_to="homing_or_cutting", values_to="count") %>% 
#   ggplot(aes(x=homing_or_cutting, y=count, group=id, colour=Cas9_parent, shape=homing_or_cutting))+
#   geom_point(aes(size=win+loss), fill="white", alpha=0.6)+
#   geom_line(alpha=0.4, size=0.1)+
#  geom_errorbar(data=total_CI, aes(min=(lower.CL*100), max=(upper.CL*100), y=response, group=NA),width=0,  size=1.2, position=position_nudge(x=0.2))+
#  geom_point(data=total_CI, aes(y=response*100, group=NA, fill=after_scale(desaturate(lighten(colour, .6), .6))), size=2, position=position_nudge(x=0.2), stroke=0.6)+
#  scale_shape_manual(values=c(21,24), labels=c("white eyes","sgRNAs"))+
#  scale_size(range=c(0,2),
#             breaks=c(50,100,150))+
#  scale_color_manual(values=c("#a86826", "#006c89"))+
#  guides(fill=FALSE, colour=FALSE)+
#  labs(x="", 
#       y="Percentage of individuals scored",
#       size="Number of offspring",
#       shape="")+
#  scale_y_continuous(limits=c(50,100),
#    labels=scales::percent_format(scale=1) # automatic percentages
#  )+
#  facet_nested_wrap(~Line+Cas9_grandparent+Cas9_parent, 
#                    nrow=1,
#                    strip.position="bottom",
#                    scales="free_x")+ # reduce lines
 # theme_custom()+
 # theme(axis.text.x=element_blank(),
  #      legend.position=c(0.85,0.15),
   #     strip.background = element_blank(),
    #    strip.placement="outside",
     #   panel.spacing = unit(0, "lines"),
      #  strip.text.x = element_text(margin = margin(5,0,5,0),
       #                             hjust = 0.5))# center strip text


# ggsave("figures/cutting.png", dpi=600,  width=6, height=6)



### Pretty table


                  