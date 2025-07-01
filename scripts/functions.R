###### custom function to add filename when reading data from excel spreadsheets

required <- list("tidyverse", 
                 "readxl", 
                 "ggbeeswarm",
                 "ggdark",
                 "colorspace",
                 "ggdist",
                 "MuMIn",
                 "ggh4x",
                "glmmTMB",
                "emmeans",
                "gt",
                "ggVennDiagram",
                "viridis",
                "ggtext",
                "DHARMa",
                "ggcorrplot")

lapply(required, library, character.only = T)

library(brglm) # robust model for perfect convergence
font_add_google("Open Sans", "Sans")



## Functions

read_plus <- function(flnm, sheet, skip, range) {
  read_excel(flnm ,sheet=sheet, skip=skip, na=c("na","NA", "-"), range = range) %>% 
    mutate(filename = flnm)
  
}

#####

##### force bind, bind rows of data when column names do not match #####

force_bind = function(df1, df2) {
  colnames(df2) = colnames(df1)
  bind_rows(df1, df2)
}

#####


##### DHARMa_check function to simulate residuals and plot - for mixed models

DHARMa_check <- function(model){
  sim <- DHARMa::simulateResiduals(model) 
 plot(sim, asFactor=T)
}


##### binned plot to check overdispersion in binomial model

bin_plot <- function(model){
x <- predict(model)
y <- resid(model)
arm::binnedplot(x,y)
}


### Homing plots

homing_plots <- function(.df, .means, .colours="grey10", facet = TRUE){
  
  homing_plot <- .df %>% 
    mutate(homing=((win/(win+loss)*100))) %>% 
    ggplot(aes(x=gRNA_type, 
               y=homing, 
               group=id,
               colour = gRNA_type))+
    geom_quasirandom(aes(size=win+loss), 
                     fill="white", 
                     shape = 21,
                     alpha=0.8)+
    geom_errorbar(data = .means, aes(min=(asymp.LCL*100), max=(asymp.UCL*100), y=prob, group=NA),width=0,  linewidth=1.2, position=position_nudge(x=0.4))+
    geom_point(data= .means, aes(y=prob*100, group=NA, fill=after_scale(desaturate(lighten(colour, .4), 0.1))), size=3, position=position_nudge(x=0.4), stroke=1.1, shape = 21)+
    geom_label(data = .means, aes(y=115, group=NA, label = paste0("n = ",n), fill= gRNA_type), colour = "black", size=5)+
    scale_size(range=c(0,3),
               breaks=c(50,100,150))+
    guides(fill=FALSE, colour=FALSE)+
    labs(x="", 
         y="Percentage of individuals scored",
         size="Number \nof \noffspring",
         shape="")+
    scale_x_discrete()+
    scale_y_continuous(limits=c(10,120),
                       breaks = c(0,25,50,75,100),
                       labels=scales::percent_format(scale=1) # automatic percentages
    )+
    scale_colour_manual(values = .colours)+
    scale_fill_manual(values = after_scale(desaturate(lighten(.colours, .8), .2)))+
    theme_custom()+
    theme(strip.text = element_markdown(),
          axis.text.x = element_markdown())
  
  
  if(facet == TRUE){
    homing_plot <- homing_plot+
      facet_wrap(~ pre_cut, scales = "free")
  }
  print(homing_plot)
}
