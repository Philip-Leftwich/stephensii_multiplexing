###### custom function to add filename when reading data from excel spreadsheets

required <- list("tidyverse", 
                 "readxl", 
                 "ggbeeswarm",
                 "ggdark",
                 "colorspace",
                 "ggdist",
                 "lme4",
                 "lmerTest",
                 "MuMIn",
                 "ggh4x",
                "glmmTMB",
                "emmeans")

lapply(required, library, character.only = T)

library(tidyverse)
library(readxl)
library(ggbeeswarm)
library(ggdark) #desaturate()
library(colorspace) #lighten()
library(ggdist) # rainclouds
library(lme4)
library(lmerTest)
library(MuMIn) # dredge
library(ggh4x) # facet nested wrap
library(glmmTMB)


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