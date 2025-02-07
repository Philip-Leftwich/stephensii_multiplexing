
source("scripts/functions.R")
source("scripts/custom_theme.R")
library(DHARMa)
library(lmerTest)
library(glmmTMB)
library(ggcorrplot)
#library(ggseqlogo)

# Data

files <- list.files(path = "data/indels", pattern = "\\.xlsx$", full.names = TRUE)


data <- map_df(files, ~read_excel(.x, range = "A1:AM403"))


#data1 <- data %>%
#  filter(nchar(Aligned_Sequence) == 471) %>% 
#  mutate(Aligned_Sequence = substr(Aligned_Sequence, 201, nchar(Aligned_Sequence)))



# ggseqlogo( data1$Aligned_Sequence, method = 'prob', facet = "wrap")



# Check tally of indels against read count
data %>%
  select(., contains("read")) %>% 
  mutate(total = rowSums(across(4:18)), .after = `#Reads`)

# 
data <- data %>%
  select(., contains("read")) %>% 
  mutate(id = 1:nrow(.),.before = `Read_Status`) %>% 
  mutate(`guide_1` = rowSums(across(contains("1"))),
         `guide_2` = rowSums(across(contains("2"))),
         `guide_3` = rowSums(across(contains("3"))),
         `guide_4` = rowSums(across(contains("4"))),
  .after = `#Reads`) %>% 
  mutate(`guide_1` = if_else(`guide_1` > 0, 1, 0),
          `guide_2` = if_else(`guide_2` > 0, 1, 0),
          `guide_3` = if_else(`guide_3` > 0, 1, 0),
          `guide_4` = if_else(`guide_4` > 0, 1, 0),  
         .after = `#Reads`) 

data_modified <- data %>% filter(`Read_Status` == "MODIFIED")


# Minor discrepancy in data? some guides not marked for Modified? Check this
data %>%
  group_by(`Read_Status`) %>%
  summarise(sum = sum(`#Reads`))

data %>%
  filter(guide_1 == 1 | guide_2 ==1 | guide_3 == 1 | guide_4 ==1) %>% 
  summarise(sum = sum(`#Reads`))

data %>%
  filter(guide_1 == 1 & guide_2 ==1 & guide_3 == 1 & guide_4 ==1) %>% 
  summarise(sum = sum(`#Reads`))



# Binomial model

data_bin <- data %>% 
  group_by(guide_1, guide_2, guide_3, guide_4) %>% 
  summarise(`#Reads` = sum(`#Reads`))


# Ggplot number of reads by missing guide frequency

Perc_bar <- data_bin %>%  
  ungroup %>% 
  mutate(sum = rowSums(across(starts_with("guide")))) %>% 
  mutate(sum = factor(sum)) %>% 
  group_by(sum) %>% 
  filter(sum != 0) %>% 
  summarise(reads = sum(`#Reads`)) %>% 
  mutate(perc = reads/sum(reads)) %>% 
  #  mutate(sum = fct_reorder(sum, reads)) %>% 
  ggplot(aes(x = sum, y = perc))+
  geom_col(fill = "lightblue", colour = "darkgrey")+ 
  coord_flip()+
  geom_label(aes(label = scales::percent(perc)), nudge_y = .05)+
  scale_y_continuous(labels = scales::percent, limits = c(0,1))+
  theme_custom()+
  labs(x = "No. of gRNA sites lost",
       y = "Percentage of reads")

ggsave("figures/Percentages.png", dpi = 900, width = 9, height = 6, units = "in")

bin_model <- glm(cbind(`#Reads`, (23958 - `#Reads`)) ~ guide_1 * guide_2 * guide_3 * guide_4, family = binomial, data = data_bin)


emmeans::emmeans(bin_model, specs = ~ guide_1 + guide_2 + guide_3 + guide_4, type = "response") %>% as_tibble() %>% arrange(desc(response))



# Modified reads only

data_bin <- data_modified %>% 
  group_by(guide_1, guide_2, guide_3, guide_4) %>% 
  summarise(`#Reads` = sum(`#Reads`))

bin_model_mod <- glm(cbind(`#Reads`, (10529 - `#Reads`)) ~ as.factor(guide_1) * as.factor(guide_2) * as.factor(guide_3) * as.factor(guide_4), family = binomial, data = data_bin, na.action = "na.fail")

emmeans::emmeans(bin_model_mod, specs = ~ guide_1 + guide_2 + guide_3 + guide_4, type = "response") %>% as_tibble() %>% arrange(desc(response))

dd <- dredge(bin_model_mod)
subset(dd, delta < 4)

#new_model <- model.avg(subset(dd, delta < 4))

#emmeans::emmeans(new_model, specs = ~ guide_1 + guide_2 + guide_3 + guide_4, type = "response") %>% as_tibble() %>% arrange(desc(response))

### Attempting new analytical methods

cor_matrix <- cor(data[, c("guide_1", "guide_2", "guide_3", "guide_4")])

cor_matrix %>%   
  ggcorrplot(hc.order = TRUE, type = "upper",
                            lab = TRUE)

guide_cols <- c("guide_1", "guide_2", "guide_3", "guide_4")
read_col <- "#Reads"  # Assuming the column for the number of reads is named "#Reads"
df_subset <- data_modified[, guide_cols]
weights <- data_modified[,read_col]

weighted_means <- data_modified %>% mutate(across(guide_1:guide_4,  ~ . * `#Reads`)) %>% 
  summarise(across(guide_1:guide_4, ~(sum(.)/sum(`#Reads`))))

library(weights)
weighted_correlation_matrix <- wtd.cors(df_subset, weight = weights)

weighted_correlation_matrix %>%   
  ggcorrplot(hc.order = TRUE, type = "upper", 
             lab = TRUE)

## Venn



data_long <- data_bin %>% 
  uncount(weights = `#Reads`) %>% 
  rename("sgRNA 338 \n 8578" = "guide_1",
         "sgRNA 347 \n 4055"= "guide_2",
        "sgRNA 362\n 4516" = "guide_3",
         "sgRNA 384 \n 4204" = "guide_4")

ggVennDiagram(lapply(data_long, function(x) which(x == 1))) + 
  scale_fill_viridis(option = "cividis")+
  labs(fill="Read count") +
  theme_void(base_size = 10)

# ggVennDiagram(lapply(data_long, function(x) which(x == 1)), force_upset = T, sets.bar.color = c("lightblue", "blue", "blue", "blue"))

ggsave("figures/Venn.png", dpi = 900, width = 12, height = 8, units = "in") 

