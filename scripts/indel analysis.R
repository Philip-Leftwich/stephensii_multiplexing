
source("scripts/functions.R")
source("scripts/custom_theme.R")



# Read in data====

files <- list.files(path = "data/indels", pattern = "\\.xlsx$", full.names = TRUE)


data <- map_df(files[2], ~read_excel(.x, range = "A1:AM29440")) %>% 
  mutate(id = row_number()+1)

data %>%
  group_by(`Read_Status`) %>%
  summarise(sum = sum(`#Reads`))

data_modified <- data %>% filter(`Read_Status` == "MODIFIED")

data_modified %>%
  select(., ,id, Aligned_Sequence, Reference_Sequence, Reference_Name, contains("read")) %>% 
  mutate(total = rowSums(across(8:22)), .after = `#Reads`) %>% 
  mutate(match = if_else(`#Reads`!=`total`,"mismatch", "ok" ), .after = `total`) %>% 
  filter(match == "mismatch")

# one row of mismatched data
data_modfied <- data_modified %>%
  select(., contains("read")) %>% 
  mutate(total = rowSums(across(4:18)), .after = `#Reads`) %>% 
  mutate(match = if_else(`#Reads`!=`total`,"mismatch", "ok" ), .after = `total`) %>% 
  filter(match != "mismatch")

data_modified %>%
  group_by(`Read_Status`) %>%
  summarise(sum = sum(`#Reads`))

data <- data_modified %>%
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



# Venn====

data_summary <- data %>% 
  summarise(
    sum_1 = sum(guide_1),
    sum_2 = sum(guide_2),
    sum_3 = sum(guide_3),
    sum_4 = sum(guide_4)
  )

## Calculate number of reads

# Precompute the new names based on your data_summary vector values:
new_names <- c(
  paste0("sgRNA 384 \n", data_summary$sum_1),
  paste0("sgRNA 362 \n", data_summary$sum_2),
  paste0("sgRNA 347\n",  data_summary$sum_3),
  paste0("sgRNA 338 \n", data_summary$sum_4)
)

# Process the data and apply the renaming with rename_at()
data_long <- data %>%
  select(`#Reads`:guide_4) %>% 
  uncount(weights = `#Reads`) %>% 
  rename_at(vars(guide_1, guide_2, guide_3, guide_4), ~ new_names)

indices_list <- lapply(data_long, function(x) which(x == 1))

ggVennDiagram(lapply(data_long, function(x) which(x == 1))) + 
  scale_fill_viridis(option = "cividis")+
  labs(fill="Read count") +
  theme_void(base_size = 10)

# ggVennDiagram(lapply(data_long, function(x) which(x == 1)), force_upset = T, sets.bar.color = c("lightblue", "blue", "blue", "blue"))

ggsave("figures/Venn.png", dpi = 900, width = 12, height = 8, units = "in") 

