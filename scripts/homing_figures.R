
source("scripts/cutting_analysis.R")


#===================================

## Figures ====

homing_summary <- emmeans::emmeans(homing_model, specs= ~ gRNA_type*pre_cut, type="response") %>% as_tibble() 

count_data <- homing_data %>% 
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
`g338-384` <- "#A1AA11"
`g384` <- "#91A6D2"
`g225` <- "#BF9900"


colours_fig_1b <- c(`g338-384`, `g384_del`)
colours_fig_2c <- c(`g338-384`, `g384`)
colours_fig_3b <- c(`g384`)
colours_fig_3d <- c(`g225`)

fig_1b <- homing_data %>% 
  filter(gRNA_type == "cd<sup><i>g384_del</i></sup>" | gRNA_type == "cd<sup><i>g338-384</i></sup>") %>% 
  filter(pre_cut == "wild-type")

fig_1b_mean <- homing_summary %>% 
  filter(gRNA_type == "cd<sup><i>g384_del</i></sup>" | gRNA_type == "cd<sup><i>g338-384</i></sup>") %>% 
  filter(pre_cut == "wild-type")

fig_2c <- homing_data %>% 
  filter(gRNA_type == "cd<sup><i>g384</i></sup>" | gRNA_type == "cd<sup><i>g338-384</i></sup>") %>% 
  filter(pre_cut == "cd<sup><i>384</i></sup>")

fig_2c_mean <- homing_summary %>% 
  filter(gRNA_type == "cd<sup><i>g384</i></sup>" | gRNA_type == "cd<sup><i>g338-384</i></sup>") %>% 
  filter(pre_cut == "cd<sup><i>384</i></sup>")

fig_3b <- homing_data %>% 
  filter(gRNA_type == "cd<sup><i>g384</i></sup>") %>% 
  filter(pre_cut == "cd<sup><i>225</i></sup>" | pre_cut == "wild-type")

fig_3b_mean <- homing_summary %>% 
  filter(gRNA_type == "cd<sup><i>g384</i></sup>") %>% 
  filter(pre_cut == "cd<sup><i>225</i></sup>" | pre_cut == "wild-type")

fig_3d <- homing_data %>% 
  filter(gRNA_type == "cd<sup><i>g225</i></sup>") %>% 
  filter(pre_cut == "cd<sup><i>384</i></sup>" | pre_cut == "wild-type")

fig_3d_mean <- homing_summary %>% 
  filter(gRNA_type == "cd<sup><i>g225</i></sup>") %>% 
  filter(pre_cut == "cd<sup><i>384</i></sup>" | pre_cut == "wild-type")

pwalk(
  tibble(
    data = list(fig_1b, fig_2c, fig_3b, fig_3d),
    fig_mean = list(fig_1b_mean, fig_2c_mean, fig_3b_mean, fig_3d_mean),
    colours = list(colours_fig_1b, colours_fig_2c, colours_fig_3b, colours_fig_3d),
    figures = list("fig 1b", "fig 2c", "fig 3b", "fig 3d")
  ),
  ~ {
   p <-  homing_plots(..1,..2,..3)
   ggsave(glue::glue("figures/{..4}_homing_plot.png"),
          p,
          dpi = 900, width = 9, height = 6, units = "in")
  }
)

