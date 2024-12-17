library(showtext)

font_add_google("Open Sans", "Sans")
plot.new()

theme_custom <- function(base_size=18, base_family="Sans"){
  theme_minimal(base_size = base_size, 
                base_family = base_family,
  ) %+replace%
    theme(axis.ticks = element_line(color = "grey92"),
          axis.ticks.length = unit(.5, "lines"),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white",
                                          colour="white"),
          legend.title = element_text(size=rel(1)),
          legend.text = element_text(color = "grey30"),
          plot.title = element_text(size = plot_title, face = "bold", hjust=0),
          plot.subtitle = element_text(size = base_size, color = "grey30"),
          plot.caption = element_text(size = rel(1.2), margin = margin(t = 15)),
          strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "cm"), hjust=0.5),
          strip.text=element_text(size=rel(1.5)))
}
