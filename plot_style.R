# ==============================================================================
# plot_style.R
# Shared publication plotting style
# ==============================================================================
library(ggplot2)

fig_width_1col <- 89
fig_width_1.5col <- 120
fig_width_2col <- 183
fig_max_height <- 247

# Lineage palette
pal_lineage <- c(
  "PH908" = "#D55E00", 
  "R1a"   = "#0072B2"
)

# Branch-class palette
pal_class <- c(
  "Relic"   = "grey75",  
  "Middle"  = "grey45", 
  "Founder" = "black"   
)

pal_qc <- c(
  "accepted"            = "grey60", 
  "flagged_problem"     = "#E41A1C",
  "missing_coordinates" = "grey90",
  "other_reviewed"      = "grey30"
)

theme_scirep <- function(base_size = 7, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
      
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size + 1, hjust = 0),
      
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      
      plot.tag = element_text(face = "bold", size = 8),
      plot.tag.position = c(0, 1)
    )
}

save_scirep_fig <- function(filename, plot, width_mm, height_mm) {
  ggsave(filename = paste0(filename, ".pdf"), plot = plot, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(filename = paste0(filename, ".png"), plot = plot, width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
}