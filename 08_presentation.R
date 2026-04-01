# ==============================================================================
# 08_presentation.R
# Figure and table rendering layer
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(sf)
  library(rnaturalearth)
})

# ------------------------------------------------------------------------------
locate_config_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- NULL
  idx <- grep(paste0("^", file_arg), args)
  if (length(idx) > 0) {
    script_path <- normalizePath(sub(file_arg, "", args[idx[1]]), winslash = "/", mustWork = FALSE)
  }
  candidates <- unique(c(
    file.path(getwd(), "config_analysis.R"),
    if (!is.null(script_path)) file.path(dirname(script_path), "config_analysis.R") else NULL,
    if (!is.null(script_path)) file.path(dirname(script_path), "..", "config_analysis.R") else NULL
  ))
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit) || length(hit) == 0) {
    stop("Could not locate config_analysis.R. Place this script inside the project tree.", call. = FALSE)
  }
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

CONFIG_FILE <- locate_config_file()
source(CONFIG_FILE)
source(file.path(DIR_ROOT, "plot_style.R"))

setwd(DIR_ROOT)
dir.create(DIR_FIGURES, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# SUPPLEMENTARY FIGURE S1: Geocoding QC Map
# ==============================================================================
cat("Rendering Figure S1: Preprocessing Geocodes...\n")

df_geo <- readr::read_csv(FILE_FIGUREDATA_PREPROCESSING_GEOCODES, show_col_types = FALSE)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Identify flagged points for halo overlay
df_flagged <- df_geo %>% dplyr::filter(qc_group == "flagged_problem")

p_s1 <- ggplot() +
  geom_sf(data = world, fill = "#F5F5F5", color = "white", linewidth = 0.4) +
  coord_sf(
    xlim = c(BBOX_LON_MIN, BBOX_LON_MAX), 
    ylim = c(BBOX_LAT_MIN, BBOX_LAT_MAX), 
    expand = FALSE
  ) +
  # Base points
  geom_point(
    data = df_geo,
    aes(x = long, y = lat, color = qc_group, alpha = qc_group),
    size = 0.8, stroke = 0
  ) +
  # Halo overlay for flagged points
  geom_point(
    data = df_flagged,
    aes(x = long, y = lat),
    shape = 21, color = "#E41A1C", fill = NA, size = 3, stroke = 0.8
  ) +
  scale_color_manual(values = pal_qc, name = "QC Status", labels = c("Accepted", "Flagged/Excluded", "Missing Coords", "Other")) +
  scale_alpha_manual(values = c("accepted" = 0.4, "flagged_problem" = 1, "missing_coordinates" = 0, "other_reviewed" = 0.5), guide = "none") +
  theme_scirep() +
  theme(
    axis.title = element_blank(),
    legend.position = "top",
    legend.justification = "left",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    panel.background = element_rect(fill = "#E0F3F8")
  )

save_scirep_fig(file.path(DIR_FIGURES, "Figure_S1_Geocodes"), p_s1, width_mm = fig_width_1.5col, height_mm = 90)


# ==============================================================================
# MAIN FIGURE: Demographic Burden (Sci Rep Results Block 1)
# ==============================================================================
cat("Rendering Main Figure: Demographic Burden (Pillar 1)...\n")

df_main_p1 <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p1_main.csv"), show_col_types = FALSE)

labels_df <- df_main_p1 %>% 
  dplyr::filter(min_start_size == 1) %>% 
  dplyr::mutate(label = paste0(hg_group, "\n(n=", n_branches, ")")) %>%
  dplyr::select(hg_group, label)
hg_labels <- setNames(labels_df$label, labels_df$hg_group)

df_shares_long <- df_main_p1 %>%
  dplyr::filter(min_start_size == 1) %>%
  dplyr::select(hg_group, prop_relic, prop_middle, prop_founder) %>%
  tidyr::pivot_longer(cols = starts_with("prop_"), names_to = "branch_class", values_to = "mass_share") %>%
  dplyr::mutate(
    branch_class = stringr::str_to_title(stringr::str_replace(branch_class, "prop_", "")),
    branch_class = factor(branch_class, levels = c("Founder", "Middle", "Relic"))
  )

# Panel a: Stacked Bar
p_burden_a <- ggplot(df_shares_long, aes(x = hg_group, y = mass_share, fill = branch_class)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = pal_class, name = "Branch Class") +
  scale_x_discrete(labels = hg_labels) +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Lineage", y = "Proportion of terminal branches") +
  theme_scirep()

# Panel b: Feasibility Margins with Danger Zone
p_burden_b <- ggplot(df_main_p1, aes(x = min_start_size, y = worst_margin, color = hg_group, group = hg_group)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#E41A1C", alpha = 0.08) +
  annotate("text", x = 1.5, y = -15, label = "Structurally Infeasible", color = "#B22222", fontface = "italic", size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2, shape = 16) + 
  scale_x_continuous(breaks = c(1, 2, 3)) + 
  scale_color_manual(values = pal_lineage, name = "Lineage") +
  labs(
    x = "Minimum surviving lineages per branch",
    y = "Worst-case feasibility margin"
  ) +
  theme_scirep()

p_fig_burden <- p_burden_a + p_burden_b + 
  patchwork::plot_annotation(tag_levels = 'a') +
  patchwork::plot_layout(widths = c(1, 1.5))

save_scirep_fig(file.path(DIR_FIGURES, "Figure_Burden_Main"), p_fig_burden, width_mm = fig_width_2col, height_mm = 80)


# ==============================================================================
# SUPPLEMENTARY FIGURE S2: Feasibility Robustness Grid
# ==============================================================================
cat("Rendering Figure S2: Robustness Grid...\n")

df_grid <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p1_s2.csv"), show_col_types = FALSE)

p_s2 <- ggplot(df_grid, aes(x = as.factor(min_start_size), y = capacity_margin, color = hg_group)) +
  annotate("rect", xmin = 0, xmax = 4, ymin = -Inf, ymax = 0, fill = "#E41A1C", alpha = 0.08) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_jitter(alpha = 0.4, width = 0.15, size = 0.8, stroke = 0) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = NA, color = "black", linewidth = 0.4) +
  facet_wrap(~hg_group, scales = "free_y") +
  scale_color_manual(values = pal_lineage, guide = "none") +
  labs(
    x = "Minimum surviving lineages per branch",
    y = "Feasibility margin across demographic grid"
  ) +
  theme_scirep()

save_scirep_fig(file.path(DIR_FIGURES, "Figure_S2_Robustness_Grid"), p_s2, width_mm = fig_width_2col, height_mm = 90)


# ==============================================================================
# SUPPLEMENTARY FIGURE S3: Terminal Branch Size Distribution
# ==============================================================================
cat("Rendering Figure S3: Terminal Size Distribution...\n")

df_sizes <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p1_s3.csv"), show_col_types = FALSE)

df_sizes <- df_sizes %>%
  dplyr::arrange(hg_group, desc(branch_size)) %>%
  dplyr::group_by(hg_group) %>%
  dplyr::mutate(
    rank = dplyr::row_number(),
    cum_pct = cumsum(branch_size) / sum(branch_size)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(class_now = factor(class_now, levels = c("Relic", "Middle", "Founder")))

annot_df <- df_sizes %>%
  dplyr::filter(rank == 5) %>%
  dplyr::select(hg_group, cum_pct) %>%
  dplyr::mutate(label = sprintf("Top 5 = %.1f%% of mass", cum_pct * 100))

p_s3 <- ggplot(df_sizes, aes(x = rank, y = branch_size)) +
  geom_col(aes(fill = class_now), width = 1) +
  geom_text(data = annot_df, aes(x = Inf, y = Inf, label = label), 
            hjust = 1.1, vjust = 1.5, size = 2.5, fontface = "italic") +
  facet_wrap(~hg_group, scales = "free_x") +
  scale_fill_manual(values = pal_class, name = "Branch Class") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    x = "Terminal branches (ranked by size)",
    y = "Observed branch mass (cluster count)"
  ) +
  theme_scirep() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(1.5, "lines")
  )

save_scirep_fig(file.path(DIR_FIGURES, "Figure_S3_Terminal_Sizes"), p_s3, width_mm = fig_width_2col, height_mm = 90)

cat("Rendering complete! Check your /figures folder.\n")


# ==============================================================================
# MAIN FIGURE: Cross-Region Topology (Sci Rep Results Block 1, Part 2)
# ==============================================================================
cat("Rendering Main Figure: Cross-Region Topology (Pillar 2)...\n")

df_main_p2 <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p2_main.csv"), show_col_types = FALSE)

df_p2_long <- df_main_p2 %>%
  dplyr::select(rule_label, PH908 = cross_region_share_PH908, R1a = cross_region_share_R1a) %>%
  tidyr::pivot_longer(cols = c(PH908, R1a), names_to = "hg_group", values_to = "cross_region_share") %>%
  dplyr::mutate(
    rule_label = stringr::str_wrap(rule_label, width = 18),
    rule_label = factor(rule_label, levels = unique(rule_label)) 
  )

df_p2_annot <- df_main_p2 %>%
  dplyr::mutate(
    rule_label = stringr::str_wrap(rule_label, width = 18),
    rule_label = factor(rule_label, levels = unique(rule_label)),
    y_max = pmax(cross_region_share_PH908, cross_region_share_R1a),
    y_bracket = y_max + 0.03,
    y_text = y_bracket + 0.02,
    label = dplyr::case_when(
      p_value_one_sided_PH908_less < 0.001 ~ "p < 0.001",
      TRUE ~ sprintf("p = %.3f", p_value_one_sided_PH908_less)
    )
  )

p_topology <- ggplot(df_p2_long, aes(x = rule_label, y = cross_region_share, fill = hg_group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black", linewidth = 0.2) + 
  
  geom_segment(data = df_p2_annot, aes(x = as.numeric(rule_label) - 0.175, xend = as.numeric(rule_label) + 0.175, 
                                       y = y_bracket, yend = y_bracket), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = df_p2_annot, aes(x = as.numeric(rule_label) - 0.175, xend = as.numeric(rule_label) - 0.175, 
                                       y = y_bracket - 0.02, yend = y_bracket), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = df_p2_annot, aes(x = as.numeric(rule_label) + 0.175, xend = as.numeric(rule_label) + 0.175, 
                                       y = y_bracket - 0.02, yend = y_bracket), inherit.aes = FALSE, linewidth = 0.4) +
  
  geom_text(data = df_p2_annot, aes(x = rule_label, y = y_text, label = label), 
            inherit.aes = FALSE, size = 2.8, fontface = "italic") +
  
  scale_fill_manual(values = pal_lineage, name = "Lineage") +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Branch-link reconstruction rule",
    y = "Inferred cross-region linkage share"
  ) +
  theme_scirep() +
  theme(axis.text.x = element_text(vjust = 1))

save_scirep_fig(file.path(DIR_FIGURES, "Figure_Topology_Main"), p_topology, width_mm = fig_width_1.5col, height_mm = 80)


# ==============================================================================
# SUPPLEMENTARY FIGURE S4: Topology Permutation Nulls
# ==============================================================================
cat("Rendering Figure S4: Topology Null Distributions...\n")

df_null_p2 <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p2_null.csv"), show_col_types = FALSE)

df_null_primary <- df_null_p2 %>%
  dplyr::filter(null_model == "age_size_class") %>%
  dplyr::mutate(
    rule_label = stringr::str_wrap(rule_label, width = 25),
    rule_label = factor(rule_label, levels = unique(rule_label))
  )

df_obs_diff <- df_main_p2 %>%
  dplyr::mutate(
    rule_label = stringr::str_wrap(rule_label, width = 25),
    rule_label = factor(rule_label, levels = unique(rule_label))
  ) %>%
  dplyr::select(rule_label, cross_region_share_diff)

p_s4_null <- ggplot(df_null_primary, aes(x = cross_region_share_diff)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_histogram(fill = "grey75", color = "white", bins = 30, linewidth = 0.2) +
  geom_vline(data = df_obs_diff, aes(xintercept = cross_region_share_diff), 
             color = "#E41A1C", linewidth = 0.8) +
  
  facet_wrap(~rule_label, scales = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  labs(
    x = "Permuted difference in cross-region linkage (PH908 - R1a)",
    y = "Permutation count",
    subtitle = "Dashed line indicates zero difference; Red line indicates observed difference"
  ) +
  theme_scirep() +
  theme(plot.subtitle = element_text(face = "italic", color = "grey30", size = 6))

save_scirep_fig(file.path(DIR_FIGURES, "Figure_S4_Topology_Nulls"), p_s4_null, width_mm = fig_width_2col, height_mm = 70)

cat("Rendering complete for Pillar 2.\n")


# ==============================================================================
# MAIN FIGURE: Anchor-Packet Occupancy (Sci Rep Results Block 2, Part 1)
# ==============================================================================
cat("Rendering Main Figure: Anchor-Packet Occupancy (Pillar 3)...\n")

df_main_p3 <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p3_main.csv"), show_col_types = FALSE)
df_sens_p3 <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p3_sens.csv"), show_col_types = FALSE)
df_supp_p3 <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p3_supp.csv"), show_col_types = FALSE)

# ---------------------------------------------------------
# Panel A: Occupancy across Branch Panels
# ---------------------------------------------------------
df_p3_long <- df_main_p3 %>%
  dplyr::select(panel, PH908 = anchor_share_focal, R1a = anchor_share_control) %>%
  tidyr::pivot_longer(cols = c(PH908, R1a), names_to = "hg_group", values_to = "anchor_share") %>%
  dplyr::mutate(panel = stringr::str_to_title(panel),
                panel = factor(panel, levels = c("Relic", "Middle", "All", "Founder")),
                hg_group = factor(hg_group, levels = c("PH908", "R1a"))) # Force order

df_p3_annot <- df_main_p3 %>%
  dplyr::mutate(
    panel = stringr::str_to_title(panel),
    panel = factor(panel, levels = c("Relic", "Middle", "All", "Founder")),
    y_max = pmax(anchor_share_focal, anchor_share_control, na.rm = TRUE),
    y_bracket = y_max + 0.05,
    y_text = y_bracket + 0.04,
    label = dplyr::case_when(
      p_one_sided_anchor_share_diff < 0.001 ~ "p < 0.001",
      p_one_sided_anchor_share_diff < 0.05 ~ sprintf("p = %.3f", p_one_sided_anchor_share_diff),
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(label))

p_anchor_a <- ggplot(df_p3_long, aes(x = panel, y = anchor_share, fill = hg_group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8, color = "black", linewidth = 0.2) +
  
  geom_segment(data = df_p3_annot, aes(x = as.numeric(panel) - 0.2, xend = as.numeric(panel) + 0.2, 
                                       y = y_bracket, yend = y_bracket), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = df_p3_annot, aes(x = as.numeric(panel) - 0.2, xend = as.numeric(panel) - 0.2, 
                                       y = y_bracket - 0.02, yend = y_bracket), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = df_p3_annot, aes(x = as.numeric(panel) + 0.2, xend = as.numeric(panel) + 0.2, 
                                       y = y_bracket - 0.02, yend = y_bracket), inherit.aes = FALSE, linewidth = 0.4) +
  
  geom_text(data = df_p3_annot, aes(x = panel, y = y_text, label = label), 
            inherit.aes = FALSE, size = 2.8, fontface = "italic") +
  
  scale_fill_manual(values = pal_lineage, name = "Lineage") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.05), expand = expansion(mult = c(0, 0))) +
  labs(x = "Branch layer", y = "Retention inside primary anchor packet") +
  theme_scirep()

# ---------------------------------------------------------
# Panel B: Packet Composition (Horizontal Bar)
# ---------------------------------------------------------
df_comp <- df_supp_p3 %>%
  dplyr::filter(figure_component == "packet_composition") %>%
  tidyr::complete(dominant_region, hg_group, fill = list(branch_share = 0, branch_n = 0)) %>%
  dplyr::mutate(
    dominant_region = forcats::fct_reorder(dominant_region, branch_share, .fun = max),
    hg_group = factor(hg_group, levels = c("PH908", "R1a"))
  )

p_anchor_b <- ggplot(df_comp, aes(y = dominant_region, x = branch_share, fill = hg_group)) +
  geom_col(position = position_dodge2(width = 0.6, preserve = "single", padding = 0), 
           width = 0.45, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = pal_lineage, guide = "none") +
  scale_x_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Share of retained relic branches", y = "Anchor Region") +
  theme_scirep() +
  theme(
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 6) 
  )

# ---------------------------------------------------------
# Panel C: Sensitivity Across Thresholds
# ---------------------------------------------------------
df_sens_long <- df_sens_p3 %>%
  dplyr::select(anchor_id, PH908 = anchor_share_focal, R1a = anchor_share_control) %>%
  tidyr::pivot_longer(cols = c(PH908, R1a), names_to = "hg_group", values_to = "anchor_share") %>%
  dplyr::mutate(
    anchor_id = factor(anchor_id, levels = c("cum50", "cum60", "top4", "cum70", "top5")),
    hg_group = factor(hg_group, levels = c("PH908", "R1a"))
  )

p_anchor_c <- ggplot(df_sens_long, aes(x = anchor_id, y = anchor_share, color = hg_group, group = hg_group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5, shape = 21, fill = "white", stroke = 1) +
  scale_color_manual(values = pal_lineage, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Anchor definition threshold", y = "Relic retention share") +
  theme_scirep()

# ---------------------------------------------------------
# Combine P3 with Patchwork
# ---------------------------------------------------------
p_fig_anchor <- (p_anchor_a | p_anchor_c) / p_anchor_b + 
  patchwork::plot_annotation(tag_levels = 'a') +
  patchwork::plot_layout(heights = c(1.2, 1), guides = "collect") & theme(legend.position = "bottom")

save_scirep_fig(file.path(DIR_FIGURES, "Figure_Anchor_Main"), p_fig_anchor, width_mm = fig_width_2col, height_mm = 140)


# ==============================================================================
# SUPPLEMENTARY FIGURE S5: Relic Age Extension
# ==============================================================================
cat("Rendering Figure S5: Relic Age Extension...\n")

df_age_ext <- df_supp_p3 %>%
  dplyr::filter(figure_component == "relic_age_extension") %>%
  dplyr::select(age_mode_pretty, age_stratum_pretty, PH908 = anchor_share_focal, 
                R1a = anchor_share_control, p_one_sided_anchor_share_diff) %>%
  tidyr::pivot_longer(cols = c(PH908, R1a), names_to = "hg_group", values_to = "anchor_share")

df_age_ext$age_mode_pretty <- factor(df_age_ext$age_mode_pretty, 
                                     levels = c("All relics", "Old vs rest (median split)", "Within-lineage tertiles"))
df_age_ext$age_stratum_pretty <- factor(df_age_ext$age_stratum_pretty, 
                                        levels = c("All relics", "Old", "Rest", "Oldest tertile", "Middle tertile", "Youngest tertile"))
df_age_ext$hg_group <- factor(df_age_ext$hg_group, levels = c("PH908", "R1a"))

df_s5_annot <- df_age_ext %>%
  dplyr::group_by(age_mode_pretty, age_stratum_pretty) %>%
  dplyr::summarise(
    y_max = max(anchor_share, na.rm = TRUE),
    p_val = dplyr::first(p_one_sided_anchor_share_diff),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    y_text = y_max + 0.06,
    label = dplyr::case_when(
      p_val < 0.001 ~ "p < 0.001",
      p_val < 0.05 ~ sprintf("p = %.3f", p_val),
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(label))

p_s5_age <- ggplot(df_age_ext, aes(x = age_stratum_pretty, y = anchor_share, fill = hg_group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black", linewidth = 0.2) +
  
  geom_text(data = df_s5_annot, aes(x = age_stratum_pretty, y = y_text, label = label), 
            inherit.aes = FALSE, size = 2.8, fontface = "italic") +
  
  facet_grid(~age_mode_pretty, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal_lineage, name = "Lineage") +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Relic branch age stratum",
    y = "Retention inside primary anchor packet"
  ) +
  theme_scirep() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1) 
  )

save_scirep_fig(file.path(DIR_FIGURES, "Figure_S5_Relic_Age_Extension"), p_s5_age, width_mm = fig_width_2col, height_mm = 80)

# ==============================================================================
# MAIN FIGURE: Relic Concentration Contrast (Sci Rep Results Block 2, Part 2)
# ==============================================================================
cat("Rendering Main Figure: Relic Concentration (Pillar 4)...\n")

df_main_p4 <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p4_main.csv"), show_col_types = FALSE)
df_sens_p4 <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p4_sens.csv"), show_col_types = FALSE)

# ---------------------------------------------------------
# Panel A: Primary Contrast in Balanced Arena
# ---------------------------------------------------------
metric_levels_a <- c(
  "Shannon Entropy\n(Lower = more concentrated)", 
  "Effective # Regions\n(Lower = more concentrated)", 
  "Top-3 Region Share\n(Higher = more concentrated)"
)

df_main_p4 <- df_main_p4 %>%
  dplyr::mutate(
    metric_pretty = dplyr::case_when(
      metric == "shannon_entropy" ~ metric_levels_a[1],
      metric == "effective_regions" ~ metric_levels_a[2],
      metric == "top3_region_share" ~ metric_levels_a[3]
    ),
    metric_pretty = factor(metric_pretty, levels = metric_levels_a)
  )

df_p4_long <- df_main_p4 %>%
  dplyr::select(metric_pretty, PH908 = focal_value, R1a = control_value) %>%
  tidyr::pivot_longer(cols = c(PH908, R1a), names_to = "hg_group", values_to = "value") %>%
  dplyr::mutate(hg_group = factor(hg_group, levels = c("PH908", "R1a")))

df_p4_annot <- df_main_p4 %>%
  dplyr::mutate(
    y_max = pmax(focal_value, control_value),
    y_bracket = y_max + (y_max * 0.08),  
    y_text = y_bracket + (y_max * 0.04),
    label = dplyr::case_when(
      p_one_sided < 0.001 ~ "p < 0.001",
      p_one_sided < 0.05 ~ sprintf("p = %.3f", p_one_sided),
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(label))

p_concent_a <- ggplot(df_p4_long, aes(x = hg_group, y = value, fill = hg_group)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  
  geom_segment(data = df_p4_annot, aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = df_p4_annot, aes(x = 1, xend = 1, y = y_bracket - (y_max*0.02), yend = y_bracket), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = df_p4_annot, aes(x = 2, xend = 2, y = y_bracket - (y_max*0.02), yend = y_bracket), inherit.aes = FALSE, linewidth = 0.4) +
  
  geom_text(data = df_p4_annot, aes(x = 1.5, y = y_text, label = label), inherit.aes = FALSE, size = 2.8, fontface = "italic") +
  
  facet_wrap(~metric_pretty, scales = "free_y") +
  scale_fill_manual(values = pal_lineage, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = "Lineage", y = "Concentration metric value") +
  theme_scirep()

# ---------------------------------------------------------
# Panel B: Sensitivity Summary
# ---------------------------------------------------------
metric_levels_b <- c(
  "Difference: Shannon Entropy\n(< 0 means PH908 more conc.)", 
  "Difference: Effective # Regions\n(< 0 means PH908 more conc.)", 
  "Difference: Top-3 Region Share\n(> 0 means PH908 more conc.)"
)

df_sens_p4 <- df_sens_p4 %>%
  dplyr::mutate(
    metric_pretty = dplyr::case_when(
      metric == "shannon_entropy" ~ metric_levels_b[1],
      metric == "effective_regions" ~ metric_levels_b[2],
      metric == "top3_region_share" ~ metric_levels_b[3]
    ),
    metric_pretty = factor(metric_pretty, levels = metric_levels_b),
    panel_label = factor(panel_label, levels = c("PH908-R1a min1", "PH908-R1a min2", "I2-R1a min1", "I2-R1a min2")),
    fill_color = ifelse(direction_consistent, "PH908", "Other")
  )

p_concent_b <- ggplot(df_sens_p4, aes(x = panel_label, y = observed_diff, group = metric_pretty)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_line(color = "grey50", linewidth = 0.4) +
  geom_point(aes(fill = fill_color), shape = 21, size = 2.5, stroke = 0.7, color = "black") +
  facet_wrap(~metric_pretty, scales = "free_y") +
  scale_fill_manual(values = c("PH908" = pal_lineage[["PH908"]], "Other" = "white"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  labs(x = "Balanced shared-region arena variant", y = "Observed difference (Focal - Control)") +
  theme_scirep() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# ---------------------------------------------------------
# Combine P4 with Patchwork
# ---------------------------------------------------------
p_fig_concent <- p_concent_a / p_concent_b + 
  patchwork::plot_annotation(tag_levels = 'a') +
  patchwork::plot_layout(heights = c(1, 1))

save_scirep_fig(file.path(DIR_FIGURES, "Figure_Concentration_Main"), p_fig_concent, width_mm = fig_width_2col, height_mm = 130)

cat("Rendering complete for Pillar 4.\n")


# ==============================================================================
# MAIN FIGURE: Depth-Structured Coherence (Sci Rep Results Block 2, Part 3)
# ==============================================================================
cat("Rendering Main Figure: Depth-Structured Coherence (Pillar 5)...\n")

df_main_p5 <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p5_main.csv"), show_col_types = FALSE)
df_inset_p5 <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p5_inset.csv"), show_col_types = FALSE)

rule_levels <- c("Nearest older branch", "Nearest older within age window", "Distance + age cost")

# ---------------------------------------------------------
# Panel A: Depth-Profile Curves
# ---------------------------------------------------------
df_main_p5 <- df_main_p5 %>%
  dplyr::mutate(
    rule_label = stringr::str_wrap(rule_label, width = 22),
    rule_label = factor(rule_label, levels = stringr::str_wrap(rule_levels, width = 22)),
    hg_group = factor(hg_group, levels = c("PH908", "R1a"))
  )

p_depth_a <- ggplot(df_main_p5, aes(x = depth_index, y = same_region_rate, color = hg_group)) +
  geom_vline(xintercept = c(1, 2, 3, 4), color = "grey95", linewidth = 0.5) +
  geom_line(linewidth = 0.8, position = position_dodge(width = 0.1)) +
  geom_point(shape = 21, fill = "white", size = 2.5, stroke = 1, position = position_dodge(width = 0.1)) +
  facet_wrap(~rule_label) +
  scale_color_manual(values = pal_lineage, name = "Lineage") +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c("Tier 1\n(Shallow)", "Tier 2", "Tier 3", "Tier 4\n(Deep)")) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.05), expand = expansion(mult = c(0, 0))) +
  labs(
    x = "Inferred parent-age depth tier",
    y = "Same-region linkage rate"
  ) +
  theme_scirep() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(2, "lines")
  )

# ---------------------------------------------------------
# Panel B: Profile AUC Difference (Bar Chart Inset)
# ---------------------------------------------------------
df_inset_p5 <- df_inset_p5 %>%
  dplyr::mutate(
    rule_label = stringr::str_wrap(rule_label, width = 18),
    rule_label = factor(rule_label, levels = stringr::str_wrap(rule_levels, width = 18)),
    fill_color = ifelse(profile_auc_diff > 0, "PH908", "Other"),
    p_label = dplyr::case_when(
      profile_auc_diff_p_one_sided < 0.001 ~ "p < 0.001",
      profile_auc_diff_p_one_sided < 0.05 ~ sprintf("p = %.3f", profile_auc_diff_p_one_sided),
      TRUE ~ "n.s."
    ),
    y_nudge = ifelse(profile_auc_diff > 0, 0.02, -0.02),
    v_just = ifelse(profile_auc_diff > 0, 0, 1)
  )

p_depth_b <- ggplot(df_inset_p5, aes(x = rule_label, y = profile_auc_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  
  geom_col(aes(fill = fill_color), width = 0.5, color = "black", linewidth = 0.2) +
  geom_text(aes(y = profile_auc_diff + y_nudge, label = p_label, vjust = v_just), 
            size = 2.8, fontface = "italic") +
  
  scale_fill_manual(values = c("PH908" = pal_lineage[["PH908"]], "Other" = "grey70"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(
    x = "Branch-link reconstruction rule", 
    y = "Profile-wide Coherence Difference\n(AUC: > 0 means PH908 more coherent)"
  ) +
  theme_scirep()

# ---------------------------------------------------------
# Combine P5 with Patchwork
# ---------------------------------------------------------
p_fig_depth <- p_depth_a / p_depth_b + 
  patchwork::plot_annotation(tag_levels = 'a') +
  patchwork::plot_layout(heights = c(1.2, 1))

save_scirep_fig(file.path(DIR_FIGURES, "Figure_Depth_Main"), p_fig_depth, width_mm = fig_width_2col, height_mm = 140)


# ==============================================================================
# SUPPLEMENTARY FIGURE S6: Depth Coherence Permutation Nulls
# ==============================================================================
cat("Rendering Figure S6: Depth Coherence Nulls...\n")

df_null_p5 <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p5_null.csv"), show_col_types = FALSE)

df_null_p5 <- df_null_p5 %>%
  dplyr::mutate(
    rule_label = stringr::str_wrap(rule_label, width = 25),
    rule_label = factor(rule_label, levels = stringr::str_wrap(rule_levels, width = 25))
  )

# Extract observed values to draw the red lines safely
df_obs_diff_p5 <- df_null_p5 %>%
  dplyr::distinct(rule_label, observed_value)

p_s6_null <- ggplot(df_null_p5, aes(x = stat_value)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_histogram(fill = "grey75", color = "white", bins = 30, linewidth = 0.2) +
  geom_vline(data = df_obs_diff_p5, aes(xintercept = observed_value), 
             color = "#E41A1C", linewidth = 0.8) +
  facet_wrap(~rule_label, scales = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Permuted Profile AUC Difference (PH908 - R1a)",
    y = "Permutation count",
    subtitle = "Dashed line indicates zero difference; Red line indicates observed difference"
  ) +
  theme_scirep() +
  theme(plot.subtitle = element_text(face = "italic", color = "grey30", size = 6))

save_scirep_fig(file.path(DIR_FIGURES, "Figure_S6_Depth_Nulls"), p_s6_null, width_mm = fig_width_2col, height_mm = 70)

cat("Rendering complete for Pillar 5.\n")

# ==============================================================================
# MAIN FIGURE: Deep-Core Source-Sink Inversion (Sci Rep Results Block 3)
# ==============================================================================
cat("Rendering Main Figure: Source-Sink Inversion (Pillar 6)...\n")

df_carry <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p6_foundercarry.csv"), show_col_types = FALSE)
df_topk <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p6_topk.csv"), show_col_types = FALSE)
df_jack <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p6_jackknife.csv"), show_col_types = FALSE)
df_grid <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p6_threshold_grid.csv"), show_col_types = FALSE)

# ---------------------------------------------------------
# Panel A: Expected vs Observed Founders (The Slopegraph)
# ---------------------------------------------------------
df_carry_long <- df_carry %>%
  dplyr::select(lineage, Expected = expected_core_founders, Observed = observed_core_founders) %>%
  tidyr::pivot_longer(cols = c(Expected, Observed), names_to = "State", values_to = "Count") %>%
  dplyr::mutate(
    State = factor(State, levels = c("Expected", "Observed")),
    lineage = factor(lineage, levels = c("PH908", "R1a")),
    label_text = sprintf("%.1f", Count),
    label_text = stringr::str_replace(label_text, "\\.0$", ""),
    x_pos = ifelse(State == "Expected", 1 - 0.15, 2 + 0.15),
    h_align = ifelse(State == "Expected", 1, 0)
  )

p_sink_a <- ggplot(df_carry_long, aes(x = State, y = Count, color = lineage, group = lineage)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4, shape = 21, fill = "white", stroke = 1.5) +
  geom_text(aes(x = x_pos, label = label_text, hjust = h_align), 
            fontface = "bold", size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = pal_lineage, name = "Lineage") +
  scale_x_discrete(expand = expansion(add = c(0.4, 0.4))) + 
  scale_y_continuous(limits = c(0, max(df_carry_long$Count) + 1), expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Founder branches in Top-4 Core",
    y = "Count",
    title = "Founder Carry Inversion"
  ) +
  theme_scirep() +
  theme(
    panel.grid.major.y = element_blank(), 
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10))
  )

# ---------------------------------------------------------
# Panel B: Top-K Core Definitions (Forest Plot)
# ---------------------------------------------------------
df_topk <- df_topk %>%
  dplyr::mutate(
    core_label = paste0("Top ", topk),
    core_label = factor(core_label, levels = c("Top 5", "Top 4", "Top 3")), 
    p_label = dplyr::case_when(
      p_one_sided < 0.001 ~ "p < 0.001",
      p_one_sided < 0.05 ~ sprintf("p = %.3f", p_one_sided),
      TRUE ~ "n.s."
    ),
    text_x = pmax(anti_survivor_contrast, ci_upper, na.rm = TRUE) + 0.08
  )

p_sink_b <- ggplot(df_topk, aes(y = core_label, x = anti_survivor_contrast)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_linerange(aes(xmin = ci_lower, xmax = ci_upper), color = "black", linewidth = 0.6) +
  geom_point(size = 3.5, shape = 21, fill = pal_lineage[["PH908"]], stroke = 1, color = "black") +
  geom_text(aes(label = p_label, x = text_x), size = 2.8, fontface = "italic", hjust = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.35))) +
  labs(
    x = "Deep-Core Capture Skew Difference\n(> 0 means PH908 traps more relics)", 
    y = "Core Definition",
    title = "Core Sensitivity"
  ) +
  theme_scirep()

# ---------------------------------------------------------
# Panel C: Jackknife Leave-One-Out (Forest Plot)
# ---------------------------------------------------------
df_jack <- df_jack %>%
  dplyr::mutate(
    dropped_region = forcats::fct_reorder(dropped_region, observed),
    p_label = dplyr::case_when(
      p_one_sided < 0.001 ~ "p < 0.001",
      p_one_sided < 0.05 ~ sprintf("p = %.3f", p_one_sided),
      p_one_sided < 0.10 ~ sprintf("p = %.3f", p_one_sided),
      TRUE ~ "n.s."
    ),
    text_x = pmax(observed, ci_high, na.rm = TRUE) + 0.08
  )

p_sink_c <- ggplot(df_jack, aes(y = dropped_region, x = observed)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), color = "black", linewidth = 0.6) +
  geom_point(size = 3.5, shape = 21, fill = pal_lineage[["PH908"]], stroke = 1, color = "black") +
  geom_text(aes(label = p_label, x = text_x), size = 2.8, fontface = "italic", hjust = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.35))) +
  labs(
    x = "Deep-Core Capture Skew Difference", 
    y = "Region dropped from Top-4 Core",
    title = "Jackknife Robustness"
  ) +
  theme_scirep()

# ---------------------------------------------------------
# Panel D: Threshold Grid (Heatmap)
# ---------------------------------------------------------
df_grid <- df_grid %>%
  dplyr::filter(adequacy_flag == TRUE) %>%
  dplyr::mutate(
    p_label = dplyr::case_when(
      p_one_sided_top4 < 0.001 ~ "p < 0.001",
      TRUE ~ sprintf("p = %.3f", p_one_sided_top4)
    ),
    label_text = sprintf("%.3f\n(%s)", anti_survivor_contrast_top4, p_label)
  )

p_sink_d <- ggplot(df_grid, aes(x = as.factor(founder_min_n), y = as.factor(relic_max_n), fill = anti_survivor_contrast_top4)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = label_text), size = 2.8, fontface = "bold") +
  scale_fill_gradient(low = "#FFF5EB", high = pal_lineage[["PH908"]], limits = c(0, max(df_grid$anti_survivor_contrast_top4)), guide = "none") +
  labs(
    x = "Founder threshold (≥ clusters)", 
    y = "Relic threshold (≤ clusters)",
    title = "Threshold Stability"
  ) +
  theme_scirep() +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

# ---------------------------------------------------------
# Combine P6 with Patchwork
# ---------------------------------------------------------
p_fig_sink <- (p_sink_a | p_sink_b) / (p_sink_c | p_sink_d) + 
  patchwork::plot_annotation(tag_levels = 'a') +
  patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")

save_scirep_fig(file.path(DIR_FIGURES, "Figure_SourceSink_Main"), p_fig_sink, width_mm = fig_width_2col, height_mm = 160)


# ==============================================================================
# SUPPLEMENTARY FIGURE S7: Source-Sink Permutation Nulls
# ==============================================================================
cat("Rendering Figure S7: Source-Sink Nulls...\n")

df_null_p6 <- readr::read_csv(file.path(DIR_FIGUREDATA, "figuredata_p6_null.csv"), show_col_types = FALSE)

df_null_top4 <- df_null_p6 %>%
  dplyr::filter(topk == 4, statistic == "anti_survivor_contrast")

obs_val_top4 <- df_topk %>% dplyr::filter(topk == 4) %>% dplyr::pull(anti_survivor_contrast)

p_s7_null <- ggplot(df_null_top4, aes(x = stat_value)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_histogram(fill = "grey75", color = "white", bins = 30, linewidth = 0.2) +
  geom_vline(xintercept = obs_val_top4, color = "#E41A1C", linewidth = 0.8) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Permuted Deep-Core Capture Skew Difference (PH908 - R1a)",
    y = "Permutation count",
    subtitle = "Dashed line indicates zero difference; Red line indicates observed difference (Top-4 Core)"
  ) +
  theme_scirep() +
  theme(plot.subtitle = element_text(face = "italic", color = "grey30", size = 6))

save_scirep_fig(file.path(DIR_FIGURES, "Figure_S7_SourceSink_Nulls"), p_s7_null, width_mm = fig_width_1.5col, height_mm = 70)

cat("============================================================\n")
cat("ALL FIGURES RENDERED SUCCESSFULLY.\n")
cat("============================================================\n")


# ==============================================================================
# Supplementary table rendering
# ==============================================================================

cat("Rendering supplementary tables...\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(gt)
})

# ------------------------------------------------------------------------------
# Output folders
# ------------------------------------------------------------------------------
DIR_TABLE_PREVIEW_HTML <- file.path(DIR_ROOT, "supplement", "tables_preview_html")
DIR_TABLE_PREVIEW_CSV  <- file.path(DIR_ROOT, "supplement", "tables_preview_csv")
dir.create(DIR_TABLE_PREVIEW_HTML, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_TABLE_PREVIEW_CSV,  showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------
fmt_num <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "", format(round(x, digits), nsmall = digits, trim = TRUE))
}

fmt_num1 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "", format(round(x, 1), nsmall = 1, trim = TRUE))
}

fmt_int <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "", as.character(as.integer(round(x))))
}

fmt_pct01 <- function(x, digits = 1) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "", paste0(format(round(100 * x, digits), nsmall = digits, trim = TRUE), "%"))
}

fmt_p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  dplyr::case_when(
    is.na(x) ~ "",
    x < 0.001 ~ "<0.001",
    TRUE ~ format(round(x, 3), nsmall = 3, trim = TRUE)
  )
}

yesno <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "",
    x %in% c(TRUE, "TRUE", "true", 1, "1") ~ "Yes",
    x %in% c(FALSE, "FALSE", "false", 0, "0") ~ "No",
    TRUE ~ as.character(x)
  )
}

clean_label <- function(x) {
  x %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_squish() %>%
    stringr::str_to_sentence()
}

ensure_cols <- function(df, cols, fill = "") {
  for (nm in cols) {
    if (!nm %in% names(df)) df[[nm]] <- fill
  }
  df
}

write_preview_csv <- function(df, filename) {
  out <- file.path(DIR_TABLE_PREVIEW_CSV, filename)
  readr::write_excel_csv(df, out, na = "")
}

theme_gt_nature <- function(gt_obj) {
  gt_obj %>%
    gt::tab_options(
      table.width = gt::pct(100),
      heading.align = "left",
      table.font.size = gt::px(11),
      heading.title.font.size = gt::px(12),
      data_row.padding = gt::px(4),
      column_labels.font.weight = "bold",
      row_group.font.weight = "bold",
      source_notes.font.size = gt::px(10),
      footnotes.font.size = gt::px(10),
      table.border.top.width = gt::px(1.2),
      table.border.bottom.width = gt::px(1.2),
      column_labels.border.top.width = gt::px(1),
      column_labels.border.bottom.width = gt::px(1),
      row_group.border.top.width = gt::px(1),
      table_body.hlines.width = gt::px(0.4)
    ) %>%
    gt::opt_table_lines()
}

save_gt_html <- function(gt_obj, filename_base) {
  gt::gtsave(
    data = gt_obj,
    filename = file.path(DIR_TABLE_PREVIEW_HTML, paste0(filename_base, ".html"))
  )
}

save_table_bundle <- function(df, gt_obj, filename_base) {
  write_preview_csv(df, paste0(filename_base, ".csv"))
  save_gt_html(gt_obj, filename_base)
}

# ------------------------------------------------------------------------------
# Table file paths
# ------------------------------------------------------------------------------
fp_table_pre    <- file.path(DIR_TABLES, "table_preprocessing_supp.csv")
fp_table_p1_rx  <- file.path(DIR_TABLES, "table_p1_rx.csv")
fp_table_p2_rx  <- file.path(DIR_TABLES, "table_p2_rx.csv")
fp_table_p3_rx  <- file.path(DIR_TABLES, "table_p3_rx.csv")
fp_table_p6_threshold <- file.path(DIR_TABLES, "table_p6_threshold_sensitivity.csv")

# ------------------------------------------------------------------------------
# Figure-data file paths
# ------------------------------------------------------------------------------
fp_fd_p1_main <- file.path(DIR_FIGUREDATA, "figuredata_p1_main.csv")
fp_fd_p2_main <- file.path(DIR_FIGUREDATA, "figuredata_p2_main.csv")
fp_fd_p3_main <- file.path(DIR_FIGUREDATA, "figuredata_p3_main.csv")
fp_fd_p3_sens <- file.path(DIR_FIGUREDATA, "figuredata_p3_sens.csv")
fp_fd_p3_supp <- file.path(DIR_FIGUREDATA, "figuredata_p3_supp.csv")
fp_fd_p4_main <- file.path(DIR_FIGUREDATA, "figuredata_p4_main.csv")
fp_fd_p4_sens <- file.path(DIR_FIGUREDATA, "figuredata_p4_sens.csv")
fp_fd_p5_inset <- file.path(DIR_FIGUREDATA, "figuredata_p5_inset.csv")
fp_fd_p6_topk <- file.path(DIR_FIGUREDATA, "figuredata_p6_topk.csv")
fp_fd_p6_fc   <- file.path(DIR_FIGUREDATA, "figuredata_p6_foundercarry.csv")
fp_fd_p6_jk   <- file.path(DIR_FIGUREDATA, "figuredata_p6_jackknife.csv")

# ==============================================================================
# Table S1. Preprocessing summary metrics
# ==============================================================================
tab_s1 <- readr::read_csv(fp_table_pre, show_col_types = FALSE) %>%
  dplyr::rename(
    Metric = metric,
    Value = value,
    Note = note
  )

gt_s1 <- tab_s1 %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Table S1. Preprocessing summary metrics"
  ) %>%
  theme_gt_nature()

save_table_bundle(tab_s1, gt_s1, "Table_S1_preprocessing_summary_metrics")

# ==============================================================================
# Table S2. Burden architecture and feasibility summary
# built from figuredata_p1_main.csv
# ==============================================================================
raw_s2 <- readr::read_csv(fp_fd_p1_main, show_col_types = FALSE)

# ---- branch architecture block ----
s2_floor1 <- raw_s2 %>%
  dplyr::filter(min_start_size == 1) %>%
  dplyr::select(
    hg_group, n_branches, total_branch_mass, n_relic, n_middle, n_founder
  )

ph908_row <- s2_floor1 %>% dplyr::filter(hg_group == "PH908")
r1a_row   <- s2_floor1 %>% dplyr::filter(hg_group == "R1a")

tab_s2_arch <- tibble::tibble(
  Section = "Branch architecture",
  Metric  = c("Branches", "Total mass", "Relic branches", "Middle branches", "Founder branches"),
  PH908   = c(
    fmt_int(ph908_row$n_branches),
    fmt_int(ph908_row$total_branch_mass),
    fmt_int(ph908_row$n_relic),
    fmt_int(ph908_row$n_middle),
    fmt_int(ph908_row$n_founder)
  ),
  R1a     = c(
    fmt_int(r1a_row$n_branches),
    fmt_int(r1a_row$total_branch_mass),
    fmt_int(r1a_row$n_relic),
    fmt_int(r1a_row$n_middle),
    fmt_int(r1a_row$n_founder)
  )
)

# ---- feasibility block ----
tab_s2_feas <- raw_s2 %>%
  dplyr::transmute(
    Section = "Feasibility summary",
    Metric = paste0("Floor ", min_start_size),
    Lineage = hg_group,
    Value = paste0(
      "Worst margin ", fmt_num1(worst_margin),
      "; feasible models ", fmt_int(n_feasible_models), "/", fmt_int(n_total_models),
      "; grid feasible ", fmt_pct01(prop_grid_feasible, 1)
    )
  ) %>%
  tidyr::pivot_wider(names_from = Lineage, values_from = Value) %>%
  ensure_cols(c("PH908", "R1a")) %>%
  dplyr::select(Section, Metric, PH908, R1a)

tab_s2 <- dplyr::bind_rows(tab_s2_arch, tab_s2_feas)

gt_s2 <- tab_s2 %>%
  gt::gt(groupname_col = "Section") %>%
  gt::tab_header(
    title = "Table S2. Burden architecture and feasibility summary"
  ) %>%
  theme_gt_nature()

save_table_bundle(tab_s2, gt_s2, "Table_S2_burden_architecture_feasibility")

# ==============================================================================
# Table S3. Cross-region topology summary
# built from figuredata_p2_main.csv
# ==============================================================================
raw_s3 <- readr::read_csv(fp_fd_p2_main, show_col_types = FALSE)

tab_s3 <- raw_s3 %>%
  dplyr::transmute(
    Rule = rule_label,
    `PH908 edges` = fmt_int(n_edges_PH908),
    `R1a edges` = fmt_int(n_edges_R1a),
    `PH908 cross-region share` = fmt_num(cross_region_share_PH908, 3),
    `R1a cross-region share` = fmt_num(cross_region_share_R1a, 3),
    `PH908-R1a difference` = fmt_num(cross_region_share_diff, 3),
    `One-sided p-value` = fmt_p(p_value_one_sided_PH908_less),
    Support = support_label
  )

gt_s3 <- tab_s3 %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Table S3. Cross-region topology summary"
  ) %>%
  theme_gt_nature()

save_table_bundle(tab_s3, gt_s3, "Table_S3_cross_region_topology_summary")

# ==============================================================================
# Table S4. Anchor-packet retention summary
# curated from figuredata_p3_main / p3_sens / p3_supp
# ==============================================================================
raw_s4_main <- readr::read_csv(fp_fd_p3_main, show_col_types = FALSE)
raw_s4_sens <- readr::read_csv(fp_fd_p3_sens, show_col_types = FALSE)
raw_s4_supp <- readr::read_csv(fp_fd_p3_supp, show_col_types = FALSE)

# ---- primary anchor row ----
s4_primary_row <- raw_s4_main %>%
  dplyr::filter(panel == "relic", anchor_id == "cum60")

tab_s4_primary <- tibble::tibble(
  Section = "Primary anchor result",
  Metric = c(
    "Anchor regions",
    "PH908 relic branches in anchor",
    "R1a relic branches in anchor",
    "PH908 anchor share",
    "R1a anchor share",
    "Anchor-share difference",
    "Log2 anchor odds ratio",
    "One-sided p-value"
  ),
  PH908 = c(
    s4_primary_row$anchor_regions_concat,
    paste0(fmt_int(s4_primary_row$focal_anchor_n), "/", fmt_int(s4_primary_row$focal_n)),
    "",
    fmt_num(s4_primary_row$anchor_share_focal, 3),
    "",
    fmt_num(s4_primary_row$anchor_share_diff, 3),
    fmt_num(s4_primary_row$anchor_log2_odds_ratio, 3),
    fmt_p(s4_primary_row$p_one_sided_anchor_share_diff)
  ),
  R1a = c(
    "",
    "",
    paste0(fmt_int(s4_primary_row$control_anchor_n), "/", fmt_int(s4_primary_row$control_n)),
    "",
    fmt_num(s4_primary_row$anchor_share_control, 3),
    "",
    "",
    ""
  ),
  Difference = c("", "", "", "", "", fmt_num(s4_primary_row$anchor_share_diff, 3), "", ""),
  `One-sided p-value` = c("", "", "", "", "", "", "", fmt_p(s4_primary_row$p_one_sided_anchor_share_diff)),
  Support = c("", "", "", "", "", "", "", s4_primary_row$support_label)
)

# ---- neighboring anchor sensitivity ----
tab_s4_sens <- raw_s4_sens %>%
  dplyr::filter(panel == "relic", anchor_id %in% c("cum50", "cum60", "top4", "cum70", "top5")) %>%
  dplyr::transmute(
    Section = "Neighboring anchor sensitivity",
    Metric = paste0(anchor_id, " (", fmt_int(anchor_region_count), " regions)"),
    PH908 = fmt_num(anchor_share_focal, 3),
    R1a = fmt_num(anchor_share_control, 3),
    Difference = fmt_num(anchor_share_diff, 3),
    `One-sided p-value` = fmt_p(p_one_sided_anchor_share_diff),
    Support = support_label
  )

# ---- relic-age extension ----
tab_s4_age <- raw_s4_supp %>%
  dplyr::filter(figure_component == "relic_age_extension") %>%
  dplyr::transmute(
    Section = "Relic-age extension",
    Metric = age_stratum_pretty,
    PH908 = fmt_num(anchor_share_focal, 3),
    R1a = fmt_num(anchor_share_control, 3),
    Difference = fmt_num(anchor_share_diff, 3),
    `One-sided p-value` = fmt_p(p_one_sided_anchor_share_diff),
    Support = support_label
  )

tab_s4 <- dplyr::bind_rows(
  tab_s4_primary %>% dplyr::select(Section, Metric, PH908, R1a, Difference, `One-sided p-value`, Support),
  tab_s4_sens,
  tab_s4_age
)

gt_s4 <- tab_s4 %>%
  gt::gt(groupname_col = "Section") %>%
  gt::tab_header(
    title = "Table S4. Anchor-packet retention summary"
  ) %>%
  theme_gt_nature()

save_table_bundle(tab_s4, gt_s4, "Table_S4_anchor_packet_retention_summary")
# ==============================================================================
# Table S5. Balanced shared-region concentration summary
# curated from figuredata_p4_main / p4_sens
# ==============================================================================
raw_s5_main <- readr::read_csv(fp_fd_p4_main, show_col_types = FALSE)
raw_s5_sens <- readr::read_csv(fp_fd_p4_sens, show_col_types = FALSE)

tab_s5_primary <- raw_s5_main %>%
  dplyr::filter(panel_id == "ph908_r1a_min1") %>%
  dplyr::transmute(
    Section = "Primary balanced arena (PH908-R1a min1)",
    Metric = metric_label,
    PH908 = fmt_num(focal_value, 3),
    R1a = fmt_num(control_value, 3),
    Difference = fmt_num(observed_diff, 3),
    `One-sided p-value` = fmt_p(p_one_sided),
    Support = support_label
  )

tab_s5_sens <- raw_s5_sens %>%
  dplyr::transmute(
    Section = "Sensitivity panels",
    Metric = paste0(panel_label, " — ", metric_label),
    PH908 = ifelse(focal_group == "PH908", fmt_num(focal_value, 3), ""),
    R1a = ifelse(control_group == "R1a", fmt_num(control_value, 3), ""),
    Difference = fmt_num(observed_diff, 3),
    `One-sided p-value` = fmt_p(p_one_sided),
    Support = support_label
  )

tab_s5 <- dplyr::bind_rows(tab_s5_primary, tab_s5_sens)

gt_s5 <- tab_s5 %>%
  gt::gt(groupname_col = "Section") %>%
  gt::tab_header(
    title = "Table S5. Balanced shared-region concentration summary"
  ) %>%
  theme_gt_nature()

save_table_bundle(tab_s5, gt_s5, "Table_S5_balanced_shared_region_concentration")

# ==============================================================================
# Table S6. Depth-coherence primary summary statistics
# Table S7. Depth-coherence secondary diagnostics
# built from figuredata_p5_inset.csv
# ==============================================================================
raw_s6 <- readr::read_csv(fp_fd_p5_inset, show_col_types = FALSE)

# ---- Table S6: primary summary ----
tab_s6a <- raw_s6 %>%
  dplyr::transmute(
    Rule = rule_label,
    `PH908 edges` = fmt_int(n_edges_PH908),
    `R1a edges` = fmt_int(n_edges_R1a),
    `Shared tiers` = fmt_int(n_shared_tiers),
    `Shared tiers meeting minimum edges` = fmt_int(n_shared_tiers_meeting_min_edges),
    `Profile AUC diff` = fmt_num(profile_auc_diff, 3),
    `One-sided p (AUC)` = fmt_p(profile_auc_diff_p_one_sided),
    Support = support_label
  )

gt_s6a <- tab_s6a %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Table S6. Depth-coherence primary summary statistics",
    subtitle = "Rule-specific primary depth-coherence summary statistics for PH908 and R1a, including edge counts, shared-tier coverage, profile-wide AUC contrasts, and one-sided permutation support."
  ) %>%
  theme_gt_nature()

save_table_bundle(tab_s6a, gt_s6a, "Table_S6_depth_coherence_primary_summary_statistics")

# ---- Table S7: secondary diagnostics ----
tab_s6b <- raw_s6 %>%
  dplyr::transmute(
    Rule = rule_label,
    `Profile mean diff` = fmt_num(profile_mean_diff, 3),
    `One-sided p (mean)` = fmt_p(profile_mean_diff_p_one_sided),
    `Deep-weighted diff` = fmt_num(deep_weighted_diff, 3),
    `One-sided p (deep)` = fmt_p(deep_weighted_diff_p_one_sided),
    `Slope diff` = fmt_num(slope_diff, 3),
    `Deep minus shallow diff` = fmt_num(deep_minus_shallow_diff, 3),
    `All shared tiers PH908 > R1a` = yesno(all_shared_tiers_PH908_gt_R1a),
    `All tiers minimum-edge adequate` = yesno(all_shared_tiers_min_edges_ok)
  )

gt_s6b <- tab_s6b %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Table S7. Depth-coherence secondary diagnostics",
    subtitle = "Secondary depth-coherence summary statistics for PH908 and R1a, including profile-mean contrasts, deep-weighted contrasts, slope summaries, and tier-adequacy indicators."
  ) %>%
  theme_gt_nature()

save_table_bundle(tab_s6b, gt_s6b, "Table_S7_depth_coherence_secondary_diagnostics")

# ==============================================================================
# Table S8. Deep-core source-sink summary
# ==============================================================================
raw_s7_topk <- readr::read_csv(fp_fd_p6_topk, show_col_types = FALSE)
raw_s7_fc   <- readr::read_csv(fp_fd_p6_fc,   show_col_types = FALSE)
raw_s7_jk   <- readr::read_csv(fp_fd_p6_jk,   show_col_types = FALSE)

wrap_core_regions <- function(x) {
  x %>%
    stringr::str_replace_all(";\\s*", ";\n")
}

# ---- primary and neighboring core definitions ----
tab_s7_core <- raw_s7_topk %>%
  dplyr::transmute(
    Section = "Deep-core definitions and primary contrasts",
    Metric = paste0("Top", topk, ifelse(is_primary, " (primary)", "")),
    `Core regions` = wrap_core_regions(core_regions),
    `PH908 relic / founder` = paste0(fmt_num(core_relic_share_PH908, 3), " / ", fmt_num(core_founder_share_PH908, 3)),
    `R1a relic / founder` = paste0(fmt_num(core_relic_share_R1a, 3), " / ", fmt_num(core_founder_share_R1a, 3)),
    `Anti-survivor contrast` = fmt_num(anti_survivor_contrast, 3),
    `One-sided p-value` = fmt_p(p_one_sided),
    Support = support_label
  )

# ---- founder carry translation ----
tab_s7_fc <- raw_s7_fc %>%
  dplyr::transmute(
    Section = "Founder-carry translation",
    Metric = lineage,
    `Core regions` = "",
    `PH908 relic / founder` = "",
    `R1a relic / founder` = "",
    `Anti-survivor contrast` = paste0(
      fmt_int(observed_core_founders), " / ", fmt_num(expected_core_founders, 3)
    ),
    `One-sided p-value` = "",
    Support = paste0("Observed / expected = ", fmt_num(observed_over_expected, 3))
  )

# ---- leave-one-out jackknife ----
tab_s7_jk <- raw_s7_jk %>%
  dplyr::transmute(
    Section = "Leave-one-out jackknife",
    Metric = dropped_region,
    `Core regions` = "",
    `PH908 relic / founder` = "",
    `R1a relic / founder` = "",
    `Anti-survivor contrast` = fmt_num(observed, 3),
    `One-sided p-value` = fmt_p(p_one_sided),
    Support = interpretation
  )

tab_s7 <- dplyr::bind_rows(tab_s7_core, tab_s7_fc, tab_s7_jk)

gt_s7 <- tab_s7 %>%
  gt::gt(groupname_col = "Section") %>%
  gt::tab_header(
    title = "Table S8. Deep-core source-sink summary",
    subtitle = "Summary metrics for the PH908-defined deep-core source-sink analysis, including neighboring core contrasts, founder-carry translation, and leave-one-out jackknife results."
  ) %>%
  gt::cols_label(
    Metric = "Metric",
    `Core regions` = "Core regions",
    `PH908 relic / founder` = "PH908\nrelic / founder",
    `R1a relic / founder` = "R1a\nrelic / founder",
    `Anti-survivor contrast` = "Anti-survivor\ncontrast",
    `One-sided p-value` = "One-sided\np-value",
    Support = "Support"
  ) %>%
  gt::tab_style(
    style = gt::cell_text(whitespace = "pre-line"),
    locations = gt::cells_body(columns = c(`Core regions`))
  ) %>%
  theme_gt_nature()

save_table_bundle(tab_s7, gt_s7, "Table_S8_deep_core_source_sink_summary")

# ==============================================================================
# Table S9. Threshold sensitivity grid for the top4 deep-core analysis
# ==============================================================================
tab_sx3 <- readr::read_csv(fp_table_p6_threshold, show_col_types = FALSE) %>%
  dplyr::transmute(
    `Relic threshold (<=)` = fmt_int(relic_max_n),
    `Founder threshold (>=)` = fmt_int(founder_min_n),
    `Anti-survivor contrast (top4)` = fmt_num(anti_survivor_contrast_top4, 3),
    `One-sided p (top4)` = fmt_p(p_one_sided_top4),
    `Core minus expanse (top4)` = fmt_num(core_minus_expanse_top4, 3),
    `One-sided p (core minus expanse)` = fmt_p(p_one_sided_core_minus_expanse_top4),
    `Adequacy flag` = yesno(adequacy_flag),
    Support = support_label
  )

gt_sx3 <- tab_sx3 %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Table S9. Threshold sensitivity grid for the top4 deep-core analysis"
  ) %>%
  theme_gt_nature()

save_table_bundle(tab_sx3, gt_sx3, "Table_S9_threshold_sensitivity_grid_top4")

# ==============================================================================
# Table S10. Leave-one-region-out robustness for the anchor analysis
# ==============================================================================
tab_rx1 <- readr::read_csv(fp_table_p3_rx, show_col_types = FALSE) %>%
  dplyr::transmute(
    `Dropped region` = dropped_region,
    `Anchor ID` = anchor_id,
    `Original anchor region count` = fmt_int(original_anchor_region_count),
    `Reduced anchor region count` = fmt_int(reduced_anchor_region_count),
    `PH908 anchor share` = fmt_num(anchor_share_focal, 3),
    `R1a anchor share` = fmt_num(anchor_share_control, 3),
    `Anchor-share difference` = fmt_num(anchor_share_diff, 3),
    `Log2 anchor odds ratio` = fmt_num(anchor_log2_odds_ratio, 3),
    `Baseline anchor-share difference` = fmt_num(baseline_anchor_share_diff, 3),
    `Attenuation vs baseline` = fmt_num(attenuation_vs_baseline, 3),
    `Sign positive` = yesno(sign_positive)
  )

gt_rx1 <- tab_rx1 %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Table S10. Leave-one-region-out robustness for the anchor analysis"
  ) %>%
  theme_gt_nature()

save_table_bundle(tab_rx1, gt_rx1, "Table_S10_leave_one_region_out_anchor")

# ==============================================================================
# Table S11. Influential-unit robustness for burden feasibility
# ==============================================================================
tab_rx2a <- readr::read_csv(fp_table_p1_rx, show_col_types = FALSE) %>%
  dplyr::transmute(
    Perturbation = clean_label(perturbation),
    Lineage = hg_group,
    `Minimum floor` = fmt_int(min_start_size),
    `Worst margin` = fmt_num1(worst_margin),
    `Feasible models` = paste0(fmt_int(n_feasible_models), "/", fmt_int(n_total_models)),
    `Branches after perturbation` = fmt_int(n_branches_after),
    `Total mass after perturbation` = fmt_int(total_mass_after),
    `Baseline worst margin` = fmt_num1(baseline_worst_margin),
    `Sign flip vs baseline` = yesno(sign_flip_vs_baseline)
  ) %>%
  dplyr::arrange(Lineage, Perturbation, `Minimum floor`)

gt_rx2a <- tab_rx2a %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Table S11. Influential-unit robustness for burden feasibility"
  ) %>%
  theme_gt_nature()

save_table_bundle(tab_rx2a, gt_rx2a, "Table_S11_influential_unit_burden")

# ==============================================================================
# Table S12. Influential-unit robustness for cross-region topology
# ==============================================================================
tab_rx2b <- readr::read_csv(fp_table_p2_rx, show_col_types = FALSE) %>%
  dplyr::transmute(
    Perturbation = clean_label(perturbation),
    Rule = rule,
    `Cross-region share difference` = fmt_num(cross_region_share_diff, 3),
    `Negative sign retained` = yesno(sign_negative),
    `Baseline negative sign` = yesno(baseline_sign_negative),
    `Sign matches baseline` = yesno(sign_matches_baseline)
  ) %>%
  dplyr::arrange(Rule, Perturbation)

gt_rx2b <- tab_rx2b %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Table S12. Influential-unit robustness for cross-region topology"
  ) %>%
  theme_gt_nature()

save_table_bundle(tab_rx2b, gt_rx2b, "Table_S12_influential_unit_topology")

cat("Done. Final Nature-style HTML tables written to:\n")
cat(DIR_TABLE_PREVIEW_HTML, "\n")
cat("Curated Excel-friendly CSV copies written to:\n")
cat(DIR_TABLE_PREVIEW_CSV, "\n")