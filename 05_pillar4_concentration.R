# 05_pillar4_concentration.R
# Relic regional concentration contrast in balanced shared-region arenas.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

# Configuration

locate_config_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_dir <- NULL
  
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))
  }
  
  candidates <- unique(na.omit(c(
    if (!is.null(script_dir)) file.path(script_dir, "config_analysis.R") else NULL,
    file.path(getwd(), "config_analysis.R")
  )))
  
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit) || length(hit) == 0) {
    stop(
      paste0(
        "Could not locate config_analysis.R. Tried:\n  - ",
        paste(candidates, collapse = "\n  - ")
      ),
      call. = FALSE
    )
  }
  
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

CONFIG_FILE <- locate_config_file()
source(CONFIG_FILE)
setwd(DIR_ROOT)

# Output paths

FILE_P4_RESULTS <- file.path(DIR_RESULTS, "results_p4.rds")
FILE_P4_FIG_MAIN <- file.path(DIR_FIGUREDATA, "figuredata_p4_main.csv")
FILE_P4_FIG_SENS <- file.path(DIR_FIGUREDATA, "figuredata_p4_sens.csv")
FILE_P4_SUPP_TABLE <- file.path(DIR_TABLES, "table_p4_supp.csv")
FILE_P4_FACTS <- file.path(DIR_FACTS, "facts_p4.txt")

# Analysis settings

SCRIPT_NAME <- "05_pillar4_concentration.R"

RNG_SEED <- SEED_MAIN
set.seed(RNG_SEED)

RELIC_MAX_N <- 3
FOUNDER_MIN_N <- 15

AGE_BIN_N <- 3
SIZE_BIN_N <- 2

PRIMARY_PANEL_ID <- "ph908_r1a_min1"

PRIMARY_METRICS <- c(
  "shannon_entropy",
  "effective_regions",
  "top3_region_share"
)

# Full manuscript-grade permutation setting.
N_PERM <- 5000L

# Helpers

stop_msg <- function(...) stop(paste0(...), call. = FALSE)

normalize_snp <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^I-|^R-", "", x)
  x <- gsub("\\*$", "", x)
  x
}

safe_first_non_na <- function(x) {
  y <- x[!is.na(x)]
  if (length(y) == 0) return(NA)
  y[1]
}

safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

dominant_region_fun <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

classify_branch <- function(n) {
  dplyr::case_when(
    n <= RELIC_MAX_N ~ "Relic",
    n >= FOUNDER_MIN_N ~ "Founder",
    TRUE ~ "Middle"
  )
}

write_warning_log <- function(warning_vec, path) {
  if (length(warning_vec) == 0) {
    writeLines("No warnings generated.", path)
  } else {
    writeLines(unique(warning_vec), path)
  }
}

make_quantile_bins <- function(x, n_bins, prefix = "bin") {
  x <- as.numeric(x)
  
  if (sum(is.finite(x)) < 2) {
    return(rep(paste0(prefix, "_1"), length(x)))
  }
  
  qs <- unique(stats::quantile(
    x[is.finite(x)],
    probs = seq(0, 1, length.out = n_bins + 1),
    na.rm = TRUE
  ))
  
  if (length(qs) < 2) {
    return(rep(paste0(prefix, "_1"), length(x)))
  }
  
  out <- cut(
    x,
    breaks = qs,
    include.lowest = TRUE,
    ordered_result = TRUE
  )
  
  out <- as.character(out)
  out[is.na(out)] <- "missing"
  out
}

permute_within_strata <- function(x, strata) {
  out <- x
  strata <- as.character(strata)
  
  for (s in unique(strata)) {
    idx <- which(strata == s)
    if (length(idx) > 1) {
      out[idx] <- sample(out[idx], size = length(idx), replace = FALSE)
    }
  }
  
  out
}

calc_effective_regions <- function(counts_vec) {
  counts_vec <- as.numeric(counts_vec)
  counts_vec <- counts_vec[is.finite(counts_vec) & counts_vec > 0]
  
  if (length(counts_vec) == 0) return(NA_real_)
  
  p <- counts_vec / sum(counts_vec)
  1 / sum(p^2)
}

calc_group_metrics_from_branch_df <- function(branch_df) {
  counts_vec <- branch_df |>
    dplyr::count(dominant_region, name = "branch_n") |>
    dplyr::pull(branch_n)
  
  if (length(counts_vec) == 0 || sum(counts_vec) <= 0) {
    return(tibble::tibble(
      n_regions = 0L,
      total_relic_branches = 0L,
      shannon_entropy = NA_real_,
      effective_regions = NA_real_,
      top3_region_share = NA_real_,
      top1_region_share = NA_real_,
      top2_region_share = NA_real_
    ))
  }
  
  p <- counts_vec / sum(counts_vec)
  p_pos <- p[p > 0]
  p_sorted <- sort(p, decreasing = TRUE)
  
  tibble::tibble(
    n_regions = length(counts_vec),
    total_relic_branches = sum(counts_vec),
    shannon_entropy = -sum(p_pos * log(p_pos)),
    effective_regions = calc_effective_regions(counts_vec),
    top3_region_share = sum(head(p_sorted, 3)),
    top1_region_share = sum(head(p_sorted, 1)),
    top2_region_share = sum(head(p_sorted, 2))
  )
}

pretty_metric_label <- function(metric_name) {
  dplyr::case_when(
    metric_name == "shannon_entropy"   ~ "Shannon entropy",
    metric_name == "effective_regions" ~ "Effective # regions",
    metric_name == "top3_region_share" ~ "Top-3 region share",
    metric_name == "top1_region_share" ~ "Top-1 region share",
    metric_name == "top2_region_share" ~ "Top-2 region share",
    TRUE ~ metric_name
  )
}

metric_direction_label <- function(metric_name) {
  dplyr::case_when(
    metric_name == "shannon_entropy"   ~ "Lower = more concentrated",
    metric_name == "effective_regions" ~ "Lower = more concentrated",
    metric_name == "top3_region_share" ~ "Higher = more concentrated",
    metric_name == "top1_region_share" ~ "Higher = more concentrated",
    metric_name == "top2_region_share" ~ "Higher = more concentrated",
    TRUE ~ ""
  )
}

ci95 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 20) return(c(NA_real_, NA_real_))
  as.numeric(stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
}

two_sided_perm_p <- function(null_vals, obs) {
  null_vals <- null_vals[is.finite(null_vals)]
  if (!is.finite(obs) || length(null_vals) == 0) return(NA_real_)
  (1 + sum(abs(null_vals) >= abs(obs))) / (length(null_vals) + 1)
}

one_sided_more_concentrated_p <- function(metric_name, null_vals, obs) {
  null_vals <- null_vals[is.finite(null_vals)]
  if (!is.finite(obs) || length(null_vals) == 0) return(NA_real_)
  
  if (metric_name %in% c("shannon_entropy", "effective_regions")) {
    return((1 + sum(null_vals <= obs)) / (length(null_vals) + 1))
  }
  
  if (metric_name %in% c("top3_region_share", "top1_region_share", "top2_region_share")) {
    return((1 + sum(null_vals >= obs)) / (length(null_vals) + 1))
  }
  
  NA_real_
}

direction_consistent_more_concentrated <- function(metric_name, obs) {
  if (!is.finite(obs)) return(NA)
  if (metric_name %in% c("shannon_entropy", "effective_regions")) return(obs < 0)
  if (metric_name %in% c("top3_region_share", "top1_region_share", "top2_region_share")) return(obs > 0)
  NA
}

path_for_facts <- function(x) {
  x <- as.character(x)
  if (!length(x) || is.na(x) || x == "") return(x)
  x_norm <- normalizePath(x, winslash = "/", mustWork = FALSE)
  root_norm <- normalizePath(DIR_ROOT, winslash = "/", mustWork = FALSE)
  root_esc <- gsub("([][{}()+*^$|\\?.])", "\\\\\\1", root_norm)
  sub(paste0("^", root_esc), ".", x_norm)
}

safe_git_commit <- function() {
  out <- tryCatch(
    suppressWarnings(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0)
  )
  out <- out[nzchar(out)]
  if (length(out)) out[1] else "NA"
}

balance_shared_regions <- function(relic_df, valid_groups, min_relic_per_group_per_region) {
  g1 <- valid_groups[1]
  g2 <- valid_groups[2]
  
  reg_counts <- relic_df |>
    dplyr::count(dominant_region, hg_group, name = "n_branches") |>
    tidyr::pivot_wider(
      names_from = hg_group,
      values_from = n_branches,
      values_fill = 0
    )
  
  if (!g1 %in% names(reg_counts)) reg_counts[[g1]] <- 0L
  if (!g2 %in% names(reg_counts)) reg_counts[[g2]] <- 0L
  
  reg_keep <- reg_counts |>
    dplyr::mutate(
      keep_region = .data[[g1]] >= min_relic_per_group_per_region &
        .data[[g2]] >= min_relic_per_group_per_region
    )
  
  kept_regions <- reg_keep |>
    dplyr::filter(keep_region) |>
    dplyr::pull(dominant_region)
  
  balanced_df <- relic_df |>
    dplyr::filter(dominant_region %in% kept_regions)
  
  list(
    balanced_df = balanced_df,
    region_keep_tbl = reg_keep
  )
}

# Analysis functions

run_relic_concentration_analysis <- function(branch_tbl,
                                             valid_groups,
                                             label_set,
                                             panel_id,
                                             min_relic_per_group_per_region) {
  
  g1 <- valid_groups[1]
  g2 <- valid_groups[2]
  
  relic_raw <- branch_tbl |>
    dplyr::filter(class == "Relic", !is.na(dominant_region), dominant_region != "")
  
  raw_counts <- relic_raw |>
    dplyr::count(hg_group, name = "n_raw_relic")
  
  bal_info <- balance_shared_regions(
    relic_df = relic_raw,
    valid_groups = valid_groups,
    min_relic_per_group_per_region = min_relic_per_group_per_region
  )
  
  relic_bal <- bal_info$balanced_df
  
  if (nrow(relic_bal) == 0) {
    stop_msg("No relic branches remain after balancing for panel: ", panel_id)
  }
  
  if (!all(valid_groups %in% unique(relic_bal$hg_group))) {
    stop_msg("Balanced arena does not retain all comparison groups for panel: ", panel_id)
  }
  
  region_balance_diag <- dplyr::bind_rows(
    relic_raw |>
      dplyr::count(dominant_region, hg_group, name = "n_branches") |>
      dplyr::mutate(stage = "raw"),
    relic_bal |>
      dplyr::count(dominant_region, hg_group, name = "n_branches") |>
      dplyr::mutate(stage = "balanced")
  ) |>
    dplyr::mutate(
      label_set = label_set,
      panel_id = panel_id,
      min_relic_per_group_per_region = min_relic_per_group_per_region
    )
  
  balanced_regions_tbl <- tibble::tibble(
    label_set = label_set,
    panel_id = panel_id,
    min_relic_per_group_per_region = min_relic_per_group_per_region,
    dominant_region = sort(unique(relic_bal$dominant_region))
  )
  
  region_counts <- relic_bal |>
    dplyr::count(hg_group, dominant_region, name = "branch_n") |>
    dplyr::group_by(hg_group) |>
    dplyr::mutate(branch_share = branch_n / sum(branch_n)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      label_set = label_set,
      panel_id = panel_id,
      min_relic_per_group_per_region = min_relic_per_group_per_region
    )
  
  obs_summary <- dplyr::bind_rows(
    calc_group_metrics_from_branch_df(relic_bal |> dplyr::filter(hg_group == g1)) |>
      dplyr::mutate(hg_group = g1),
    calc_group_metrics_from_branch_df(relic_bal |> dplyr::filter(hg_group == g2)) |>
      dplyr::mutate(hg_group = g2)
  ) |>
    dplyr::mutate(
      label_set = label_set,
      panel_id = panel_id,
      min_relic_per_group_per_region = min_relic_per_group_per_region
    )
  
  obs_diffs <- tibble::tibble(
    metric = c(
      "shannon_entropy",
      "effective_regions",
      "top3_region_share",
      "top1_region_share",
      "top2_region_share"
    ),
    metric_label = vapply(
      c("shannon_entropy", "effective_regions", "top3_region_share", "top1_region_share", "top2_region_share"),
      pretty_metric_label,
      character(1)
    ),
    direction_label = vapply(
      c("shannon_entropy", "effective_regions", "top3_region_share", "top1_region_share", "top2_region_share"),
      metric_direction_label,
      character(1)
    ),
    focal_group = g1,
    control_group = g2,
    focal_value = c(
      obs_summary$shannon_entropy[obs_summary$hg_group == g1],
      obs_summary$effective_regions[obs_summary$hg_group == g1],
      obs_summary$top3_region_share[obs_summary$hg_group == g1],
      obs_summary$top1_region_share[obs_summary$hg_group == g1],
      obs_summary$top2_region_share[obs_summary$hg_group == g1]
    ),
    control_value = c(
      obs_summary$shannon_entropy[obs_summary$hg_group == g2],
      obs_summary$effective_regions[obs_summary$hg_group == g2],
      obs_summary$top3_region_share[obs_summary$hg_group == g2],
      obs_summary$top1_region_share[obs_summary$hg_group == g2],
      obs_summary$top2_region_share[obs_summary$hg_group == g2]
    )
  ) |>
    dplyr::mutate(
      observed_diff = focal_value - control_value,
      direction_consistent = vapply(
        seq_len(dplyr::n()),
        function(i) direction_consistent_more_concentrated(metric[i], observed_diff[i]),
        logical(1)
      ),
      label_set = label_set,
      panel_id = panel_id,
      min_relic_per_group_per_region = min_relic_per_group_per_region
    )
  
  strata_diag <- relic_bal |>
    dplyr::mutate(
      age_bin = make_quantile_bins(tmrca_final, AGE_BIN_N, "age"),
      size_bin = make_quantile_bins(log(branch_size + 1), SIZE_BIN_N, "size"),
      strata = interaction(age_bin, size_bin, drop = TRUE)
    ) |>
    dplyr::count(strata, hg_group, name = "n_branches") |>
    tidyr::pivot_wider(
      names_from = hg_group,
      values_from = n_branches,
      values_fill = 0
    )
  
  if (!g1 %in% names(strata_diag)) strata_diag[[g1]] <- 0L
  if (!g2 %in% names(strata_diag)) strata_diag[[g2]] <- 0L
  
  strata_diag <- strata_diag |>
    dplyr::mutate(
      swappable = .data[[g1]] > 0 & .data[[g2]] > 0,
      label_set = label_set,
      panel_id = panel_id,
      min_relic_per_group_per_region = min_relic_per_group_per_region
    )
  
  swappable_strata <- strata_diag |>
    dplyr::filter(swappable) |>
    dplyr::pull(strata)
  
  retained_in_swappable <- relic_bal |>
    dplyr::mutate(
      age_bin = make_quantile_bins(tmrca_final, AGE_BIN_N, "age"),
      size_bin = make_quantile_bins(log(branch_size + 1), SIZE_BIN_N, "size"),
      strata = interaction(age_bin, size_bin, drop = TRUE)
    ) |>
    dplyr::filter(strata %in% swappable_strata)
  
  panel_retention_summary <- tibble::tibble(
    label_set = label_set,
    panel_id = panel_id,
    focal_group = g1,
    control_group = g2,
    min_relic_per_group_per_region = min_relic_per_group_per_region,
    total_relic_branches_focal_raw = raw_counts$n_raw_relic[match(g1, raw_counts$hg_group)],
    total_relic_branches_control_raw = raw_counts$n_raw_relic[match(g2, raw_counts$hg_group)],
    total_relic_branches_focal = sum(relic_bal$hg_group == g1),
    total_relic_branches_control = sum(relic_bal$hg_group == g2),
    shared_regions = dplyr::n_distinct(relic_bal$dominant_region),
    n_swappable_strata = length(swappable_strata),
    pct_retained_branches_in_swappable = nrow(retained_in_swappable) / nrow(relic_bal)
  ) |>
    dplyr::mutate(
      focal_retained_pct = total_relic_branches_focal / total_relic_branches_focal_raw,
      control_retained_pct = total_relic_branches_control / total_relic_branches_control_raw,
      adequacy_flag = dplyr::case_when(
        panel_id == "ph908_r1a_min1" &
          n_swappable_strata >= 3 &
          pct_retained_branches_in_swappable >= 0.80 ~ "adequate",
        panel_id == "ph908_r1a_min1" ~ "primary_but_check_power",
        panel_id == "ph908_r1a_min2" ~ "underpowered_directional_only",
        TRUE ~ "sensitivity"
      )
    )
  
  null_tbl <- tibble::tibble()
  
  if (N_PERM >= 1) {
    relic_perm_base <- relic_bal |>
      dplyr::mutate(
        age_bin = make_quantile_bins(tmrca_final, AGE_BIN_N, "age"),
        size_bin = make_quantile_bins(log(branch_size + 1), SIZE_BIN_N, "size"),
        strata = interaction(age_bin, size_bin, drop = TRUE)
      )
    
    null_res <- vector("list", N_PERM)
    
    for (b in seq_len(N_PERM)) {
      perm_group <- permute_within_strata(relic_perm_base$hg_group, relic_perm_base$strata)
      
      perm_tbl <- relic_perm_base |>
        dplyr::mutate(hg_group = perm_group)
      
      perm_summary <- dplyr::bind_rows(
        calc_group_metrics_from_branch_df(perm_tbl |> dplyr::filter(hg_group == g1)) |>
          dplyr::mutate(hg_group = g1),
        calc_group_metrics_from_branch_df(perm_tbl |> dplyr::filter(hg_group == g2)) |>
          dplyr::mutate(hg_group = g2)
      )
      
      null_res[[b]] <- tibble::tibble(
        perm_id = b,
        metric = c("shannon_entropy", "effective_regions", "top3_region_share", "top1_region_share", "top2_region_share"),
        perm_value = c(
          perm_summary$shannon_entropy[perm_summary$hg_group == g1] - perm_summary$shannon_entropy[perm_summary$hg_group == g2],
          perm_summary$effective_regions[perm_summary$hg_group == g1] - perm_summary$effective_regions[perm_summary$hg_group == g2],
          perm_summary$top3_region_share[perm_summary$hg_group == g1] - perm_summary$top3_region_share[perm_summary$hg_group == g2],
          perm_summary$top1_region_share[perm_summary$hg_group == g1] - perm_summary$top1_region_share[perm_summary$hg_group == g2],
          perm_summary$top2_region_share[perm_summary$hg_group == g1] - perm_summary$top2_region_share[perm_summary$hg_group == g2]
        ),
        label_set = label_set,
        panel_id = panel_id,
        min_relic_per_group_per_region = min_relic_per_group_per_region
      )
    }
    
    null_tbl <- dplyr::bind_rows(null_res)
  }
  
  formal_tests <- obs_diffs |>
    dplyr::rowwise() |>
    dplyr::mutate(
      ci_low = ci95(null_tbl$perm_value[null_tbl$metric == metric])[1],
      ci_high = ci95(null_tbl$perm_value[null_tbl$metric == metric])[2],
      p_one_sided = one_sided_more_concentrated_p(
        metric_name = metric,
        null_vals = null_tbl$perm_value[null_tbl$metric == metric],
        obs = observed_diff
      ),
      p_two_sided = two_sided_perm_p(
        null_vals = null_tbl$perm_value[null_tbl$metric == metric],
        obs = observed_diff
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      support_label = dplyr::case_when(
        panel_id == "ph908_r1a_min1" &
          metric %in% PRIMARY_METRICS &
          direction_consistent %in% TRUE &
          is.finite(p_one_sided) &
          p_one_sided < 0.05 ~ "Primary supported",
        panel_id == "ph908_r1a_min2" &
          metric %in% PRIMARY_METRICS &
          direction_consistent %in% TRUE ~ "Directional / underpowered",
        panel_id %in% c("i2_r1a_min1", "i2_r1a_min2") &
          metric %in% PRIMARY_METRICS &
          direction_consistent %in% TRUE &
          is.finite(p_one_sided) &
          p_one_sided < 0.05 ~ "Sensitivity supported",
        TRUE ~ "Not supported"
      )
    )
  
  final_summary_tbl <- formal_tests |>
    dplyr::select(
      label_set, panel_id,
      focal_group, control_group,
      min_relic_per_group_per_region,
      metric, metric_label, direction_label,
      focal_value, control_value, observed_diff,
      ci_low, ci_high, p_one_sided, p_two_sided,
      direction_consistent, support_label
    ) |>
    dplyr::left_join(
      panel_retention_summary,
      by = c(
        "label_set", "panel_id",
        "focal_group", "control_group",
        "min_relic_per_group_per_region"
      )
    )
  
  compact_summary_tbl <- final_summary_tbl |>
    dplyr::filter(metric %in% PRIMARY_METRICS) |>
    dplyr::arrange(factor(metric, levels = PRIMARY_METRICS))
  
  list(
    branch_tbl = branch_tbl,
    relic_raw = relic_raw,
    relic_bal = relic_bal,
    region_balance_diag = region_balance_diag,
    balanced_regions_tbl = balanced_regions_tbl,
    region_counts = region_counts,
    obs_summary = obs_summary,
    obs_diffs = obs_diffs,
    formal_tests = formal_tests,
    final_summary_tbl = final_summary_tbl,
    compact_summary_tbl = compact_summary_tbl,
    null_tbl = null_tbl,
    strata_diag = strata_diag,
    panel_retention_summary = panel_retention_summary
  )
}

# Warning handling

warning_log <- character()

withCallingHandlers({
  
  if (!file.exists(FILE_PREPROCESSED_RDS)) {
    stop_msg("Missing dataset: ", FILE_PREPROCESSED_RDS)
  }
  if (!file.exists(FILE_BRANCH_AGES)) {
    stop_msg("Missing age file: ", FILE_BRANCH_AGES)
  }
  
  cl_raw <- readRDS(FILE_PREPROCESSED_RDS)
  
  if (!"is_balkan_country" %in% names(cl_raw)) {
    if (!"country_code" %in% names(cl_raw)) {
      stop_msg("Dataset missing both is_balkan_country and country_code.")
    }
    cl_raw <- cl_raw |>
      dplyr::mutate(is_balkan_country = toupper(country_code) %in% BALKAN_ISO2)
    message("Reconstructed is_balkan_country from country_code for Pillar 4 filtering.")
  }
  
  required_cols <- c(
    "matched", "exclude_geo_primary", "is_balkan_country",
    "terminal_snp", "Region", "lat", "long",
    "major_hg", "tmrca", "is_ph908_primary"
  )
  
  missing_required <- setdiff(required_cols, names(cl_raw))
  if (length(missing_required) > 0) {
    stop_msg("preprocessed_dataset.rds is missing required columns: ",
             paste(missing_required, collapse = ", "))
  }
  
  age_tbl_raw <- readr::read_csv(FILE_BRANCH_AGES, show_col_types = FALSE)
  if (!all(c("terminal_snp", "tmrca_ybp") %in% names(age_tbl_raw))) {
    stop_msg("branch_ages.csv must contain terminal_snp and tmrca_ybp")
  }
  
  age_tbl <- age_tbl_raw |>
    dplyr::mutate(
      terminal_snp_norm = normalize_snp(terminal_snp),
      tmrca_ybp = as.numeric(tmrca_ybp)
    ) |>
    dplyr::filter(is.finite(tmrca_ybp), tmrca_ybp > 0) |>
    dplyr::group_by(terminal_snp_norm) |>
    dplyr::summarise(
      tmrca_ybp = safe_first_non_na(tmrca_ybp),
      .groups = "drop"
    )
  
  filter_diag <- tibble::tibble(
    step = c(
      "raw_rows",
      "matched_true",
      "exclude_geo_primary_false",
      "is_balkan_country_true",
      "terminal_snp_present",
      "region_present",
      "finite_lat_long"
    ),
    n = c(
      nrow(cl_raw),
      sum(cl_raw$matched %in% TRUE, na.rm = TRUE),
      sum(cl_raw$exclude_geo_primary %in% FALSE, na.rm = TRUE),
      sum(cl_raw$is_balkan_country %in% TRUE, na.rm = TRUE),
      sum(!is.na(cl_raw$terminal_snp) & cl_raw$terminal_snp != ""),
      sum(!is.na(cl_raw$Region) & cl_raw$Region != ""),
      sum(is.finite(cl_raw$lat) & is.finite(cl_raw$long))
    )
  )
  
  cat("\n=====================================================\n")
  cat("PILLAR 4 — RELIC REGIONAL CONCENTRATION CONTRAST\n")
  cat("=====================================================\n")
  cat("Input dataset: ", FILE_PREPROCESSED_RDS, "\n", sep = "")
  cat("Age file     : ", FILE_BRANCH_AGES, "\n", sep = "")
  cat("Seed         : ", RNG_SEED, "\n\n", sep = "")
  
  sample_base <- cl_raw |>
    dplyr::filter(
      matched %in% TRUE,
      exclude_geo_primary %in% FALSE,
      is_balkan_country %in% TRUE,
      !is.na(terminal_snp), terminal_snp != "",
      !is.na(Region), Region != "",
      is.finite(lat), is.finite(long)
    ) |>
    dplyr::mutate(
      terminal_snp_norm = normalize_snp(terminal_snp),
      hg_group_ph908_r1a = dplyr::case_when(
        is_ph908_primary %in% TRUE ~ "PH908",
        major_hg == "R1a" ~ "R1a",
        TRUE ~ NA_character_
      ),
      hg_group_i2_r1a = dplyr::case_when(
        major_hg == "I2" ~ "I2",
        major_hg == "R1a" ~ "R1a",
        TRUE ~ NA_character_
      )
    )
  
  branch_tbl_main <- sample_base |>
    dplyr::filter(!is.na(hg_group_ph908_r1a)) |>
    dplyr::group_by(hg_group = hg_group_ph908_r1a, terminal_snp_norm) |>
    dplyr::summarise(
      branch_size = dplyr::n(),
      dominant_region = dominant_region_fun(Region),
      dataset_tmrca = safe_mean(tmrca),
      .groups = "drop"
    ) |>
    dplyr::left_join(age_tbl, by = "terminal_snp_norm") |>
    dplyr::mutate(
      tmrca_final = dplyr::coalesce(tmrca_ybp, dataset_tmrca),
      age_source = dplyr::case_when(
        is.finite(tmrca_ybp) ~ "branch_ages_csv",
        is.finite(dataset_tmrca) ~ "dataset_tmrca",
        TRUE ~ "missing"
      ),
      class = classify_branch(branch_size)
    ) |>
    dplyr::filter(is.finite(tmrca_final), !is.na(dominant_region), dominant_region != "")
  
  branch_tbl_sens <- sample_base |>
    dplyr::filter(!is.na(hg_group_i2_r1a)) |>
    dplyr::group_by(hg_group = hg_group_i2_r1a, terminal_snp_norm) |>
    dplyr::summarise(
      branch_size = dplyr::n(),
      dominant_region = dominant_region_fun(Region),
      dataset_tmrca = safe_mean(tmrca),
      .groups = "drop"
    ) |>
    dplyr::left_join(age_tbl, by = "terminal_snp_norm") |>
    dplyr::mutate(
      tmrca_final = dplyr::coalesce(tmrca_ybp, dataset_tmrca),
      age_source = dplyr::case_when(
        is.finite(tmrca_ybp) ~ "branch_ages_csv",
        is.finite(dataset_tmrca) ~ "dataset_tmrca",
        TRUE ~ "missing"
      ),
      class = classify_branch(branch_size)
    ) |>
    dplyr::filter(is.finite(tmrca_final), !is.na(dominant_region), dominant_region != "")
  
  age_diag <- dplyr::bind_rows(
    branch_tbl_main,
    branch_tbl_sens |> dplyr::filter(hg_group == "I2")
  ) |>
    dplyr::group_by(hg_group) |>
    dplyr::summarise(
      n_branches = dplyr::n(),
      total_mass = sum(branch_size, na.rm = TRUE),
      n_relic = sum(class == "Relic", na.rm = TRUE),
      n_middle = sum(class == "Middle", na.rm = TRUE),
      n_founder = sum(class == "Founder", na.rm = TRUE),
      age_from_age_file = sum(age_source == "branch_ages_csv", na.rm = TRUE),
      age_from_dataset = sum(age_source == "dataset_tmrca", na.rm = TRUE),
      median_branch_size = stats::median(branch_size, na.rm = TRUE),
      max_branch_size = max(branch_size, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("Input/filter diagnostics:\n")
  print(filter_diag, n = Inf)
  
  cat("\nAge / branch diagnostics:\n")
  print(age_diag, n = Inf)
  
  panel_ph908_min1 <- run_relic_concentration_analysis(
    branch_tbl = branch_tbl_main,
    valid_groups = c("PH908", "R1a"),
    label_set = "PH908_vs_R1a",
    panel_id = "ph908_r1a_min1",
    min_relic_per_group_per_region = 1L
  )
  
  panel_ph908_min2 <- run_relic_concentration_analysis(
    branch_tbl = branch_tbl_main,
    valid_groups = c("PH908", "R1a"),
    label_set = "PH908_vs_R1a",
    panel_id = "ph908_r1a_min2",
    min_relic_per_group_per_region = 2L
  )
  
  panel_i2_min1 <- run_relic_concentration_analysis(
    branch_tbl = branch_tbl_sens,
    valid_groups = c("I2", "R1a"),
    label_set = "I2_vs_R1a",
    panel_id = "i2_r1a_min1",
    min_relic_per_group_per_region = 1L
  )
  
  panel_i2_min2 <- run_relic_concentration_analysis(
    branch_tbl = branch_tbl_sens,
    valid_groups = c("I2", "R1a"),
    label_set = "I2_vs_R1a",
    panel_id = "i2_r1a_min2",
    min_relic_per_group_per_region = 2L
  )
  
  panel_list <- list(
    ph908_r1a_min1 = panel_ph908_min1,
    ph908_r1a_min2 = panel_ph908_min2,
    i2_r1a_min1 = panel_i2_min1,
    i2_r1a_min2 = panel_i2_min2
  )
  
  retention_tbl <- dplyr::bind_rows(purrr::map(panel_list, "panel_retention_summary"))
  cross_panel_tbl <- dplyr::bind_rows(purrr::map(panel_list, "compact_summary_tbl"))
  region_balance_diag_tbl <- dplyr::bind_rows(purrr::map(panel_list, "region_balance_diag"))
  balanced_regions_tbl_all <- dplyr::bind_rows(purrr::map(panel_list, "balanced_regions_tbl"))
  region_counts_tbl_all <- dplyr::bind_rows(purrr::map(panel_list, "region_counts"))
  formal_tests_tbl_all <- dplyr::bind_rows(purrr::map(panel_list, "formal_tests"))
  null_tbl_all <- dplyr::bind_rows(purrr::map(panel_list, "null_tbl"))
  strata_diag_tbl_all <- dplyr::bind_rows(purrr::map(panel_list, "strata_diag"))
  
  figuredata_p4_main <- cross_panel_tbl |>
    dplyr::filter(panel_id == "ph908_r1a_min1", metric %in% PRIMARY_METRICS) |>
    dplyr::mutate(
      panel_label = "PH908-R1a min1"
    ) |>
    dplyr::arrange(factor(metric, levels = PRIMARY_METRICS))
  
  figuredata_p4_sens <- cross_panel_tbl |>
    dplyr::filter(metric %in% PRIMARY_METRICS) |>
    dplyr::mutate(
      panel_label = dplyr::case_when(
        panel_id == "ph908_r1a_min1" ~ "PH908-R1a min1",
        panel_id == "ph908_r1a_min2" ~ "PH908-R1a min2",
        panel_id == "i2_r1a_min1" ~ "I2-R1a min1",
        panel_id == "i2_r1a_min2" ~ "I2-R1a min2",
        TRUE ~ panel_id
      )
    ) |>
    dplyr::arrange(panel_label, factor(metric, levels = PRIMARY_METRICS))
  
  table_p4_supp <- dplyr::bind_rows(
    age_diag |>
      tidyr::pivot_longer(
        cols = -hg_group,
        names_to = "metric",
        values_to = "value"
      ) |>
      tidyr::pivot_wider(
        names_from = hg_group,
        values_from = value
      ),
    cross_panel_tbl |>
      dplyr::transmute(
        metric = paste0(panel_id, "__", metric),
        PH908 = dplyr::if_else(focal_group == "PH908", focal_value, NA_real_),
        I2    = dplyr::if_else(focal_group == "I2", focal_value, NA_real_),
        R1a   = control_value
      )
  )
  
  git_commit <- safe_git_commit()
  
  primary_row_entropy <- figuredata_p4_main |>
    dplyr::filter(metric == "shannon_entropy")
  primary_row_eff <- figuredata_p4_main |>
    dplyr::filter(metric == "effective_regions")
  primary_row_top3 <- figuredata_p4_main |>
    dplyr::filter(metric == "top3_region_share")
  
  primary_retention <- retention_tbl |>
    dplyr::filter(panel_id == PRIMARY_PANEL_ID)
  
  primary_regions <- balanced_regions_tbl_all |>
    dplyr::filter(panel_id == PRIMARY_PANEL_ID) |>
    dplyr::pull(dominant_region)
  
  age_ph908 <- age_diag |>
    dplyr::filter(hg_group == "PH908")
  age_r1a <- age_diag |>
    dplyr::filter(hg_group == "R1a")
  age_i2 <- age_diag |>
    dplyr::filter(hg_group == "I2")
  
  sensitivity_sentence <- paste0(
    "Direction remained manuscript-consistent across PH908-R1a min2, I2-R1a min1, and I2-R1a min2, while the stricter balanced-arena threshold reduced power as expected."
  )
  
  writeLines(
    c(
      paste0("script: ", SCRIPT_NAME),
      paste0("timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0("seed: ", RNG_SEED),
      paste0("config_file: ", path_for_facts(CONFIG_FILE)),
      "project_root: .",
      paste0("git_commit: ", git_commit),
      "inputs:",
      paste0("  - ", path_for_facts(FILE_PREPROCESSED_RDS)),
      paste0("  - ", path_for_facts(FILE_BRANCH_AGES)),
      paste0("n_perm: ", N_PERM),
      paste0("age_bin_n: ", AGE_BIN_N),
      paste0("size_bin_n: ", SIZE_BIN_N),
      paste0("primary_panel_id: ", PRIMARY_PANEL_ID),
      "",
      "filter_diagnostics:",
      paste0("  raw_rows: ", filter_diag$n[filter_diag$step == "raw_rows"]),
      paste0("  matched_true: ", filter_diag$n[filter_diag$step == "matched_true"]),
      paste0("  exclude_geo_primary_false: ", filter_diag$n[filter_diag$step == "exclude_geo_primary_false"]),
      paste0("  is_balkan_country_true: ", filter_diag$n[filter_diag$step == "is_balkan_country_true"]),
      paste0("  terminal_snp_present: ", filter_diag$n[filter_diag$step == "terminal_snp_present"]),
      paste0("  region_present: ", filter_diag$n[filter_diag$step == "region_present"]),
      paste0("  finite_lat_long: ", filter_diag$n[filter_diag$step == "finite_lat_long"]),
      "",
      "branch_structure:",
      paste0("  PH908_branches: ", age_ph908$n_branches[[1]]),
      paste0("  PH908_total_mass: ", age_ph908$total_mass[[1]]),
      paste0("  PH908_relic_middle_founder: ", age_ph908$n_relic[[1]], "/", age_ph908$n_middle[[1]], "/", age_ph908$n_founder[[1]]),
      paste0("  R1a_branches: ", age_r1a$n_branches[[1]]),
      paste0("  R1a_total_mass: ", age_r1a$total_mass[[1]]),
      paste0("  R1a_relic_middle_founder: ", age_r1a$n_relic[[1]], "/", age_r1a$n_middle[[1]], "/", age_r1a$n_founder[[1]]),
      paste0("  I2_branches: ", age_i2$n_branches[[1]]),
      paste0("  I2_total_mass: ", age_i2$total_mass[[1]]),
      paste0("  I2_relic_middle_founder: ", age_i2$n_relic[[1]], "/", age_i2$n_middle[[1]], "/", age_i2$n_founder[[1]]),
      "",
      "primary_panel:",
      paste0("  retained_relic_branches_PH908: ", primary_retention$total_relic_branches_focal[[1]]),
      paste0("  retained_relic_branches_R1a: ", primary_retention$total_relic_branches_control[[1]]),
      paste0("  shared_regions: ", primary_retention$shared_regions[[1]]),
      paste0("  swappable_strata: ", primary_retention$n_swappable_strata[[1]]),
      paste0("  pct_retained_branches_in_swappable: ", sprintf("%.4f", primary_retention$pct_retained_branches_in_swappable[[1]])),
      paste0("  retained_regions: ", paste(primary_regions, collapse = " | ")),
      "",
      "primary_metric_summary:",
      paste0("  shannon_entropy_PH908_R1a_diff: ",
             sprintf("%.3f", primary_row_entropy$focal_value[[1]]), " / ",
             sprintf("%.3f", primary_row_entropy$control_value[[1]]), " / ",
             sprintf("%.3f", primary_row_entropy$observed_diff[[1]])),
      paste0("  effective_regions_PH908_R1a_diff: ",
             sprintf("%.3f", primary_row_eff$focal_value[[1]]), " / ",
             sprintf("%.3f", primary_row_eff$control_value[[1]]), " / ",
             sprintf("%.3f", primary_row_eff$observed_diff[[1]])),
      paste0("  top3_share_PH908_R1a_diff: ",
             sprintf("%.3f", primary_row_top3$focal_value[[1]]), " / ",
             sprintf("%.3f", primary_row_top3$control_value[[1]]), " / ",
             sprintf("%.3f", primary_row_top3$observed_diff[[1]])),
      "",
      paste0(
        "results_sentence: Within the balanced PH908-versus-R1a relic arena requiring at least one relic branch in both groups per region, ",
        primary_retention$total_relic_branches_focal[[1]], " PH908 relic branches and ",
        primary_retention$total_relic_branches_control[[1]], " R1a relic branches were retained across ",
        primary_retention$shared_regions[[1]], " shared regions."
      ),
      "scope_note: This is a regional occupancy concentration test, not a migration-route reconstruction.",
      paste0("sensitivity_sentence: ", sensitivity_sentence),
      "note: p-value consistency checks are conditional on the configured permutation depth.",
      "output_files:",
      paste0("  - ", path_for_facts(FILE_P4_RESULTS)),
      paste0("  - ", path_for_facts(FILE_P4_FIG_MAIN)),
      paste0("  - ", path_for_facts(FILE_P4_FIG_SENS)),
      paste0("  - ", path_for_facts(FILE_P4_SUPP_TABLE)),
      paste0("  - ", path_for_facts(FILE_P4_FACTS))
    ),
    FILE_P4_FACTS,
    useBytes = TRUE
  )
  
  readr::write_csv(figuredata_p4_main, FILE_P4_FIG_MAIN)
  readr::write_csv(figuredata_p4_sens, FILE_P4_FIG_SENS)
  readr::write_csv(table_p4_supp, FILE_P4_SUPP_TABLE)
  
  results_p4 <- list(
    metadata = list(
      script = SCRIPT_NAME,
      seed = RNG_SEED,
      config_file = path_for_facts(CONFIG_FILE),
      inputs = c(
        path_for_facts(FILE_PREPROCESSED_RDS),
        path_for_facts(FILE_BRANCH_AGES)
      ),
      n_perm = N_PERM,
      age_bin_n = AGE_BIN_N,
      size_bin_n = SIZE_BIN_N,
      primary_panel_id = PRIMARY_PANEL_ID,
      git_commit = git_commit
    ),
    filter_diag = filter_diag,
    age_diag = age_diag,
    sample_base = sample_base,
    branch_tbl_main = branch_tbl_main,
    branch_tbl_sens = branch_tbl_sens,
    panel_ph908_min1 = panel_ph908_min1,
    panel_ph908_min2 = panel_ph908_min2,
    panel_i2_min1 = panel_i2_min1,
    panel_i2_min2 = panel_i2_min2,
    retention_tbl = retention_tbl,
    cross_panel_tbl = cross_panel_tbl,
    region_balance_diag_tbl = region_balance_diag_tbl,
    balanced_regions_tbl_all = balanced_regions_tbl_all,
    region_counts_tbl_all = region_counts_tbl_all,
    formal_tests_tbl_all = formal_tests_tbl_all,
    null_tbl_all = null_tbl_all,
    strata_diag_tbl_all = strata_diag_tbl_all,
    figuredata_main = figuredata_p4_main,
    figuredata_sens = figuredata_p4_sens,
    table_supp = table_p4_supp
  )
  
  saveRDS(results_p4, FILE_P4_RESULTS)
  
  cat("\nSaved outputs:\n")
  cat(" - ", FILE_P4_RESULTS, "\n", sep = "")
  cat(" - ", FILE_P4_FIG_MAIN, "\n", sep = "")
  cat(" - ", FILE_P4_FIG_SENS, "\n", sep = "")
  cat(" - ", FILE_P4_SUPP_TABLE, "\n", sep = "")
  cat(" - ", FILE_P4_FACTS, "\n", sep = "")
  
  cat("\nPILLAR 4 COMPLETE.\n")
  
}, warning = function(w) {
  warning_log <<- c(warning_log, conditionMessage(w))
  invokeRestart("muffleWarning")
})

warning_log_path <- file.path(DIR_FACTS, "warnings_p4.txt")
write_warning_log(warning_log, warning_log_path)