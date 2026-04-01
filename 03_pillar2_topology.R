# 03_pillar2_topology.R
# Cross-region branch-link topology contrast for PH908 versus R1a.

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
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
        paste(candidates, collapse = "\n  - "),
        "\n\nRun the script from the project root or place it beside config_analysis.R."
      ),
      call. = FALSE
    )
  }
  
  normalizePath(hit, winslash = "/", mustWork = TRUE)
}

CONFIG_FILE <- locate_config_file()
source(CONFIG_FILE)
setwd(DIR_ROOT)
set.seed(SEED_MAIN)

# Output paths

FILE_P2_RESULTS <- file.path(DIR_RESULTS, "results_p2.rds")
FILE_P2_FIG_MAIN <- file.path(DIR_FIGUREDATA, "figuredata_p2_main.csv")
FILE_P2_FIG_NULL <- file.path(DIR_FIGUREDATA, "figuredata_p2_null.csv")
FILE_P2_SUPP_TABLE <- file.path(DIR_TABLES, "table_p2_supp.csv")
FILE_P2_RX_TABLE <- file.path(DIR_TABLES, "table_p2_rx.csv")
FILE_P2_FACTS <- file.path(DIR_FACTS, "facts_p2.txt")

# Analysis settings

RULE_LEVELS <- c("nearest_older", "nearest_older_window", "distance_age_cost")
RULE_LABELS <- c(
  nearest_older = "Nearest older branch",
  nearest_older_window = "Nearest older within age window",
  distance_age_cost = "Distance + age cost"
)

RELIC_MAX_N <- 3L
FOUNDER_MIN_N <- 15L
PH908_GATING_MODE <- "primary"

MAX_AGE_GAP_WINDOW <- 1200
LAMBDA_AGE_COST <- 1.0

AGE_BIN_N <- 4L
SIZE_BIN_N <- 4L
GEO_LAT_BIN_N <- 3L
GEO_LON_BIN_N <- 3L

# Full manuscript-grade permutation settings.
N_PERM_PRIMARY <- 5000L
N_PERM_GEO_SENSITIVITY <- 5000L

ROBUSTNESS_SPECS <- tibble::tribble(
  ~perturbation,
  "baseline",
  "drop_largest_PH908",
  "drop_largest_R1a",
  "drop_top2_cumulative"
)

required_inputs <- c(FILE_PREPROCESSED_RDS, FILE_BRANCH_AGES)
missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs) > 0) {
  stop(
    paste0(
      "Missing required input file(s):\n  - ",
      paste(missing_inputs, collapse = "\n  - ")
    ),
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("03 PILLAR 2 — CROSS-REGION BRANCH-LINK TOPOLOGY\n")
cat("============================================================\n")
cat("Root: ", DIR_ROOT, "\n", sep = "")
cat("Seed: ", SEED_MAIN, "\n", sep = "")
cat("Inputs:\n")
cat("  - ", FILE_PREPROCESSED_RDS, "\n", sep = "")
cat("  - ", FILE_BRANCH_AGES, "\n", sep = "")
cat("Primary permutations: ", N_PERM_PRIMARY, "\n", sep = "")
cat("Geo sensitivity permutations: ", N_PERM_GEO_SENSITIVITY, "\n\n", sep = "")

# Helpers

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

git_commit <- safe_git_commit()

normalize_snp <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^I-|^R-", "", x)
  x <- gsub("\\*$", "", x)
  x
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

standardize_age_table <- function(age_raw) {
  nm <- names(age_raw)
  
  snp_col_candidates <- c("terminal_snp", "terminal_snp_norm", "snp", "branch", "node")
  age_col_candidates <- c("tmrca_ybp", "tmrca", "age_ybp", "ybp")
  
  snp_col <- snp_col_candidates[snp_col_candidates %in% nm][1]
  age_col <- age_col_candidates[age_col_candidates %in% nm][1]
  
  if (is.na(snp_col) || is.na(age_col)) {
    stop(
      paste0(
        "branch_ages.csv must contain a terminal-SNP column and an age column. ",
        "Accepted SNP columns: ", paste(snp_col_candidates, collapse = ", "),
        ". Accepted age columns: ", paste(age_col_candidates, collapse = ", "),
        ". Found columns: ", paste(nm, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  age_raw |>
    dplyr::transmute(
      terminal_snp_norm = normalize_snp(.data[[snp_col]]),
      tmrca_ybp = suppressWarnings(as.numeric(.data[[age_col]]))
    ) |>
    dplyr::filter(!is.na(terminal_snp_norm), terminal_snp_norm != "", is.finite(tmrca_ybp)) |>
    dplyr::group_by(terminal_snp_norm) |>
    dplyr::summarise(
      tmrca_ybp = dplyr::first(tmrca_ybp),
      .groups = "drop"
    )
}

assert_true <- function(flag, msg) {
  if (!isTRUE(flag)) stop(msg, call. = FALSE)
}

assert_equal <- function(actual, expected, msg) {
  actual_cmp <- suppressWarnings(as.numeric(actual))
  expected_cmp <- suppressWarnings(as.numeric(expected))
  
  ok <- length(actual) == 1 &&
    length(expected) == 1 &&
    !is.na(actual_cmp) &&
    !is.na(expected_cmp) &&
    actual_cmp == expected_cmp
  
  if (!ok) {
    stop(
      paste0(
        msg,
        "\n  actual: ", paste(actual, collapse = ", "),
        "\n  expected: ", paste(expected, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

assert_close <- function(actual, expected, tol, msg) {
  if (length(actual) != 1 || !is.finite(actual) || abs(actual - expected) > tol) {
    stop(
      paste0(
        msg,
        "\n  actual: ", format(actual, digits = 10),
        "\n  expected: ", format(expected, digits = 10),
        "\n  tol: ", tol
      ),
      call. = FALSE
    )
  }
}

make_quantile_bins <- function(x, n_bins, prefix) {
  x <- as.numeric(x)
  
  if (sum(is.finite(x)) < 2) {
    return(rep(paste0(prefix, "_1"), length(x)))
  }
  
  probs <- seq(0, 1, length.out = n_bins + 1)
  qs <- unique(stats::quantile(x[is.finite(x)], probs = probs, na.rm = TRUE))
  
  if (length(qs) < 2) {
    return(rep(paste0(prefix, "_1"), length(x)))
  }
  
  out <- cut(x, breaks = qs, include.lowest = TRUE, ordered_result = TRUE)
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
      out[idx] <- sample(out[idx], length(idx), replace = FALSE)
    }
  }
  
  out
}

haversine_km <- function(lat1, lon1, lat2, lon2) {
  r <- 6371
  to_rad <- pi / 180
  
  lat1 <- lat1 * to_rad
  lon1 <- lon1 * to_rad
  lat2 <- lat2 * to_rad
  lon2 <- lon2 * to_rad
  
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  2 * r * asin(pmin(1, sqrt(a)))
}

choose_parent_index <- function(current_row, older_df, rule) {
  if (nrow(older_df) == 0) return(NA_integer_)
  
  dist_vec <- haversine_km(
    current_row$lat_c[[1]],
    current_row$lon_c[[1]],
    older_df$lat_c,
    older_df$lon_c
  )
  
  age_gap_vec <- older_df$tmrca_final - current_row$tmrca_final[[1]]
  valid <- is.finite(dist_vec) & is.finite(age_gap_vec) & age_gap_vec >= 0
  
  if (!any(valid)) return(NA_integer_)
  
  dist_vec_valid <- dist_vec[valid]
  age_gap_valid <- age_gap_vec[valid]
  valid_idx <- which(valid)
  
  if (rule == "nearest_older") {
    return(valid_idx[which.min(dist_vec_valid)])
  }
  
  if (rule == "nearest_older_window") {
    in_window <- age_gap_valid <= MAX_AGE_GAP_WINDOW
    if (any(in_window)) {
      return(valid_idx[in_window][which.min(dist_vec_valid[in_window])])
    }
    return(valid_idx[which.min(dist_vec_valid)])
  }
  
  if (rule == "distance_age_cost") {
    dist_scale <- stats::median(dist_vec_valid, na.rm = TRUE)
    age_scale <- stats::median(age_gap_valid, na.rm = TRUE)
    
    if (!is.finite(dist_scale) || dist_scale <= 0) dist_scale <- 1
    if (!is.finite(age_scale) || age_scale <= 0) age_scale <- 1
    
    score <- (dist_vec_valid / dist_scale) + LAMBDA_AGE_COST * (age_gap_valid / age_scale)
    return(valid_idx[which.min(score)])
  }
  
  stop(paste0("Unknown rule: ", rule), call. = FALSE)
}

write_warning_log <- function(warning_vec, path) {
  if (length(warning_vec) == 0) {
    writeLines("No warnings generated.", path)
  } else {
    writeLines(unique(warning_vec), path)
  }
}

get_main_row <- function(tbl, rule_name) {
  tbl |>
    dplyr::filter(rule == rule_name)
}

build_branch_table <- function(base_rows, age_tbl) {
  base_rows |>
    dplyr::group_by(hg_group, terminal_snp_norm) |>
    dplyr::summarise(
      branch_size = dplyr::n(),
      lat_c = mean(lat, na.rm = TRUE),
      lon_c = mean(long, na.rm = TRUE),
      dominant_region = dominant_region_fun(Region),
      dataset_tmrca = safe_mean(tmrca),
      .groups = "drop"
    ) |>
    dplyr::left_join(age_tbl, by = "terminal_snp_norm") |>
    dplyr::mutate(
      tmrca_final = dplyr::coalesce(tmrca_ybp, dataset_tmrca),
      age_source = dplyr::case_when(
        !is.na(tmrca_ybp) ~ "external_age",
        is.na(tmrca_ybp) & !is.na(dataset_tmrca) ~ "dataset_age",
        TRUE ~ "missing"
      ),
      class_now = classify_branch(branch_size)
    ) |>
    dplyr::filter(is.finite(tmrca_final))
}

build_network_by_rule <- function(df_group, rule) {
  df_group <- df_group |>
    dplyr::filter(
      is.finite(tmrca_final),
      is.finite(lat_c),
      is.finite(lon_c),
      !is.na(dominant_region),
      dominant_region != ""
    ) |>
    dplyr::arrange(
      dplyr::desc(tmrca_final),
      dplyr::desc(branch_size),
      terminal_snp_norm
    )
  
  if (nrow(df_group) < 2) {
    return(tibble::tibble(
      parent_snp = character(),
      child_snp = character(),
      parent_tmrca = numeric(),
      child_tmrca = numeric(),
      age_gap = numeric(),
      edge_km = numeric(),
      parent_region = character(),
      child_region = character(),
      same_region = logical(),
      cross_region = logical(),
      rule = character(),
      hg_group = character()
    ))
  }
  
  edges <- vector("list", nrow(df_group) - 1)
  
  for (i in seq.int(2, nrow(df_group))) {
    current <- df_group[i, , drop = FALSE]
    older <- df_group[seq_len(i - 1), , drop = FALSE]
    
    j <- choose_parent_index(current, older, rule)
    
    if (!is.na(j)) {
      edge_km <- haversine_km(
        current$lat_c[[1]], current$lon_c[[1]],
        older$lat_c[[j]], older$lon_c[[j]]
      )
      
      edges[[i - 1]] <- tibble::tibble(
        parent_snp = older$terminal_snp_norm[[j]],
        child_snp = current$terminal_snp_norm[[1]],
        parent_tmrca = older$tmrca_final[[j]],
        child_tmrca = current$tmrca_final[[1]],
        age_gap = older$tmrca_final[[j]] - current$tmrca_final[[1]],
        edge_km = edge_km,
        parent_region = older$dominant_region[[j]],
        child_region = current$dominant_region[[1]],
        same_region = older$dominant_region[[j]] == current$dominant_region[[1]],
        cross_region = older$dominant_region[[j]] != current$dominant_region[[1]],
        rule = rule,
        hg_group = current$hg_group[[1]]
      )
    }
  }
  
  dplyr::bind_rows(edges)
}

build_all_rules_observed <- function(branch_tbl) {
  out <- vector("list", length(RULE_LEVELS) * 2)
  k <- 1
  
  for (rr in RULE_LEVELS) {
    for (gg in c("PH908", "R1a")) {
      out[[k]] <- build_network_by_rule(
        branch_tbl |>
          dplyr::filter(hg_group == gg),
        rr
      )
      k <- k + 1
    }
  }
  
  dplyr::bind_rows(out)
}

compute_rule_metrics <- function(edges_df) {
  edges_df |>
    dplyr::group_by(rule, hg_group) |>
    dplyr::summarise(
      n_edges = dplyr::n(),
      cross_region_share = mean(cross_region, na.rm = TRUE),
      same_region_share = mean(same_region, na.rm = TRUE),
      median_edge_km = stats::median(edge_km, na.rm = TRUE),
      median_age_gap = stats::median(age_gap, na.rm = TRUE),
      .groups = "drop"
    )
}

compute_diff_table <- function(metrics_df) {
  metrics_df |>
    dplyr::select(rule, hg_group, cross_region_share, n_edges) |>
    tidyr::pivot_wider(
      names_from = hg_group,
      values_from = c(cross_region_share, n_edges)
    ) |>
    dplyr::mutate(
      cross_region_share_diff = cross_region_share_PH908 - cross_region_share_R1a
    ) |>
    dplyr::arrange(factor(rule, levels = RULE_LEVELS))
}

make_primary_permutation_table <- function(branch_tbl_main) {
  perm_base <- branch_tbl_main |>
    dplyr::mutate(
      age_bin = make_quantile_bins(tmrca_final, AGE_BIN_N, "age"),
      size_bin = make_quantile_bins(log(branch_size + 1), SIZE_BIN_N, "size"),
      strata = interaction(age_bin, size_bin, class_now, drop = TRUE)
    )
  
  strata_diag <- perm_base |>
    dplyr::count(strata, hg_group, name = "n_branches") |>
    dplyr::arrange(strata, hg_group)
  
  null_list <- vector("list", N_PERM_PRIMARY)
  
  for (b in seq_len(N_PERM_PRIMARY)) {
    perm_group <- permute_within_strata(perm_base$hg_group, perm_base$strata)
    
    perm_tbl <- perm_base |>
      dplyr::mutate(hg_group = perm_group) |>
      dplyr::select(dplyr::all_of(names(branch_tbl_main)))
    
    perm_edges <- build_all_rules_observed(perm_tbl)
    perm_metrics <- compute_rule_metrics(perm_edges)
    perm_diff <- compute_diff_table(perm_metrics) |>
      dplyr::mutate(
        perm_id = b,
        null_model = "age_size_class"
      )
    
    null_list[[b]] <- perm_diff
  }
  
  list(
    strata_diag = strata_diag,
    null_tbl = dplyr::bind_rows(null_list)
  )
}

make_geo_permutation_table <- function(branch_tbl_main) {
  perm_geo_base <- branch_tbl_main |>
    dplyr::mutate(
      age_bin = make_quantile_bins(tmrca_final, AGE_BIN_N, "age"),
      size_bin = make_quantile_bins(log(branch_size + 1), SIZE_BIN_N, "size"),
      lat_bin = make_quantile_bins(lat_c, GEO_LAT_BIN_N, "lat"),
      lon_bin = make_quantile_bins(lon_c, GEO_LON_BIN_N, "lon"),
      strata = interaction(age_bin, size_bin, class_now, lat_bin, lon_bin, drop = TRUE)
    )
  
  strata_geo_diag <- perm_geo_base |>
    dplyr::count(strata, hg_group, name = "n_branches") |>
    dplyr::arrange(strata, hg_group)
  
  null_geo_list <- vector("list", N_PERM_GEO_SENSITIVITY)
  
  for (b in seq_len(N_PERM_GEO_SENSITIVITY)) {
    perm_group <- permute_within_strata(perm_geo_base$hg_group, perm_geo_base$strata)
    
    perm_tbl <- perm_geo_base |>
      dplyr::mutate(hg_group = perm_group) |>
      dplyr::select(dplyr::all_of(names(branch_tbl_main)))
    
    perm_edges <- build_all_rules_observed(perm_tbl)
    perm_metrics <- compute_rule_metrics(perm_edges)
    perm_diff <- compute_diff_table(perm_metrics) |>
      dplyr::mutate(
        perm_id = b,
        null_model = "age_size_class_geography"
      )
    
    null_geo_list[[b]] <- perm_diff
  }
  
  list(
    strata_geo_diag = strata_geo_diag,
    null_geo_tbl = dplyr::bind_rows(null_geo_list)
  )
}

compute_test_table <- function(obs_diff_tbl, null_tbl, null_model_name) {
  dplyr::bind_rows(lapply(RULE_LEVELS, function(rr) {
    obs_val <- obs_diff_tbl |>
      dplyr::filter(rule == rr) |>
      dplyr::pull(cross_region_share_diff)
    
    null_vals <- null_tbl |>
      dplyr::filter(rule == rr) |>
      dplyr::pull(cross_region_share_diff)
    
    tibble::tibble(
      rule = rr,
      null_model = null_model_name,
      test = "Cross_region_share_diff",
      statistic = obs_val,
      p_value_one_sided_PH908_less = (1 + sum(null_vals <= obs_val, na.rm = TRUE)) / (length(null_vals) + 1),
      p_value_two_sided = (1 + sum(abs(null_vals) >= abs(obs_val), na.rm = TRUE)) / (length(null_vals) + 1)
    )
  }))
}

make_robustness_branch_table <- function(branch_tbl_main, perturbation) {
  tbl <- branch_tbl_main
  
  if (perturbation == "baseline") {
    return(tbl)
  }
  
  if (perturbation == "drop_largest_PH908") {
    drop_row <- tbl |>
      dplyr::filter(hg_group == "PH908") |>
      dplyr::arrange(dplyr::desc(branch_size), terminal_snp_norm) |>
      dplyr::slice(1) |>
      dplyr::mutate(.drop_flag = TRUE) |>
      dplyr::select(hg_group, terminal_snp_norm, .drop_flag)
    
    return(
      tbl |>
        dplyr::left_join(drop_row, by = c("hg_group", "terminal_snp_norm")) |>
        dplyr::filter(is.na(.drop_flag)) |>
        dplyr::select(-.drop_flag)
    )
  }
  
  if (perturbation == "drop_largest_R1a") {
    drop_row <- tbl |>
      dplyr::filter(hg_group == "R1a") |>
      dplyr::arrange(dplyr::desc(branch_size), terminal_snp_norm) |>
      dplyr::slice(1) |>
      dplyr::mutate(.drop_flag = TRUE) |>
      dplyr::select(hg_group, terminal_snp_norm, .drop_flag)
    
    return(
      tbl |>
        dplyr::left_join(drop_row, by = c("hg_group", "terminal_snp_norm")) |>
        dplyr::filter(is.na(.drop_flag)) |>
        dplyr::select(-.drop_flag)
    )
  }
  
  if (perturbation == "drop_top2_cumulative") {
    drop_rows <- tbl |>
      dplyr::arrange(dplyr::desc(branch_size), hg_group, terminal_snp_norm) |>
      dplyr::slice_head(n = 2) |>
      dplyr::mutate(.drop_flag = TRUE) |>
      dplyr::select(hg_group, terminal_snp_norm, .drop_flag)
    
    return(
      tbl |>
        dplyr::left_join(drop_rows, by = c("hg_group", "terminal_snp_norm")) |>
        dplyr::filter(is.na(.drop_flag)) |>
        dplyr::select(-.drop_flag)
    )
  }
  
  stop(paste0("Unknown perturbation: ", perturbation), call. = FALSE)
}

# Input preparation

warning_log <- character()

withCallingHandlers({
  
  cl_raw <- readRDS(FILE_PREPROCESSED_RDS)
  
  if (!"is_balkan_country" %in% names(cl_raw)) {
    if (!"country_code" %in% names(cl_raw)) {
      stop(
        "Dataset is missing both is_balkan_country and country_code, so the Pillar 2 Balkan-country filter cannot be reconstructed.",
        call. = FALSE
      )
    }
    cl_raw <- cl_raw |>
      dplyr::mutate(
        is_balkan_country = toupper(country_code) %in% BALKAN_ISO2
      )
    message("Reconstructed is_balkan_country from country_code for Pillar 2 filtering.")
  }
  
  required_cols <- c(
    "matched", "exclude_geo_primary", "is_balkan_country",
    "terminal_snp", "major_hg",
    "is_ph908_primary", "is_ph908_accurate", "is_ph908_wide",
    "tmrca", "Region", "lat", "long"
  )
  missing_cols <- setdiff(required_cols, names(cl_raw))
  if (length(missing_cols) > 0) {
    stop(
      paste0("Dataset missing required columns after Pillar 2 input preparation: ", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  
  age_tbl <- standardize_age_table(readr::read_csv(FILE_BRANCH_AGES, show_col_types = FALSE))
  
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
      sum(is.finite(cl_raw$lat) & is.finite(cl_raw$long), na.rm = TRUE)
    )
  )
  
  base_rows <- cl_raw |>
    dplyr::filter(
      matched == TRUE,
      exclude_geo_primary == FALSE,
      is_balkan_country == TRUE,
      !is.na(terminal_snp), terminal_snp != "",
      !is.na(Region), Region != "",
      is.finite(lat), is.finite(long)
    ) |>
    dplyr::mutate(
      terminal_snp_norm = normalize_snp(terminal_snp),
      ph908_flag = dplyr::case_when(
        PH908_GATING_MODE == "primary"  ~ is_ph908_primary,
        PH908_GATING_MODE == "accurate" ~ is_ph908_accurate,
        PH908_GATING_MODE == "wide"     ~ is_ph908_wide,
        TRUE ~ FALSE
      ),
      hg_group = dplyr::case_when(
        ph908_flag ~ "PH908",
        major_hg == "R1a" ~ "R1a",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(hg_group))
  
  branch_tbl_main <- build_branch_table(base_rows, age_tbl)
  
  age_diag <- branch_tbl_main |>
    dplyr::group_by(hg_group) |>
    dplyr::summarise(
      n_branches = dplyr::n(),
      total_mass = sum(branch_size, na.rm = TRUE),
      n_relic = sum(class_now == "Relic", na.rm = TRUE),
      n_middle = sum(class_now == "Middle", na.rm = TRUE),
      n_founder = sum(class_now == "Founder", na.rm = TRUE),
      age_from_age_file = sum(age_source == "external_age", na.rm = TRUE),
      age_from_dataset = sum(age_source == "dataset_age", na.rm = TRUE),
      pct_age_from_age_file = 100 * mean(age_source == "external_age", na.rm = TRUE),
      pct_age_from_dataset = 100 * mean(age_source == "dataset_age", na.rm = TRUE),
      median_branch_size = stats::median(branch_size, na.rm = TRUE),
      max_branch_size = max(branch_size, na.rm = TRUE),
      .groups = "drop"
    )
  
  age_ph908 <- age_diag |>
    dplyr::filter(hg_group == "PH908")
  age_r1a <- age_diag |>
    dplyr::filter(hg_group == "R1a")
  
  branch_geo_diag <- branch_tbl_main |>
    dplyr::group_by(hg_group) |>
    dplyr::summarise(
      n_branches = dplyr::n(),
      mean_lat = mean(lat_c, na.rm = TRUE),
      mean_lon = mean(lon_c, na.rm = TRUE),
      median_lat = stats::median(lat_c, na.rm = TRUE),
      median_lon = stats::median(lon_c, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Observed branch-link networks
  
  edges_obs <- build_all_rules_observed(branch_tbl_main)
  
  network_diag <- edges_obs |>
    dplyr::group_by(rule, hg_group) |>
    dplyr::summarise(
      n_edges_output = dplyr::n(),
      median_edge_km = stats::median(edge_km, na.rm = TRUE),
      median_age_gap = stats::median(age_gap, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      branch_tbl_main |>
        dplyr::count(hg_group, name = "n_branches_input"),
      by = "hg_group"
    ) |>
    dplyr::mutate(
      n_unmatched_branches = n_branches_input - 1L - n_edges_output
    ) |>
    dplyr::arrange(factor(rule, levels = RULE_LEVELS), hg_group)
  
  metrics_obs <- compute_rule_metrics(edges_obs) |>
    dplyr::mutate(
      rule_label = unname(RULE_LABELS[rule])
    ) |>
    dplyr::arrange(factor(rule, levels = RULE_LEVELS), hg_group)
  
  obs_diff_tbl <- compute_diff_table(metrics_obs) |>
    dplyr::mutate(
      rule_label = unname(RULE_LABELS[rule])
    )
  
  # Null-model permutations
  
  primary_null <- make_primary_permutation_table(branch_tbl_main)
  strata_diag <- primary_null$strata_diag
  null_tbl <- primary_null$null_tbl
  
  geo_null <- make_geo_permutation_table(branch_tbl_main)
  strata_geo_diag <- geo_null$strata_geo_diag
  null_geo_tbl <- geo_null$null_geo_tbl
  
  tests_primary <- compute_test_table(obs_diff_tbl, null_tbl, "age_size_class")
  tests_geo <- compute_test_table(obs_diff_tbl, null_geo_tbl, "age_size_class_geography")
  
  tests_tbl <- dplyr::bind_rows(tests_primary, tests_geo) |>
    dplyr::arrange(factor(rule, levels = RULE_LEVELS), null_model)
  
  decision_tbl <- tests_tbl |>
    dplyr::filter(
      null_model == "age_size_class",
      test == "Cross_region_share_diff"
    ) |>
    dplyr::mutate(
      supports_PH908_less_cross_region = p_value_one_sided_PH908_less < 0.05
    ) |>
    dplyr::select(
      rule,
      statistic,
      p_value_one_sided_PH908_less,
      supports_PH908_less_cross_region
    )
  
  geo_sensitivity_tbl <- tests_tbl |>
    dplyr::filter(
      null_model == "age_size_class_geography",
      test == "Cross_region_share_diff"
    ) |>
    dplyr::select(
      rule,
      statistic,
      p_value_one_sided_PH908_less,
      p_value_two_sided
    )
  
  final_summary_tbl <- metrics_obs |>
    dplyr::select(
      rule,
      hg_group,
      cross_region_share,
      n_edges
    ) |>
    tidyr::pivot_wider(
      names_from = hg_group,
      values_from = c(cross_region_share, n_edges)
    ) |>
    dplyr::left_join(
      tests_tbl |>
        dplyr::filter(
          null_model == "age_size_class",
          test == "Cross_region_share_diff"
        ) |>
        dplyr::select(
          rule,
          cross_region_share_diff = statistic,
          p_value_one_sided_PH908_less,
          p_value_two_sided
        ),
      by = "rule"
    ) |>
    dplyr::left_join(
      decision_tbl |>
        dplyr::select(rule, supports_PH908_less_cross_region),
      by = "rule"
    ) |>
    dplyr::mutate(
      rule_label = dplyr::case_when(
        rule == "nearest_older" ~ "Nearest older branch",
        rule == "nearest_older_window" ~ "Nearest older within age window",
        rule == "distance_age_cost" ~ "Distance + age cost",
        TRUE ~ rule
      ),
      support_label = dplyr::case_when(
        supports_PH908_less_cross_region %in% TRUE ~ "Supported",
        supports_PH908_less_cross_region %in% FALSE ~ "Not supported",
        TRUE ~ "NA"
      )
    ) |>
    dplyr::select(
      rule,
      rule_label,
      n_edges_PH908,
      n_edges_R1a,
      cross_region_share_PH908,
      cross_region_share_R1a,
      cross_region_share_diff,
      p_value_one_sided_PH908_less,
      p_value_two_sided,
      support_label
    ) |>
    dplyr::arrange(factor(rule, levels = RULE_LEVELS))
  
  figuredata_p2_main <- final_summary_tbl |>
    dplyr::left_join(
      geo_sensitivity_tbl |>
        dplyr::rename(
          p_value_one_sided_geo = p_value_one_sided_PH908_less,
          p_value_two_sided_geo = p_value_two_sided
        ),
      by = "rule"
    ) |>
    dplyr::mutate(
      primary_support = support_label
    )
  
  figuredata_p2_null <- dplyr::bind_rows(
    null_tbl |>
      dplyr::select(rule, perm_id, cross_region_share_diff, null_model),
    null_geo_tbl |>
      dplyr::select(rule, perm_id, cross_region_share_diff, null_model)
  ) |>
    dplyr::mutate(
      rule_label = unname(RULE_LABELS[rule])
    )
  
  table_p2_supp <- dplyr::bind_rows(
    age_diag |>
      dplyr::mutate(section = "age_diag") |>
      tidyr::pivot_longer(
        cols = -c(section, hg_group),
        names_to = "metric",
        values_to = "value"
      ),
    branch_geo_diag |>
      dplyr::mutate(section = "branch_geo_diag") |>
      tidyr::pivot_longer(
        cols = -c(section, hg_group),
        names_to = "metric",
        values_to = "value"
      ),
    network_diag |>
      dplyr::mutate(section = "network_diag") |>
      tidyr::pivot_longer(
        cols = -c(section, rule, hg_group),
        names_to = "metric",
        values_to = "value"
      ) |>
      dplyr::mutate(metric = paste0(rule, "__", metric))
  )
  
  robustness_rows <- purrr::map_dfr(ROBUSTNESS_SPECS$perturbation, function(perturbation) {
    branch_tbl_rx <- make_robustness_branch_table(branch_tbl_main, perturbation)
    edges_rx <- build_all_rules_observed(branch_tbl_rx)
    metrics_rx <- compute_rule_metrics(edges_rx)
    diff_rx <- compute_diff_table(metrics_rx)
    
    diff_rx |>
      dplyr::mutate(
        perturbation = perturbation,
        sign_negative = cross_region_share_diff < 0
      ) |>
      dplyr::select(
        perturbation,
        rule,
        cross_region_share_diff,
        sign_negative
      )
  })
  
  baseline_sign <- robustness_rows |>
    dplyr::filter(perturbation == "baseline") |>
    dplyr::select(rule, baseline_sign_negative = sign_negative)
  
  table_p2_rx <- robustness_rows |>
    dplyr::left_join(baseline_sign, by = "rule") |>
    dplyr::mutate(
      sign_matches_baseline = sign_negative == baseline_sign_negative
    ) |>
    dplyr::arrange(
      factor(perturbation, levels = ROBUSTNESS_SPECS$perturbation),
      factor(rule, levels = RULE_LEVELS)
    )
  
  # Consistency checks
  
  assert_equal(age_ph908$n_branches[[1]], 79L, "PH908 branch count drift in Pillar 2.")
  assert_equal(age_ph908$total_mass[[1]], 623L, "PH908 total branch mass drift in Pillar 2.")
  assert_equal(age_ph908$n_relic[[1]], 50L, "PH908 relic count drift in Pillar 2.")
  assert_equal(age_ph908$n_middle[[1]], 20L, "PH908 middle count drift in Pillar 2.")
  assert_equal(age_ph908$n_founder[[1]], 9L, "PH908 founder count drift in Pillar 2.")
  
  assert_equal(age_r1a$n_branches[[1]], 81L, "R1a branch count drift in Pillar 2.")
  assert_equal(age_r1a$total_mass[[1]], 729L, "R1a total branch mass drift in Pillar 2.")
  assert_equal(age_r1a$n_relic[[1]], 38L, "R1a relic count drift in Pillar 2.")
  assert_equal(age_r1a$n_middle[[1]], 29L, "R1a middle count drift in Pillar 2.")
  assert_equal(age_r1a$n_founder[[1]], 14L, "R1a founder count drift in Pillar 2.")
  
  for (rr in RULE_LEVELS) {
    nd_ph <- network_diag |>
      dplyr::filter(rule == rr, hg_group == "PH908")
    nd_r1a <- network_diag |>
      dplyr::filter(rule == rr, hg_group == "R1a")
    
    assert_equal(nd_ph$n_edges_output[[1]], 78L, paste0("PH908 edge count drift for rule ", rr, "."))
    assert_equal(nd_r1a$n_edges_output[[1]], 80L, paste0("R1a edge count drift for rule ", rr, "."))
    assert_equal(nd_ph$n_unmatched_branches[[1]], 0L, paste0("PH908 unmatched branches drift for rule ", rr, "."))
    assert_equal(nd_r1a$n_unmatched_branches[[1]], 0L, paste0("R1a unmatched branches drift for rule ", rr, "."))
  }
  
  row_nearest <- get_main_row(figuredata_p2_main, "nearest_older")
  row_window <- get_main_row(figuredata_p2_main, "nearest_older_window")
  row_cost <- get_main_row(figuredata_p2_main, "distance_age_cost")
  
  assert_close(row_nearest$cross_region_share_PH908[[1]], 0.654, tol = 0.0015, "PH908 nearest_older cross-region share drift.")
  assert_close(row_nearest$cross_region_share_R1a[[1]], 0.850, tol = 0.0015, "R1a nearest_older cross-region share drift.")
  assert_close(row_nearest$cross_region_share_diff[[1]], -0.196, tol = 0.0015, "Nearest_older cross-region difference drift.")
  
  assert_close(row_window$cross_region_share_PH908[[1]], 0.654, tol = 0.0015, "PH908 nearest_older_window cross-region share drift.")
  assert_close(row_window$cross_region_share_R1a[[1]], 0.863, tol = 0.0015, "R1a nearest_older_window cross-region share drift.")
  assert_close(row_window$cross_region_share_diff[[1]], -0.209, tol = 0.0015, "Nearest_older_window cross-region difference drift.")
  
  assert_close(row_cost$cross_region_share_PH908[[1]], 0.718, tol = 0.0015, "PH908 distance_age_cost cross-region share drift.")
  assert_close(row_cost$cross_region_share_R1a[[1]], 0.875, tol = 0.0015, "R1a distance_age_cost cross-region share drift.")
  assert_close(row_cost$cross_region_share_diff[[1]], -0.157, tol = 0.0015, "Distance_age_cost cross-region difference drift.")
  
  assert_true(
    all(table_p2_rx$cross_region_share_diff < 0, na.rm = TRUE),
    "At least one influential-unit robustness perturbation flipped the sign of the cross-region contrast."
  )
  
  if (N_PERM_PRIMARY >= 5000L) {
    assert_close(row_nearest$p_value_one_sided_PH908_less[[1]], 0.0197, tol = 0.01, "Nearest_older primary p-value drift.")
    assert_close(row_window$p_value_one_sided_PH908_less[[1]], 0.0093, tol = 0.01, "Nearest_older_window primary p-value drift.")
    assert_close(row_cost$p_value_one_sided_PH908_less[[1]], 0.0020, tol = 0.01, "Distance_age_cost primary p-value drift.")
  }
  
  if (N_PERM_GEO_SENSITIVITY >= 5000L) {
    assert_close(row_nearest$p_value_one_sided_geo[[1]], 0.189, tol = 0.03, "Nearest_older geography-aware p-value drift.")
    assert_close(row_window$p_value_one_sided_geo[[1]], 0.128, tol = 0.03, "Nearest_older_window geography-aware p-value drift.")
    assert_close(row_cost$p_value_one_sided_geo[[1]], 0.056, tol = 0.03, "Distance_age_cost geography-aware p-value drift.")
  }
  
  # Facts
  
  n_sign_flips <- sum(!table_p2_rx$sign_matches_baseline, na.rm = TRUE)
  all_signs_negative <- all(table_p2_rx$cross_region_share_diff < 0, na.rm = TRUE)
  
  results_sentence <- paste0(
    "Across all three age-aware branch-link construction rules, PH908 showed lower inferred cross-region linkage than the R1a control, with observed PH908-R1a differences of ",
    sprintf("%.3f", row_nearest$cross_region_share_diff[[1]]), ", ",
    sprintf("%.3f", row_window$cross_region_share_diff[[1]]), ", and ",
    sprintf("%.3f", row_cost$cross_region_share_diff[[1]]),
    " under nearest_older, nearest_older_window, and distance_age_cost, respectively."
  )
  
  caption_sentence <- paste0(
    "Cross-region linkage share remained lower in PH908 than in R1a under all three branch-link construction rules; observed PH908-R1a differences were ",
    sprintf("%.3f", row_nearest$cross_region_share_diff[[1]]), ", ",
    sprintf("%.3f", row_window$cross_region_share_diff[[1]]), ", and ",
    sprintf("%.3f", row_cost$cross_region_share_diff[[1]]), "."
  )
  
  sensitivity_sentence <- paste0(
    "Influential-unit robustness preserved the negative sign under all perturbations (",
    paste(ROBUSTNESS_SPECS$perturbation, collapse = ", "),
    "), with ", n_sign_flips,
    " sign mismatches relative to baseline across all rule-by-perturbation checks."
  )
  
  writeLines(
    c(
      "script: 03_pillar2_topology.R",
      paste0("timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0("seed: ", SEED_MAIN),
      paste0("config_file: ", path_for_facts(CONFIG_FILE)),
      "project_root: .",
      paste0("git_commit: ", git_commit),
      "inputs:",
      paste0("  - ", path_for_facts(FILE_PREPROCESSED_RDS)),
      paste0("  - ", path_for_facts(FILE_BRANCH_AGES)),
      "",
      paste0("ph908_gating_mode: ", PH908_GATING_MODE),
      paste0("max_age_gap_window: ", MAX_AGE_GAP_WINDOW),
      paste0("lambda_age_cost: ", LAMBDA_AGE_COST),
      paste0("n_perm_primary: ", N_PERM_PRIMARY),
      paste0("n_perm_geo_sensitivity: ", N_PERM_GEO_SENSITIVITY),
      paste0("age_bin_n: ", AGE_BIN_N),
      paste0("size_bin_n: ", SIZE_BIN_N),
      paste0("geo_lat_bin_n: ", GEO_LAT_BIN_N),
      paste0("geo_lon_bin_n: ", GEO_LON_BIN_N),
      "",
      "filter_diagnostics:",
      paste0("  raw_rows: ", filter_diag$n[filter_diag$step == "raw_rows"]),
      paste0("  matched_true: ", filter_diag$n[filter_diag$step == "matched_true"]),
      paste0("  exclude_geo_primary_false: ", filter_diag$n[filter_diag$step == "exclude_geo_primary_false"]),
      paste0("  is_balkan_country_true: ", filter_diag$n[filter_diag$step == "is_balkan_country_true"]),
      paste0("  terminal_snp_present: ", filter_diag$n[filter_diag$step == "terminal_snp_present"]),
      paste0("  region_present: ", filter_diag$n[filter_diag$step == "region_present"]),
      paste0("  finite_lat_long: ", filter_diag$n[filter_diag$step == "finite_lat_long"]),
      paste0("  usable_sample_rows: ", nrow(base_rows)),
      "",
      "branch_structure:",
      paste0("  PH908_branches: ", age_ph908$n_branches[[1]]),
      paste0("  PH908_total_mass: ", age_ph908$total_mass[[1]]),
      paste0("  PH908_relic_middle_founder: ", age_ph908$n_relic[[1]], "/", age_ph908$n_middle[[1]], "/", age_ph908$n_founder[[1]]),
      paste0("  R1a_branches: ", age_r1a$n_branches[[1]]),
      paste0("  R1a_total_mass: ", age_r1a$total_mass[[1]]),
      paste0("  R1a_relic_middle_founder: ", age_r1a$n_relic[[1]], "/", age_r1a$n_middle[[1]], "/", age_r1a$n_founder[[1]]),
      "",
      "network_edges:",
      paste0("  nearest_older_edges_PH908_R1a: ",
             network_diag$n_edges_output[network_diag$rule == "nearest_older" & network_diag$hg_group == "PH908"], "/",
             network_diag$n_edges_output[network_diag$rule == "nearest_older" & network_diag$hg_group == "R1a"]),
      paste0("  nearest_older_window_edges_PH908_R1a: ",
             network_diag$n_edges_output[network_diag$rule == "nearest_older_window" & network_diag$hg_group == "PH908"], "/",
             network_diag$n_edges_output[network_diag$rule == "nearest_older_window" & network_diag$hg_group == "R1a"]),
      paste0("  distance_age_cost_edges_PH908_R1a: ",
             network_diag$n_edges_output[network_diag$rule == "distance_age_cost" & network_diag$hg_group == "PH908"], "/",
             network_diag$n_edges_output[network_diag$rule == "distance_age_cost" & network_diag$hg_group == "R1a"]),
      "",
      "observed_cross_region_summary:",
      paste0("  nearest_older_PH908_R1a_diff: ",
             sprintf("%.3f", row_nearest$cross_region_share_PH908[[1]]), " / ",
             sprintf("%.3f", row_nearest$cross_region_share_R1a[[1]]), " / ",
             sprintf("%.3f", row_nearest$cross_region_share_diff[[1]])),
      paste0("  nearest_older_window_PH908_R1a_diff: ",
             sprintf("%.3f", row_window$cross_region_share_PH908[[1]]), " / ",
             sprintf("%.3f", row_window$cross_region_share_R1a[[1]]), " / ",
             sprintf("%.3f", row_window$cross_region_share_diff[[1]])),
      paste0("  distance_age_cost_PH908_R1a_diff: ",
             sprintf("%.3f", row_cost$cross_region_share_PH908[[1]]), " / ",
             sprintf("%.3f", row_cost$cross_region_share_R1a[[1]]), " / ",
             sprintf("%.3f", row_cost$cross_region_share_diff[[1]])),
      "",
      "primary_null_pvalues:",
      paste0("  nearest_older: ", format(row_nearest$p_value_one_sided_PH908_less[[1]], digits = 8)),
      paste0("  nearest_older_window: ", format(row_window$p_value_one_sided_PH908_less[[1]], digits = 8)),
      paste0("  distance_age_cost: ", format(row_cost$p_value_one_sided_PH908_less[[1]], digits = 8)),
      "",
      "geo_sensitivity_pvalues:",
      paste0("  nearest_older: ", format(row_nearest$p_value_one_sided_geo[[1]], digits = 8)),
      paste0("  nearest_older_window: ", format(row_window$p_value_one_sided_geo[[1]], digits = 8)),
      paste0("  distance_age_cost: ", format(row_cost$p_value_one_sided_geo[[1]], digits = 8)),
      "",
      paste0("results_sentence: ", results_sentence),
      paste0("caption_sentence: ", caption_sentence),
      paste0("sensitivity_sentence: ", sensitivity_sentence),
      paste0("sign_flip_anywhere: ", if (all_signs_negative) "no" else "yes"),
      "note: p-value consistency checks are conditional on the configured permutation depth."
    ),
    FILE_P2_FACTS
  )
  
  # Output serialization
  
  readr::write_csv(figuredata_p2_main, FILE_P2_FIG_MAIN)
  readr::write_csv(figuredata_p2_null, FILE_P2_FIG_NULL)
  readr::write_csv(table_p2_supp, FILE_P2_SUPP_TABLE)
  readr::write_csv(table_p2_rx, FILE_P2_RX_TABLE)
  
  results_p2 <- list(
    metadata = list(
      script = "03_pillar2_topology.R",
      seed = SEED_MAIN,
      config_file = path_for_facts(CONFIG_FILE),
      inputs = c(
        path_for_facts(FILE_PREPROCESSED_RDS),
        path_for_facts(FILE_BRANCH_AGES)
      ),
      ph908_gating_mode = PH908_GATING_MODE,
      max_age_gap_window = MAX_AGE_GAP_WINDOW,
      lambda_age_cost = LAMBDA_AGE_COST,
      age_bin_n = AGE_BIN_N,
      size_bin_n = SIZE_BIN_N,
      geo_lat_bin_n = GEO_LAT_BIN_N,
      geo_lon_bin_n = GEO_LON_BIN_N,
      n_perm_primary = N_PERM_PRIMARY,
      n_perm_geo_sensitivity = N_PERM_GEO_SENSITIVITY,
      git_commit = git_commit
    ),
    filter_diag = filter_diag,
    age_diag = age_diag,
    branch_geo_diag = branch_geo_diag,
    branch_table = branch_tbl_main,
    edges_observed = edges_obs,
    metrics_observed = metrics_obs,
    obs_diff = obs_diff_tbl,
    network_diag = network_diag,
    strata_diag = strata_diag,
    strata_geo_diag = strata_geo_diag,
    null_primary = null_tbl,
    null_geo_sensitivity = null_geo_tbl,
    tests = tests_tbl,
    figuredata_main = figuredata_p2_main,
    figuredata_null = figuredata_p2_null,
    table_supp = table_p2_supp,
    table_rx = table_p2_rx
  )
  
  saveRDS(results_p2, FILE_P2_RESULTS)
  
}, warning = function(w) {
  warning_log <<- c(warning_log, conditionMessage(w))
  invokeRestart("muffleWarning")
})

warning_log_path <- file.path(DIR_FACTS, "warnings_p2.txt")
write_warning_log(warning_log, warning_log_path)

cat("============================================================\n")
cat("PILLAR 2 COMPLETE\n")
cat("============================================================\n")
cat("Results object: ", FILE_P2_RESULTS, "\n", sep = "")
cat("Facts file: ", FILE_P2_FACTS, "\n", sep = "")
cat("Figuredata main: ", FILE_P2_FIG_MAIN, "\n", sep = "")
cat("Figuredata null: ", FILE_P2_FIG_NULL, "\n", sep = "")
cat("Supplement table: ", FILE_P2_SUPP_TABLE, "\n", sep = "")
cat("Robustness table: ", FILE_P2_RX_TABLE, "\n", sep = "")
cat("Warnings log: ", warning_log_path, "\n", sep = "")
cat("============================================================\n")