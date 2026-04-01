# 04_pillar3_anchor.R
# Anchor-packet occupancy, layered retention, and relic age-structured persistence.

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

FILE_P3_RESULTS <- file.path(DIR_RESULTS, "results_p3.rds")
FILE_P3_FIG_MAIN <- file.path(DIR_FIGUREDATA, "figuredata_p3_main.csv")
FILE_P3_FIG_SENS <- file.path(DIR_FIGUREDATA, "figuredata_p3_sens.csv")
FILE_P3_FIG_SUPP <- file.path(DIR_FIGUREDATA, "figuredata_p3_supp.csv")
FILE_P3_SUPP_TABLE <- file.path(DIR_TABLES, "table_p3_supp.csv")
FILE_P3_RX_TABLE <- file.path(DIR_TABLES, "table_p3_rx.csv")
FILE_P3_FACTS <- file.path(DIR_FACTS, "facts_p3.txt")

# Analysis settings

PH908_GATING_MODE <- "primary"

RELIC_MAX_N <- 3L
FOUNDER_MIN_N <- 15L

AGE_BIN_N <- 3L
SIZE_BIN_N <- 2L

# Full manuscript-grade permutation settings.
N_PERM_MAIN <- 5000L
N_PERM_RELIC_AGE <- 5000L

ANCHOR_TOP_K_VALUES <- c(4L, 5L)
ANCHOR_CUM_VALUES <- c(0.50, 0.60, 0.70)
PRIMARY_ANCHOR_ID <- "cum60"

MIN_BRANCHES_PER_GROUP_PRIMARY <- 8L
MIN_BRANCHES_PER_GROUP_MINIMAL <- 5L
EPSILON <- 0.5

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
cat("04 PILLAR 3 — ANCHOR-PACKET OCCUPANCY\n")
cat("============================================================\n")
cat("Root: ", DIR_ROOT, "\n", sep = "")
cat("Seed: ", SEED_MAIN, "\n", sep = "")
cat("Inputs:\n")
cat("  - ", FILE_PREPROCESSED_RDS, "\n", sep = "")
cat("  - ", FILE_BRANCH_AGES, "\n", sep = "")
cat("Primary permutations: ", N_PERM_MAIN, "\n", sep = "")
cat("Relic-age permutations: ", N_PERM_RELIC_AGE, "\n\n", sep = "")

# Helpers

path_for_facts <- function(x) {
  x <- as.character(x)
  if (!length(x) || is.na(x) || x == "") return(x)
  x_norm <- normalizePath(x, winslash = "/", mustWork = FALSE)
  root_norm <- normalizePath(DIR_ROOT, winslash = "/", mustWork = FALSE)
  root_esc <- gsub("([][{}()+*^$|\\?.])", "\\\\\\1", root_norm)
  sub(paste0("^", root_esc), ".", x_norm)
}

append_facts_lines <- function(path, lines) {
  cat(paste(lines, collapse = "\n"), "\n", file = path, append = TRUE, sep = "")
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

writeLines(
  c(
    "script: 04_pillar3_anchor.R",
    paste0("timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("seed: ", SEED_MAIN),
    paste0("config_file: ", path_for_facts(CONFIG_FILE)),
    "project_root: .",
    paste0("git_commit: ", git_commit),
    "inputs:",
    paste0("  - ", path_for_facts(FILE_PREPROCESSED_RDS)),
    paste0("  - ", path_for_facts(FILE_BRANCH_AGES)),
    ""
  ),
  FILE_P3_FACTS,
  useBytes = TRUE
)

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

one_sided_greater_perm_p <- function(null_vals, obs) {
  null_vals <- null_vals[is.finite(null_vals)]
  if (!is.finite(obs) || length(null_vals) == 0) return(NA_real_)
  (1 + sum(null_vals >= obs)) / (length(null_vals) + 1)
}

make_quantile_bins <- function(x, n_bins = 3L, prefix = "bin") {
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

log2_odds_ratio <- function(a, b, c, d, eps = EPSILON) {
  as.numeric(log2(((a + eps) / (b + eps)) / ((c + eps) / (d + eps))))
}

panel_support_label <- function(n1, n2) {
  dplyr::case_when(
    n1 >= MIN_BRANCHES_PER_GROUP_PRIMARY & n2 >= MIN_BRANCHES_PER_GROUP_PRIMARY ~ "primary_ready",
    n1 >= MIN_BRANCHES_PER_GROUP_MINIMAL & n2 >= MIN_BRANCHES_PER_GROUP_MINIMAL ~ "minimum_formal_counts",
    TRUE ~ "underpowered"
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

standardize_age_table <- function(tbl) {
  if (!all(c("terminal_snp", "tmrca_ybp") %in% names(tbl))) {
    stop("branch_ages.csv must contain terminal_snp and tmrca_ybp.", call. = FALSE)
  }
  
  tbl |>
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
}

make_filter_diag <- function(cl_raw) {
  tibble::tibble(
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
}

prepare_input_rows <- function(cl_raw) {
  if (!"is_balkan_country" %in% names(cl_raw)) {
    if (!"country_code" %in% names(cl_raw)) {
      stop(
        "Dataset is missing both is_balkan_country and country_code, so the Pillar 3 Balkan-country filter cannot be reconstructed.",
        call. = FALSE
      )
    }
    cl_raw <- cl_raw |>
      dplyr::mutate(is_balkan_country = toupper(country_code) %in% BALKAN_ISO2)
    message("Reconstructed is_balkan_country from country_code for Pillar 3 filtering.")
  }
  
  required_cols <- c(
    "matched", "exclude_geo_primary", "is_balkan_country",
    "terminal_snp", "Region", "lat", "long",
    "major_hg", "tmrca", "is_ph908_primary"
  )
  missing_cols <- setdiff(required_cols, names(cl_raw))
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Dataset missing required columns after Pillar 3 input preparation: ",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  cl_raw |>
    dplyr::filter(
      matched %in% TRUE,
      exclude_geo_primary %in% FALSE,
      is_balkan_country %in% TRUE,
      !is.na(terminal_snp),
      terminal_snp != "",
      !is.na(Region),
      Region != "",
      is.finite(lat),
      is.finite(long)
    ) |>
    dplyr::mutate(
      terminal_snp_norm = normalize_snp(terminal_snp),
      hg_group = dplyr::case_when(
        is_ph908_primary %in% TRUE ~ "PH908",
        major_hg == "R1a" ~ "R1a",
        TRUE ~ "Other"
      )
    ) |>
    dplyr::filter(hg_group %in% c("PH908", "R1a"))
}

build_branch_table <- function(base_rows, age_tbl) {
  branch_tbl <- base_rows |>
    dplyr::group_by(hg_group, terminal_snp_norm) |>
    dplyr::summarise(
      branch_size = dplyr::n(),
      lat_c = mean(lat, na.rm = TRUE),
      lon_c = mean(long, na.rm = TRUE),
      dominant_region = dominant_region_fun(Region),
      n_regions_in_branch = dplyr::n_distinct(Region[!is.na(Region) & Region != ""]),
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
    dplyr::filter(is.finite(tmrca_final))
  
  if (nrow(branch_tbl) == 0) {
    stop("No branches with usable TMRCA after merging age information.", call. = FALSE)
  }
  
  branch_tbl
}

make_age_diag <- function(branch_tbl) {
  branch_tbl |>
    dplyr::group_by(hg_group) |>
    dplyr::summarise(
      n_branches = dplyr::n(),
      total_mass = sum(branch_size, na.rm = TRUE),
      n_relic = sum(class == "Relic", na.rm = TRUE),
      n_middle = sum(class == "Middle", na.rm = TRUE),
      n_founder = sum(class == "Founder", na.rm = TRUE),
      age_from_age_file = sum(age_source == "branch_ages_csv", na.rm = TRUE),
      age_from_dataset = sum(age_source == "dataset_tmrca", na.rm = TRUE),
      pct_age_from_age_file = 100 * mean(age_source == "branch_ages_csv", na.rm = TRUE),
      pct_age_from_dataset = 100 * mean(age_source == "dataset_tmrca", na.rm = TRUE),
      median_branch_size = stats::median(branch_size, na.rm = TRUE),
      max_branch_size = max(branch_size, na.rm = TRUE),
      .groups = "drop"
    )
}

build_anchor_definition_tbl <- function(ph908_relic_rank_tbl) {
  topk_tbl <- purrr::map_dfr(ANCHOR_TOP_K_VALUES, function(k) {
    reg_tbl <- ph908_relic_rank_tbl |>
      dplyr::filter(region_rank <= k)
    
    tibble::tibble(
      anchor_id = paste0("top", k),
      anchor_family = "PH908_relic_rank",
      anchor_type = "top_k",
      anchor_param = as.character(k),
      anchor_region_count = nrow(reg_tbl),
      anchor_branch_share_PH908_source = sum(reg_tbl$region_share, na.rm = TRUE),
      anchor_regions_concat = paste(reg_tbl$dominant_region, collapse = "|")
    )
  })
  
  cum_tbl <- purrr::map_dfr(ANCHOR_CUM_VALUES, function(cutoff) {
    reg_tbl <- ph908_relic_rank_tbl |>
      dplyr::filter(cum_share <= cutoff)
    
    if (nrow(reg_tbl) == 0) {
      reg_tbl <- ph908_relic_rank_tbl |>
        dplyr::slice_head(n = 1)
    }
    
    if (max(reg_tbl$cum_share, na.rm = TRUE) < cutoff) {
      extra_row <- ph908_relic_rank_tbl |>
        dplyr::filter(region_rank == max(reg_tbl$region_rank) + 1L) |>
        dplyr::slice_head(n = 1)
      
      if (nrow(extra_row) > 0) {
        reg_tbl <- dplyr::bind_rows(reg_tbl, extra_row)
      }
    }
    
    tibble::tibble(
      anchor_id = paste0("cum", sprintf("%02d", as.integer(round(cutoff * 100)))),
      anchor_family = "PH908_relic_rank",
      anchor_type = "cum_share",
      anchor_param = sprintf("%.2f", cutoff),
      anchor_region_count = nrow(reg_tbl),
      anchor_branch_share_PH908_source = sum(reg_tbl$region_share, na.rm = TRUE),
      anchor_regions_concat = paste(reg_tbl$dominant_region, collapse = "|")
    )
  })
  
  dplyr::bind_rows(topk_tbl, cum_tbl) |>
    dplyr::distinct(anchor_id, .keep_all = TRUE) |>
    dplyr::mutate(
      anchor_label = dplyr::case_when(
        anchor_id == PRIMARY_ANCHOR_ID ~ paste0("Primary: ", anchor_id),
        TRUE ~ paste0("Sensitivity: ", anchor_id)
      )
    ) |>
    dplyr::arrange(
      dplyr::desc(anchor_id == PRIMARY_ANCHOR_ID),
      anchor_type,
      anchor_param
    )
}

extract_anchor_regions <- function(anchor_def_tbl, anchor_id) {
  x <- anchor_def_tbl |>
    dplyr::filter(anchor_id == !!anchor_id) |>
    dplyr::pull(anchor_regions_concat)
  
  if (length(x) == 0 || is.na(x[1]) || x[1] == "") return(character())
  strsplit(x[1], "\\|")[[1]]
}

build_panel_branch_tbl <- function(branch_tbl) {
  dplyr::bind_rows(
    branch_tbl |>
      dplyr::filter(class == "Relic") |>
      dplyr::mutate(panel = "relic"),
    branch_tbl |>
      dplyr::filter(class == "Middle") |>
      dplyr::mutate(panel = "middle"),
    branch_tbl |>
      dplyr::filter(class == "Founder") |>
      dplyr::mutate(panel = "founder"),
    branch_tbl |>
      dplyr::mutate(panel = "all")
  ) |>
    dplyr::mutate(
      panel = factor(panel, levels = c("relic", "middle", "all", "founder"))
    )
}

compute_anchor_metrics <- function(df_panel_anchor, focal_group = "PH908", control_group = "R1a") {
  counts <- df_panel_anchor |>
    dplyr::count(hg_group, inside_anchor, name = "n") |>
    tidyr::complete(
      hg_group = c(focal_group, control_group),
      inside_anchor = c(TRUE, FALSE),
      fill = list(n = 0)
    )
  
  a <- counts |>
    dplyr::filter(hg_group == focal_group, inside_anchor %in% TRUE) |>
    dplyr::pull(n)
  b <- counts |>
    dplyr::filter(hg_group == focal_group, inside_anchor %in% FALSE) |>
    dplyr::pull(n)
  c <- counts |>
    dplyr::filter(hg_group == control_group, inside_anchor %in% TRUE) |>
    dplyr::pull(n)
  d <- counts |>
    dplyr::filter(hg_group == control_group, inside_anchor %in% FALSE) |>
    dplyr::pull(n)
  
  anchor_share_focal <- if ((a + b) > 0) a / (a + b) else NA_real_
  anchor_share_control <- if ((c + d) > 0) c / (c + d) else NA_real_
  
  tibble::tibble(
    focal_group = focal_group,
    control_group = control_group,
    focal_n = a + b,
    control_n = c + d,
    focal_anchor_n = a,
    focal_non_anchor_n = b,
    control_anchor_n = c,
    control_non_anchor_n = d,
    anchor_share_focal = anchor_share_focal,
    anchor_share_control = anchor_share_control,
    anchor_share_diff = anchor_share_focal - anchor_share_control,
    anchor_log2_odds_ratio = log2_odds_ratio(a, b, c, d, eps = EPSILON)
  )
}

pretty_age_mode <- function(x) {
  dplyr::case_when(
    x == "all_relics" ~ "All relics",
    x == "old_vs_rest_median" ~ "Old vs rest (median split)",
    x == "within_lineage_tertiles" ~ "Within-lineage tertiles",
    TRUE ~ as.character(x)
  )
}

pretty_age_stratum <- function(x) {
  dplyr::case_when(
    x == "all_relics" ~ "All relics",
    x == "old" ~ "Old",
    x == "rest" ~ "Rest",
    x == "oldest" ~ "Oldest tertile",
    x == "middle" ~ "Middle tertile",
    x == "youngest" ~ "Youngest tertile",
    TRUE ~ as.character(x)
  )
}

build_age_labeled_relics <- function(df_relic) {
  all_tbl <- df_relic |>
    dplyr::mutate(
      age_mode = "all_relics",
      age_stratum = "all_relics"
    )
  
  med_tbl <- df_relic |>
    dplyr::group_by(hg_group) |>
    dplyr::mutate(
      med_tmrca = stats::median(tmrca_final, na.rm = TRUE),
      age_mode = "old_vs_rest_median",
      age_stratum = dplyr::case_when(
        tmrca_final >= med_tmrca ~ "old",
        TRUE ~ "rest"
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-med_tmrca)
  
  tert_tbl <- df_relic |>
    dplyr::group_by(hg_group) |>
    dplyr::mutate(
      q1 = stats::quantile(tmrca_final, probs = 1 / 3, na.rm = TRUE, type = 7),
      q2 = stats::quantile(tmrca_final, probs = 2 / 3, na.rm = TRUE, type = 7),
      age_mode = "within_lineage_tertiles",
      age_stratum = dplyr::case_when(
        tmrca_final >= q2 ~ "oldest",
        tmrca_final >= q1 & tmrca_final < q2 ~ "middle",
        TRUE ~ "youngest"
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-q1, -q2)
  
  dplyr::bind_rows(all_tbl, med_tbl, tert_tbl) |>
    dplyr::mutate(
      age_mode_pretty = pretty_age_mode(age_mode),
      age_stratum_pretty = pretty_age_stratum(age_stratum)
    )
}

compute_age_mode_anchor_metrics <- function(df, focal_group = "PH908", control_group = "R1a") {
  compute_anchor_metrics(df, focal_group = focal_group, control_group = control_group)
}

make_panelized_observed <- function(panel_branch_tbl, anchor_def_tbl) {
  observed_list <- list()
  packet_comp_list <- list()
  k <- 0L
  j <- 0L
  
  for (aid in anchor_def_tbl$anchor_id) {
    anchor_regions <- extract_anchor_regions(anchor_def_tbl, aid)
    
    for (pp in levels(panel_branch_tbl$panel)) {
      df_pa <- panel_branch_tbl |>
        dplyr::filter(panel == pp) |>
        dplyr::mutate(
          anchor_id = aid,
          inside_anchor = dominant_region %in% anchor_regions
        )
      
      met <- compute_anchor_metrics(df_pa, focal_group = "PH908", control_group = "R1a") |>
        dplyr::mutate(
          panel = pp,
          anchor_id = aid,
          panel_support = panel_support_label(focal_n[[1]], control_n[[1]])
        )
      
      k <- k + 1L
      observed_list[[k]] <- met
      
      comp_tbl <- df_pa |>
        dplyr::count(anchor_id, panel, hg_group, dominant_region, inside_anchor, name = "branch_n") |>
        dplyr::group_by(anchor_id, panel, hg_group) |>
        dplyr::mutate(branch_share = branch_n / sum(branch_n)) |>
        dplyr::ungroup()
      
      j <- j + 1L
      packet_comp_list[[j]] <- comp_tbl
    }
  }
  
  observed_tbl <- dplyr::bind_rows(observed_list) |>
    dplyr::left_join(
      anchor_def_tbl |>
        dplyr::select(anchor_id, anchor_type, anchor_param, anchor_region_count, anchor_regions_concat),
      by = "anchor_id"
    ) |>
    dplyr::mutate(
      panel = factor(panel, levels = c("relic", "middle", "all", "founder"))
    ) |>
    dplyr::arrange(anchor_id, panel)
  
  packet_comp_tbl <- dplyr::bind_rows(packet_comp_list)
  
  list(observed_tbl = observed_tbl, packet_comp_tbl = packet_comp_tbl)
}

make_strata_diag <- function(panel_branch_tbl) {
  panel_branch_tbl |>
    dplyr::mutate(
      age_bin = make_quantile_bins(tmrca_final, n_bins = AGE_BIN_N, prefix = "age"),
      size_bin = make_quantile_bins(log(branch_size + 1), n_bins = SIZE_BIN_N, prefix = "size"),
      strata = dplyr::case_when(
        panel %in% c("relic", "middle", "founder") ~ as.character(interaction(age_bin, size_bin, drop = TRUE)),
        panel == "all" ~ as.character(interaction(age_bin, size_bin, class, drop = TRUE)),
        TRUE ~ as.character(interaction(age_bin, size_bin, drop = TRUE))
      )
    ) |>
    dplyr::count(panel, strata, hg_group, name = "n_branches") |>
    dplyr::arrange(panel, strata, hg_group)
}

run_main_permutation_null <- function(panel_branch_tbl, anchor_def_tbl, n_perm) {
  if (n_perm < 1) return(tibble::tibble())
  
  null_list <- vector("list", n_perm)
  
  for (b in seq_len(n_perm)) {
    perm_res_list <- list()
    kk <- 0L
    
    for (pp in levels(panel_branch_tbl$panel)) {
      perm_base <- panel_branch_tbl |>
        dplyr::filter(panel == pp) |>
        dplyr::mutate(
          age_bin = make_quantile_bins(tmrca_final, n_bins = AGE_BIN_N, prefix = "age"),
          size_bin = make_quantile_bins(log(branch_size + 1), n_bins = SIZE_BIN_N, prefix = "size"),
          strata = dplyr::case_when(
            panel %in% c("relic", "middle", "founder") ~ as.character(interaction(age_bin, size_bin, drop = TRUE)),
            panel == "all" ~ as.character(interaction(age_bin, size_bin, class, drop = TRUE)),
            TRUE ~ as.character(interaction(age_bin, size_bin, drop = TRUE))
          )
        )
      
      perm_group <- permute_within_strata(perm_base$hg_group, perm_base$strata)
      
      perm_df <- perm_base |>
        dplyr::mutate(hg_group = perm_group)
      
      for (aid in anchor_def_tbl$anchor_id) {
        anchor_regions <- extract_anchor_regions(anchor_def_tbl, aid)
        
        perm_df2 <- perm_df |>
          dplyr::mutate(inside_anchor = dominant_region %in% anchor_regions)
        
        met <- compute_anchor_metrics(perm_df2, focal_group = "PH908", control_group = "R1a") |>
          dplyr::mutate(
            perm_id = b,
            panel = pp,
            anchor_id = aid
          )
        
        kk <- kk + 1L
        perm_res_list[[kk]] <- met
      }
    }
    
    null_list[[b]] <- dplyr::bind_rows(perm_res_list)
  }
  
  dplyr::bind_rows(null_list) |>
    dplyr::left_join(
      anchor_def_tbl |>
        dplyr::select(anchor_id, anchor_type, anchor_param),
      by = "anchor_id"
    )
}

build_formal_tests <- function(observed_tbl, null_tbl) {
  observed_tbl |>
    dplyr::rowwise() |>
    dplyr::mutate(
      ci_anchor_share_diff_low = ci95(
        null_tbl$anchor_share_diff[
          null_tbl$panel == panel &
            null_tbl$anchor_id == anchor_id
        ]
      )[1],
      ci_anchor_share_diff_high = ci95(
        null_tbl$anchor_share_diff[
          null_tbl$panel == panel &
            null_tbl$anchor_id == anchor_id
        ]
      )[2],
      p_one_sided_anchor_share_diff = one_sided_greater_perm_p(
        null_tbl$anchor_share_diff[
          null_tbl$panel == panel &
            null_tbl$anchor_id == anchor_id
        ],
        anchor_share_diff
      ),
      p_two_sided_anchor_share_diff = two_sided_perm_p(
        null_tbl$anchor_share_diff[
          null_tbl$panel == panel &
            null_tbl$anchor_id == anchor_id
        ],
        anchor_share_diff
      ),
      ci_log2_or_low = ci95(
        null_tbl$anchor_log2_odds_ratio[
          null_tbl$panel == panel &
            null_tbl$anchor_id == anchor_id
        ]
      )[1],
      ci_log2_or_high = ci95(
        null_tbl$anchor_log2_odds_ratio[
          null_tbl$panel == panel &
            null_tbl$anchor_id == anchor_id
        ]
      )[2],
      p_one_sided_log2_or = one_sided_greater_perm_p(
        null_tbl$anchor_log2_odds_ratio[
          null_tbl$panel == panel &
            null_tbl$anchor_id == anchor_id
        ],
        anchor_log2_odds_ratio
      ),
      p_two_sided_log2_or = two_sided_perm_p(
        null_tbl$anchor_log2_odds_ratio[
          null_tbl$panel == panel &
            null_tbl$anchor_id == anchor_id
        ],
        anchor_log2_odds_ratio
      ),
      support_label = dplyr::case_when(
        panel == "relic" &
          anchor_id == PRIMARY_ANCHOR_ID &
          panel_support == "primary_ready" &
          is.finite(anchor_share_diff) &
          anchor_share_diff > 0 &
          is.finite(p_one_sided_anchor_share_diff) &
          p_one_sided_anchor_share_diff < 0.05 ~ "Primary supported",
        panel %in% c("middle", "all") &
          panel_support != "underpowered" &
          is.finite(anchor_share_diff) &
          anchor_share_diff > 0 &
          is.finite(p_one_sided_anchor_share_diff) &
          p_one_sided_anchor_share_diff < 0.05 ~ "Secondary supported",
        panel == "founder" ~ "Descriptive only",
        panel_support == "underpowered" ~ "Underpowered",
        TRUE ~ "Not supported"
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(
      dplyr::desc(anchor_id == PRIMARY_ANCHOR_ID),
      factor(panel, levels = c("relic", "middle", "all", "founder")),
      anchor_id
    )
}

make_final_summary <- function(formal_tests) {
  formal_tests |>
    dplyr::select(
      panel,
      anchor_id,
      anchor_type,
      anchor_param,
      anchor_region_count,
      focal_group,
      control_group,
      focal_n,
      control_n,
      focal_anchor_n,
      control_anchor_n,
      anchor_share_focal,
      anchor_share_control,
      anchor_share_diff,
      ci_anchor_share_diff_low,
      ci_anchor_share_diff_high,
      p_one_sided_anchor_share_diff,
      p_two_sided_anchor_share_diff,
      anchor_log2_odds_ratio,
      ci_log2_or_low,
      ci_log2_or_high,
      p_one_sided_log2_or,
      p_two_sided_log2_or,
      panel_support,
      support_label,
      anchor_regions_concat
    ) |>
    dplyr::arrange(
      dplyr::desc(anchor_id == PRIMARY_ANCHOR_ID),
      factor(panel, levels = c("relic", "middle", "all", "founder")),
      anchor_id
    )
}

run_relic_age_extension <- function(branch_tbl, primary_anchor_regions, n_perm) {
  relic_tbl <- branch_tbl |>
    dplyr::filter(class == "Relic") |>
    dplyr::mutate(inside_anchor = dominant_region %in% primary_anchor_regions)
  
  relic_age_long <- build_age_labeled_relics(relic_tbl)
  
  relic_age_counts_tbl <- relic_age_long |>
    dplyr::count(age_mode, age_stratum, age_mode_pretty, age_stratum_pretty, hg_group, name = "n_branches") |>
    dplyr::arrange(age_mode, age_stratum, hg_group)
  
  relic_age_obs_list <- list()
  rr <- 0L
  
  for (am in unique(relic_age_long$age_mode)) {
    strata_here <- unique(relic_age_long$age_stratum[relic_age_long$age_mode == am])
    
    for (ss in strata_here) {
      df_sub <- relic_age_long |>
        dplyr::filter(age_mode == am, age_stratum == ss)
      
      met <- compute_age_mode_anchor_metrics(df_sub, focal_group = "PH908", control_group = "R1a") |>
        dplyr::mutate(
          age_mode = am,
          age_stratum = ss,
          age_mode_pretty = pretty_age_mode(am),
          age_stratum_pretty = pretty_age_stratum(ss),
          panel_support = panel_support_label(focal_n[[1]], control_n[[1]])
        )
      
      rr <- rr + 1L
      relic_age_obs_list[[rr]] <- met
    }
  }
  
  relic_age_obs_tbl <- dplyr::bind_rows(relic_age_obs_list) |>
    dplyr::arrange(age_mode, age_stratum)
  
  relic_age_null_tbl <- tibble::tibble()
  
  if (n_perm >= 1) {
    relic_age_null_list <- vector("list", n_perm)
    
    for (b in seq_len(n_perm)) {
      perm_base_relic <- relic_tbl |>
        dplyr::mutate(
          age_bin = make_quantile_bins(tmrca_final, n_bins = AGE_BIN_N, prefix = "age"),
          size_bin = make_quantile_bins(log(branch_size + 1), n_bins = SIZE_BIN_N, prefix = "size"),
          strata = as.character(interaction(age_bin, size_bin, drop = TRUE))
        )
      
      perm_group <- permute_within_strata(perm_base_relic$hg_group, perm_base_relic$strata)
      
      perm_relic <- perm_base_relic |>
        dplyr::mutate(hg_group = perm_group) |>
        dplyr::select(-age_bin, -size_bin, -strata)
      
      perm_relic_age_long <- build_age_labeled_relics(perm_relic)
      
      perm_sub_list <- list()
      qq <- 0L
      
      for (am in unique(perm_relic_age_long$age_mode)) {
        strata_here <- unique(perm_relic_age_long$age_stratum[perm_relic_age_long$age_mode == am])
        
        for (ss in strata_here) {
          df_sub <- perm_relic_age_long |>
            dplyr::filter(age_mode == am, age_stratum == ss)
          
          met <- compute_age_mode_anchor_metrics(df_sub, focal_group = "PH908", control_group = "R1a") |>
            dplyr::mutate(
              perm_id = b,
              age_mode = am,
              age_stratum = ss
            )
          
          qq <- qq + 1L
          perm_sub_list[[qq]] <- met
        }
      }
      
      relic_age_null_list[[b]] <- dplyr::bind_rows(perm_sub_list)
    }
    
    relic_age_null_tbl <- dplyr::bind_rows(relic_age_null_list) |>
      dplyr::mutate(
        age_mode_pretty = pretty_age_mode(age_mode),
        age_stratum_pretty = pretty_age_stratum(age_stratum)
      )
  }
  
  relic_age_tests_tbl <- relic_age_obs_tbl |>
    dplyr::rowwise() |>
    dplyr::mutate(
      ci_anchor_share_diff_low = ci95(
        relic_age_null_tbl$anchor_share_diff[
          relic_age_null_tbl$age_mode == age_mode &
            relic_age_null_tbl$age_stratum == age_stratum
        ]
      )[1],
      ci_anchor_share_diff_high = ci95(
        relic_age_null_tbl$anchor_share_diff[
          relic_age_null_tbl$age_mode == age_mode &
            relic_age_null_tbl$age_stratum == age_stratum
        ]
      )[2],
      p_one_sided_anchor_share_diff = one_sided_greater_perm_p(
        relic_age_null_tbl$anchor_share_diff[
          relic_age_null_tbl$age_mode == age_mode &
            relic_age_null_tbl$age_stratum == age_stratum
        ],
        anchor_share_diff
      ),
      p_two_sided_anchor_share_diff = two_sided_perm_p(
        relic_age_null_tbl$anchor_share_diff[
          relic_age_null_tbl$age_mode == age_mode &
            relic_age_null_tbl$age_stratum == age_stratum
        ],
        anchor_share_diff
      ),
      ci_log2_or_low = ci95(
        relic_age_null_tbl$anchor_log2_odds_ratio[
          relic_age_null_tbl$age_mode == age_mode &
            relic_age_null_tbl$age_stratum == age_stratum
        ]
      )[1],
      ci_log2_or_high = ci95(
        relic_age_null_tbl$anchor_log2_odds_ratio[
          relic_age_null_tbl$age_mode == age_mode &
            relic_age_null_tbl$age_stratum == age_stratum
        ]
      )[2],
      p_one_sided_log2_or = one_sided_greater_perm_p(
        relic_age_null_tbl$anchor_log2_odds_ratio[
          relic_age_null_tbl$age_mode == age_mode &
            relic_age_null_tbl$age_stratum == age_stratum
        ],
        anchor_log2_odds_ratio
      ),
      p_two_sided_log2_or = two_sided_perm_p(
        relic_age_null_tbl$anchor_log2_odds_ratio[
          relic_age_null_tbl$age_mode == age_mode &
            relic_age_null_tbl$age_stratum == age_stratum
        ],
        anchor_log2_odds_ratio
      ),
      support_label = dplyr::case_when(
        panel_support == "underpowered" ~ "Underpowered",
        is.finite(anchor_share_diff) &
          anchor_share_diff > 0 &
          is.finite(p_one_sided_anchor_share_diff) &
          p_one_sided_anchor_share_diff < 0.05 ~ "Directional persistence supported",
        TRUE ~ "Directional persistence not supported"
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(age_mode, age_stratum)
  
  list(
    relic_tbl = relic_tbl,
    relic_age_counts_tbl = relic_age_counts_tbl,
    relic_age_obs_tbl = relic_age_obs_tbl,
    relic_age_null_tbl = relic_age_null_tbl,
    relic_age_tests_tbl = relic_age_tests_tbl
  )
}

run_leave_one_region_out <- function(branch_tbl, anchor_def_tbl, formal_tests) {
  primary_regions <- extract_anchor_regions(anchor_def_tbl, PRIMARY_ANCHOR_ID)
  base_primary_relic <- formal_tests |>
    dplyr::filter(anchor_id == PRIMARY_ANCHOR_ID, panel == "relic")
  
  purrr::map_dfr(primary_regions, function(region_drop) {
    reduced_anchor <- setdiff(primary_regions, region_drop)
    
    df_relic <- branch_tbl |>
      dplyr::filter(class == "Relic") |>
      dplyr::mutate(
        inside_anchor = dominant_region %in% reduced_anchor
      )
    
    met <- compute_anchor_metrics(df_relic, focal_group = "PH908", control_group = "R1a")
    
    met |>
      dplyr::mutate(
        dropped_region = region_drop,
        anchor_id = PRIMARY_ANCHOR_ID,
        original_anchor_region_count = length(primary_regions),
        reduced_anchor_region_count = length(reduced_anchor),
        baseline_anchor_share_diff = base_primary_relic$anchor_share_diff[[1]],
        attenuation_vs_baseline = anchor_share_diff - baseline_anchor_share_diff,
        sign_positive = anchor_share_diff > 0
      )
  }) |>
    dplyr::arrange(dropped_region)
}

make_figuredata_main <- function(final_summary_tbl) {
  final_summary_tbl |>
    dplyr::filter(
      anchor_id == PRIMARY_ANCHOR_ID,
      panel %in% c("relic", "middle", "all", "founder")
    ) |>
    dplyr::mutate(
      panel = factor(panel, levels = c("relic", "middle", "all", "founder")),
      interpretation_bucket = dplyr::case_when(
        panel == "relic" ~ "Primary inference",
        panel %in% c("middle", "all") ~ "Secondary support",
        panel == "founder" ~ "Descriptive only",
        TRUE ~ "Other"
      )
    ) |>
    dplyr::arrange(panel)
}

make_figuredata_sens <- function(final_summary_tbl) {
  final_summary_tbl |>
    dplyr::filter(
      panel == "relic",
      anchor_id %in% c("cum50", "cum60", "cum70", "top4", "top5")
    ) |>
    dplyr::mutate(
      anchor_id = factor(anchor_id, levels = c("cum50", "cum60", "cum70", "top4", "top5"))
    ) |>
    dplyr::arrange(anchor_id)
}

make_figuredata_supp <- function(packet_comp_tbl, relic_age_tests_tbl, primary_anchor_id) {
  packet_tbl <- packet_comp_tbl |>
    dplyr::filter(
      anchor_id == primary_anchor_id,
      panel == "relic",
      inside_anchor %in% TRUE
    ) |>
    dplyr::mutate(
      figure_component = "packet_composition"
    )
  
  relic_age_tbl <- relic_age_tests_tbl |>
    dplyr::mutate(
      figure_component = "relic_age_extension",
      dominant_region = NA_character_,
      hg_group = NA_character_,
      branch_n = NA_integer_,
      branch_share = NA_real_,
      inside_anchor = NA
    ) |>
    dplyr::select(
      figure_component,
      age_mode,
      age_stratum,
      age_mode_pretty,
      age_stratum_pretty,
      focal_n,
      control_n,
      focal_anchor_n,
      control_anchor_n,
      anchor_share_focal,
      anchor_share_control,
      anchor_share_diff,
      anchor_log2_odds_ratio,
      p_one_sided_anchor_share_diff,
      p_one_sided_log2_or,
      support_label,
      dominant_region,
      hg_group,
      branch_n,
      branch_share,
      inside_anchor
    )
  
  packet_tbl |>
    dplyr::mutate(
      age_mode = NA_character_,
      age_stratum = NA_character_,
      age_mode_pretty = NA_character_,
      age_stratum_pretty = NA_character_,
      focal_n = NA_integer_,
      control_n = NA_integer_,
      focal_anchor_n = NA_integer_,
      control_anchor_n = NA_integer_,
      anchor_share_focal = NA_real_,
      anchor_share_control = NA_real_,
      anchor_share_diff = NA_real_,
      anchor_log2_odds_ratio = NA_real_,
      p_one_sided_anchor_share_diff = NA_real_,
      p_one_sided_log2_or = NA_real_,
      support_label = NA_character_
    ) |>
    dplyr::select(
      figure_component,
      age_mode,
      age_stratum,
      age_mode_pretty,
      age_stratum_pretty,
      focal_n,
      control_n,
      focal_anchor_n,
      control_anchor_n,
      anchor_share_focal,
      anchor_share_control,
      anchor_share_diff,
      anchor_log2_odds_ratio,
      p_one_sided_anchor_share_diff,
      p_one_sided_log2_or,
      support_label,
      dominant_region,
      hg_group,
      branch_n,
      branch_share,
      inside_anchor
    ) |>
    dplyr::bind_rows(relic_age_tbl)
}

make_supp_table <- function(age_diag, ph908_relic_rank, anchor_def_tbl, final_summary_tbl, relic_age_tests_tbl) {
  age_long <- age_diag |>
    dplyr::select(
      hg_group,
      n_branches,
      total_mass,
      n_relic,
      n_middle,
      n_founder,
      age_from_age_file,
      age_from_dataset,
      median_branch_size,
      max_branch_size
    ) |>
    tidyr::pivot_longer(
      cols = -hg_group,
      names_to = "metric",
      values_to = "value"
    ) |>
    tidyr::pivot_wider(
      names_from = hg_group,
      values_from = value
    )
  
  rank_tbl <- ph908_relic_rank |>
    dplyr::transmute(
      metric = paste0("ph908_relic_rank__", region_rank, "__", dominant_region),
      PH908 = relic_branches,
      R1a = region_share
    )
  
  anchor_tbl <- anchor_def_tbl |>
    dplyr::transmute(
      metric = paste0("anchor_def__", anchor_id),
      PH908 = anchor_region_count,
      R1a = anchor_branch_share_PH908_source
    )
  
  panel_tbl <- final_summary_tbl |>
    dplyr::filter(anchor_id %in% c(PRIMARY_ANCHOR_ID, "cum50", "cum70", "top4", "top5")) |>
    dplyr::transmute(
      metric = paste0("panel__", panel, "__", anchor_id),
      PH908 = anchor_share_focal,
      R1a = anchor_share_control
    )
  
  age_ext_tbl <- relic_age_tests_tbl |>
    dplyr::transmute(
      metric = paste0("relic_age__", age_mode, "__", age_stratum),
      PH908 = anchor_share_focal,
      R1a = anchor_share_control
    )
  
  dplyr::bind_rows(age_long, rank_tbl, anchor_tbl, panel_tbl, age_ext_tbl)
}

# Input preparation

cl_raw <- readRDS(FILE_PREPROCESSED_RDS)
age_tbl <- standardize_age_table(readr::read_csv(FILE_BRANCH_AGES, show_col_types = FALSE))

filter_diag <- make_filter_diag(
  if ("is_balkan_country" %in% names(cl_raw)) {
    cl_raw
  } else if ("country_code" %in% names(cl_raw)) {
    cl_raw |>
      dplyr::mutate(is_balkan_country = toupper(country_code) %in% BALKAN_ISO2)
  } else {
    cl_raw
  }
)

base_rows <- prepare_input_rows(cl_raw)
branch_tbl_main <- build_branch_table(base_rows, age_tbl)
age_diag <- make_age_diag(branch_tbl_main)

# Anchor definition

ph908_relic_rank <- branch_tbl_main |>
  dplyr::filter(
    hg_group == "PH908",
    class == "Relic",
    !is.na(dominant_region),
    dominant_region != ""
  ) |>
  dplyr::count(dominant_region, name = "relic_branches") |>
  dplyr::arrange(dplyr::desc(relic_branches), dominant_region) |>
  dplyr::mutate(
    region_rank = dplyr::row_number(),
    region_share = relic_branches / sum(relic_branches),
    cum_share = cumsum(region_share)
  ) |>
  dplyr::select(region_rank, dominant_region, relic_branches, region_share, cum_share)

if (nrow(ph908_relic_rank) == 0) {
  stop("No PH908 relic regions available for anchor definition.", call. = FALSE)
}

anchor_def_tbl <- build_anchor_definition_tbl(ph908_relic_rank)
assert_true(PRIMARY_ANCHOR_ID %in% anchor_def_tbl$anchor_id, "Primary anchor not found in anchor definition table.")

# Panel construction

panel_branch_tbl <- build_panel_branch_tbl(branch_tbl_main)

panel_counts <- panel_branch_tbl |>
  dplyr::count(panel, hg_group, class, name = "n_branches") |>
  dplyr::arrange(panel, hg_group, class)

observed_parts <- make_panelized_observed(panel_branch_tbl, anchor_def_tbl)
observed_tbl <- observed_parts$observed_tbl
packet_comp_tbl <- observed_parts$packet_comp_tbl
strata_diag <- make_strata_diag(panel_branch_tbl)

# Null-model permutations

null_tbl <- run_main_permutation_null(
  panel_branch_tbl = panel_branch_tbl,
  anchor_def_tbl = anchor_def_tbl,
  n_perm = N_PERM_MAIN
)

formal_tests <- build_formal_tests(observed_tbl, null_tbl)
final_summary_tbl <- make_final_summary(formal_tests)
sensitivity_tbl <- final_summary_tbl |>
  dplyr::filter(panel %in% c("relic", "middle", "all"), anchor_id != PRIMARY_ANCHOR_ID) |>
  dplyr::arrange(panel, anchor_id)

# Relic age extension analysis

primary_anchor_regions <- extract_anchor_regions(anchor_def_tbl, PRIMARY_ANCHOR_ID)
relic_age_parts <- run_relic_age_extension(
  branch_tbl = branch_tbl_main,
  primary_anchor_regions = primary_anchor_regions,
  n_perm = N_PERM_RELIC_AGE
)

relic_tbl <- relic_age_parts$relic_tbl
relic_age_counts_tbl <- relic_age_parts$relic_age_counts_tbl
relic_age_obs_tbl <- relic_age_parts$relic_age_obs_tbl
relic_age_null_tbl <- relic_age_parts$relic_age_null_tbl
relic_age_tests_tbl <- relic_age_parts$relic_age_tests_tbl

# Robustness analysis

table_p3_rx <- run_leave_one_region_out(branch_tbl_main, anchor_def_tbl, formal_tests)

# Figure data and tables

figuredata_p3_main <- make_figuredata_main(final_summary_tbl)
figuredata_p3_sens <- make_figuredata_sens(final_summary_tbl)
figuredata_p3_supp <- make_figuredata_supp(packet_comp_tbl, relic_age_tests_tbl, PRIMARY_ANCHOR_ID)
table_p3_supp <- make_supp_table(age_diag, ph908_relic_rank, anchor_def_tbl, final_summary_tbl, relic_age_tests_tbl)

# Consistency checks

age_ph908 <- age_diag |>
  dplyr::filter(hg_group == "PH908")
age_r1a <- age_diag |>
  dplyr::filter(hg_group == "R1a")

assert_equal(age_ph908$n_branches[[1]], 79L, "PH908 branch count drift in Pillar 3.")
assert_equal(age_ph908$total_mass[[1]], 623L, "PH908 total branch mass drift in Pillar 3.")
assert_equal(age_ph908$n_relic[[1]], 50L, "PH908 relic count drift in Pillar 3.")
assert_equal(age_ph908$n_middle[[1]], 20L, "PH908 middle count drift in Pillar 3.")
assert_equal(age_ph908$n_founder[[1]], 9L, "PH908 founder count drift in Pillar 3.")

assert_equal(age_r1a$n_branches[[1]], 81L, "R1a branch count drift in Pillar 3.")
assert_equal(age_r1a$total_mass[[1]], 729L, "R1a total branch mass drift in Pillar 3.")
assert_equal(age_r1a$n_relic[[1]], 38L, "R1a relic count drift in Pillar 3.")
assert_equal(age_r1a$n_middle[[1]], 29L, "R1a middle count drift in Pillar 3.")
assert_equal(age_r1a$n_founder[[1]], 14L, "R1a founder count drift in Pillar 3.")

primary_anchor <- anchor_def_tbl |>
  dplyr::filter(anchor_id == PRIMARY_ANCHOR_ID)

assert_equal(primary_anchor$anchor_region_count[[1]], 5L, "Primary anchor region count drift in Pillar 3.")

expected_regions <- c("Босанска Крајина", "Брда", "Стари Влах", "Косово", "Херцеговина")
actual_regions <- extract_anchor_regions(anchor_def_tbl, PRIMARY_ANCHOR_ID)
assert_true(
  identical(sort(actual_regions), sort(expected_regions)),
  paste0(
    "Primary anchor region set drift in Pillar 3.\n  actual: ",
    paste(actual_regions, collapse = ", "),
    "\n  expected: ",
    paste(expected_regions, collapse = ", ")
  )
)

row_relic_primary <- final_summary_tbl |>
  dplyr::filter(anchor_id == PRIMARY_ANCHOR_ID, panel == "relic")
row_middle_primary <- final_summary_tbl |>
  dplyr::filter(anchor_id == PRIMARY_ANCHOR_ID, panel == "middle")
row_all_primary <- final_summary_tbl |>
  dplyr::filter(anchor_id == PRIMARY_ANCHOR_ID, panel == "all")
row_founder_primary <- final_summary_tbl |>
  dplyr::filter(anchor_id == PRIMARY_ANCHOR_ID, panel == "founder")

assert_equal(row_relic_primary$focal_anchor_n[[1]], 30L, "Primary relic PH908 anchor count drift.")
assert_equal(row_relic_primary$focal_n[[1]], 50L, "Primary relic PH908 denominator drift.")
assert_equal(row_relic_primary$control_anchor_n[[1]], 4L, "Primary relic R1a anchor count drift.")
assert_equal(row_relic_primary$control_n[[1]], 38L, "Primary relic R1a denominator drift.")
assert_close(row_relic_primary$anchor_share_focal[[1]], 0.600, tol = 0.0015, "Primary relic PH908 anchor share drift.")
assert_close(row_relic_primary$anchor_share_control[[1]], 0.105, tol = 0.0015, "Primary relic R1a anchor share drift.")
assert_close(row_relic_primary$anchor_share_diff[[1]], 0.495, tol = 0.0015, "Primary relic anchor share diff drift.")
assert_close(row_relic_primary$anchor_log2_odds_ratio[[1]], 3.512, tol = 0.02, "Primary relic log2 OR drift.")

assert_close(row_middle_primary$anchor_share_focal[[1]], 0.550, tol = 0.0015, "Primary middle PH908 anchor share drift.")
assert_close(row_middle_primary$anchor_share_control[[1]], 0.345, tol = 0.0015, "Primary middle R1a anchor share drift.")
assert_close(row_middle_primary$anchor_share_diff[[1]], 0.205, tol = 0.0015, "Primary middle anchor share diff drift.")
assert_close(row_middle_primary$anchor_log2_odds_ratio[[1]], 1.169, tol = 0.02, "Primary middle log2 OR drift.")

assert_equal(row_all_primary$focal_anchor_n[[1]], 48L, "Primary all-branch PH908 anchor count drift.")
assert_equal(row_all_primary$focal_n[[1]], 79L, "Primary all-branch PH908 denominator drift.")
assert_equal(row_all_primary$control_anchor_n[[1]], 24L, "Primary all-branch R1a anchor count drift.")
assert_equal(row_all_primary$control_n[[1]], 81L, "Primary all-branch R1a denominator drift.")
assert_close(row_all_primary$anchor_share_focal[[1]], 0.608, tol = 0.0015, "Primary all-branch PH908 anchor share drift.")
assert_close(row_all_primary$anchor_share_control[[1]], 0.296, tol = 0.0015, "Primary all-branch R1a anchor share drift.")
assert_close(row_all_primary$anchor_share_diff[[1]], 0.311, tol = 0.0015, "Primary all-branch anchor share diff drift.")
assert_close(row_all_primary$anchor_log2_odds_ratio[[1]], 1.853, tol = 0.02, "Primary all-branch log2 OR drift.")

assert_close(row_founder_primary$anchor_share_focal[[1]], 0.778, tol = 0.0015, "Primary founder PH908 anchor share drift.")
assert_close(row_founder_primary$anchor_share_control[[1]], 0.714, tol = 0.0015, "Primary founder R1a anchor share drift.")

row_cum50 <- final_summary_tbl |>
  dplyr::filter(anchor_id == "cum50", panel == "relic")
row_cum70 <- final_summary_tbl |>
  dplyr::filter(anchor_id == "cum70", panel == "relic")
row_top4 <- final_summary_tbl |>
  dplyr::filter(anchor_id == "top4", panel == "relic")
row_top5 <- final_summary_tbl |>
  dplyr::filter(anchor_id == "top5", panel == "relic")

assert_close(row_cum50$anchor_share_focal[[1]], 0.500, tol = 0.0015, "cum50 relic PH908 share drift.")
assert_close(row_cum50$anchor_share_control[[1]], 0.053, tol = 0.0015, "cum50 relic R1a share drift.")
assert_close(row_cum50$anchor_share_diff[[1]], 0.447, tol = 0.0015, "cum50 relic diff drift.")
assert_close(row_cum50$anchor_log2_odds_ratio[[1]], 3.868, tol = 0.02, "cum50 relic log2 OR drift.")

assert_close(row_cum70$anchor_share_focal[[1]], 0.720, tol = 0.0015, "cum70 relic PH908 share drift.")
assert_close(row_cum70$anchor_share_control[[1]], 0.263, tol = 0.0015, "cum70 relic R1a share drift.")
assert_close(row_cum70$anchor_share_diff[[1]], 0.457, tol = 0.0015, "cum70 relic diff drift.")
assert_close(row_cum70$anchor_log2_odds_ratio[[1]], 2.772, tol = 0.02, "cum70 relic log2 OR drift.")

assert_close(row_top4$anchor_share_diff[[1]], 0.447, tol = 0.0015, "top4 relic diff drift.")
assert_close(row_top5$anchor_share_diff[[1]], 0.495, tol = 0.0015, "top5 relic diff drift.")

assert_true(
  all(table_p3_rx$sign_positive, na.rm = TRUE),
  "At least one leave-one-region-out run flipped the primary relic anchor contrast non-positive."
)

if (N_PERM_MAIN >= 5000L) {
  assert_close(row_relic_primary$p_one_sided_anchor_share_diff[[1]], 0.0002, tol = 0.001, "Primary relic share p-value drift.")
  assert_close(row_relic_primary$p_one_sided_log2_or[[1]], 0.0002, tol = 0.001, "Primary relic log2 OR p-value drift.")
  assert_close(row_all_primary$p_one_sided_anchor_share_diff[[1]], 0.0012, tol = 0.002, "Primary all-branch share p-value drift.")
  assert_close(row_all_primary$p_one_sided_log2_or[[1]], 0.0010, tol = 0.002, "Primary all-branch log2 OR p-value drift.")
}

# Facts output
writeLines(character(), FILE_P3_FACTS)
writeLines(
  c(
    "script: 04_pillar3_anchor.R",
    paste0("timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("seed: ", SEED_MAIN),
    paste0("config_file: ", path_for_facts(CONFIG_FILE)),
    "project_root: .",
    paste0("git_commit: ", git_commit),
    "inputs:",
    paste0("  - ", path_for_facts(FILE_PREPROCESSED_RDS)),
    paste0("  - ", path_for_facts(FILE_BRANCH_AGES)),
    ""
  ),
  FILE_P3_FACTS,
  useBytes = TRUE
)

results_sentence <- paste0(
  "Under the pre-specified cum60 PH908-defined relic anchor packet, PH908 placed ",
  row_relic_primary$focal_anchor_n[[1]], " of ", row_relic_primary$focal_n[[1]],
  " relic branches inside the packet (", sprintf("%.3f", row_relic_primary$anchor_share_focal[[1]]),
  "), whereas R1a placed ", row_relic_primary$control_anchor_n[[1]], " of ", row_relic_primary$control_n[[1]],
  " (", sprintf("%.3f", row_relic_primary$anchor_share_control[[1]]),
  "), yielding an anchor-share difference of ",
  sprintf("%.3f", row_relic_primary$anchor_share_diff[[1]]),
  " and a log2 anchor odds ratio of ",
  sprintf("%.3f", row_relic_primary$anchor_log2_odds_ratio[[1]]), "."
)

caption_sentence <- paste0(
  "Primary anchor = cum60, defined from the PH908 relic regional rank and consisting of ",
  primary_anchor$anchor_region_count[[1]],
  " regions: ",
  paste(actual_regions, collapse = ", "),
  ". The relic panel is the primary inferential target."
)

sensitivity_sentence <- paste0(
  "Relic-panel sensitivity remained positive across cum50, cum60, cum70, top4, and top5, with leave-one-region-out robustness preserving a positive primary relic contrast in all ",
  nrow(table_p3_rx),
  " reduced-anchor runs."
)

append_facts_lines(
  FILE_P3_FACTS,
  c(
    paste0("ph908_gating_mode: ", PH908_GATING_MODE),
    paste0("n_perm_main: ", N_PERM_MAIN),
    paste0("n_perm_relic_age: ", N_PERM_RELIC_AGE),
    paste0("age_bin_n: ", AGE_BIN_N),
    paste0("size_bin_n: ", SIZE_BIN_N),
    paste0("primary_anchor_id: ", PRIMARY_ANCHOR_ID),
    paste0("anchor_top_k_values: ", paste(ANCHOR_TOP_K_VALUES, collapse = ",")),
    paste0("anchor_cum_values: ", paste(sprintf("%.2f", ANCHOR_CUM_VALUES), collapse = ",")),
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
    "primary_anchor:",
    paste0("  anchor_id: ", PRIMARY_ANCHOR_ID),
    paste0("  anchor_region_count: ", primary_anchor$anchor_region_count[[1]]),
    paste0("  anchor_regions: ", paste(actual_regions, collapse = " | ")),
    paste0("  anchor_branch_share_PH908_source: ", sprintf("%.3f", primary_anchor$anchor_branch_share_PH908_source[[1]])),
    "",
    "primary_panel_summary:",
    paste0("  relic_PH908_count_share: ", row_relic_primary$focal_anchor_n[[1]], "/", row_relic_primary$focal_n[[1]], " / ", sprintf("%.3f", row_relic_primary$anchor_share_focal[[1]])),
    paste0("  relic_R1a_count_share: ", row_relic_primary$control_anchor_n[[1]], "/", row_relic_primary$control_n[[1]], " / ", sprintf("%.3f", row_relic_primary$anchor_share_control[[1]])),
    paste0("  relic_diff_log2or: ", sprintf("%.3f", row_relic_primary$anchor_share_diff[[1]]), " / ", sprintf("%.3f", row_relic_primary$anchor_log2_odds_ratio[[1]])),
    paste0("  middle_diff_log2or: ", sprintf("%.3f", row_middle_primary$anchor_share_diff[[1]]), " / ", sprintf("%.3f", row_middle_primary$anchor_log2_odds_ratio[[1]])),
    paste0("  all_diff_log2or: ", sprintf("%.3f", row_all_primary$anchor_share_diff[[1]]), " / ", sprintf("%.3f", row_all_primary$anchor_log2_odds_ratio[[1]])),
    paste0("  founder_diff: ", sprintf("%.3f", row_founder_primary$anchor_share_diff[[1]])),
    "",
    "relic_sensitivity_summary:",
    paste0("  cum50_diff_log2or: ", sprintf("%.3f", row_cum50$anchor_share_diff[[1]]), " / ", sprintf("%.3f", row_cum50$anchor_log2_odds_ratio[[1]])),
    paste0("  cum60_diff_log2or: ", sprintf("%.3f", row_relic_primary$anchor_share_diff[[1]]), " / ", sprintf("%.3f", row_relic_primary$anchor_log2_odds_ratio[[1]])),
    paste0("  cum70_diff_log2or: ", sprintf("%.3f", row_cum70$anchor_share_diff[[1]]), " / ", sprintf("%.3f", row_cum70$anchor_log2_odds_ratio[[1]])),
    paste0("  top4_diff_log2or: ", sprintf("%.3f", row_top4$anchor_share_diff[[1]]), " / ", sprintf("%.3f", row_top4$anchor_log2_odds_ratio[[1]])),
    paste0("  top5_diff_log2or: ", sprintf("%.3f", row_top5$anchor_share_diff[[1]]), " / ", sprintf("%.3f", row_top5$anchor_log2_odds_ratio[[1]])),
    "",
    "leave_one_region_out:",
    paste0("  n_runs: ", nrow(table_p3_rx)),
    paste0("  sign_flip_anywhere: ", ifelse(all(table_p3_rx$sign_positive), "no", "yes")),
    paste0("  min_reduced_anchor_diff: ", sprintf("%.3f", min(table_p3_rx$anchor_share_diff, na.rm = TRUE))),
    paste0("  max_reduced_anchor_diff: ", sprintf("%.3f", max(table_p3_rx$anchor_share_diff, na.rm = TRUE))),
    "",
    "relic_age_extension:",
    paste0(
      "  all_relics_diff_log2or: ",
      sprintf("%.3f", relic_age_tests_tbl$anchor_share_diff[relic_age_tests_tbl$age_mode == "all_relics" & relic_age_tests_tbl$age_stratum == "all_relics"]),
      " / ",
      sprintf("%.3f", relic_age_tests_tbl$anchor_log2_odds_ratio[relic_age_tests_tbl$age_mode == "all_relics" & relic_age_tests_tbl$age_stratum == "all_relics"])
    ),
    paste0(
      "  old_vs_rest_old_diff_log2or: ",
      sprintf("%.3f", relic_age_tests_tbl$anchor_share_diff[relic_age_tests_tbl$age_mode == "old_vs_rest_median" & relic_age_tests_tbl$age_stratum == "old"]),
      " / ",
      sprintf("%.3f", relic_age_tests_tbl$anchor_log2_odds_ratio[relic_age_tests_tbl$age_mode == "old_vs_rest_median" & relic_age_tests_tbl$age_stratum == "old"])
    ),
    paste0(
      "  old_vs_rest_rest_diff_log2or: ",
      sprintf("%.3f", relic_age_tests_tbl$anchor_share_diff[relic_age_tests_tbl$age_mode == "old_vs_rest_median" & relic_age_tests_tbl$age_stratum == "rest"]),
      " / ",
      sprintf("%.3f", relic_age_tests_tbl$anchor_log2_odds_ratio[relic_age_tests_tbl$age_mode == "old_vs_rest_median" & relic_age_tests_tbl$age_stratum == "rest"])
    ),
    paste0(
      "  tertiles_oldest_diff_log2or: ",
      sprintf("%.3f", relic_age_tests_tbl$anchor_share_diff[relic_age_tests_tbl$age_mode == "within_lineage_tertiles" & relic_age_tests_tbl$age_stratum == "oldest"]),
      " / ",
      sprintf("%.3f", relic_age_tests_tbl$anchor_log2_odds_ratio[relic_age_tests_tbl$age_mode == "within_lineage_tertiles" & relic_age_tests_tbl$age_stratum == "oldest"])
    ),
    paste0(
      "  tertiles_middle_diff_log2or: ",
      sprintf("%.3f", relic_age_tests_tbl$anchor_share_diff[relic_age_tests_tbl$age_mode == "within_lineage_tertiles" & relic_age_tests_tbl$age_stratum == "middle"]),
      " / ",
      sprintf("%.3f", relic_age_tests_tbl$anchor_log2_odds_ratio[relic_age_tests_tbl$age_mode == "within_lineage_tertiles" & relic_age_tests_tbl$age_stratum == "middle"])
    ),
    paste0(
      "  tertiles_youngest_diff_log2or: ",
      sprintf("%.3f", relic_age_tests_tbl$anchor_share_diff[relic_age_tests_tbl$age_mode == "within_lineage_tertiles" & relic_age_tests_tbl$age_stratum == "youngest"]),
      " / ",
      sprintf("%.3f", relic_age_tests_tbl$anchor_log2_odds_ratio[relic_age_tests_tbl$age_mode == "within_lineage_tertiles" & relic_age_tests_tbl$age_stratum == "youngest"])
    ),
    paste0("results_sentence: ", results_sentence),
    paste0("caption_sentence: ", caption_sentence),
    paste0("sensitivity_sentence: ", sensitivity_sentence),
    "note: p-value consistency checks are conditional on the configured permutation depth.",
    "output_files:",
    paste0("  - ", path_for_facts(FILE_P3_RESULTS)),
    paste0("  - ", path_for_facts(FILE_P3_FIG_MAIN)),
    paste0("  - ", path_for_facts(FILE_P3_FIG_SENS)),
    paste0("  - ", path_for_facts(FILE_P3_FIG_SUPP)),
    paste0("  - ", path_for_facts(FILE_P3_SUPP_TABLE)),
    paste0("  - ", path_for_facts(FILE_P3_RX_TABLE)),
    paste0("  - ", path_for_facts(FILE_P3_FACTS))
  )
)

# Exports

readr::write_csv(figuredata_p3_main, FILE_P3_FIG_MAIN)
readr::write_csv(figuredata_p3_sens, FILE_P3_FIG_SENS)
readr::write_csv(figuredata_p3_supp, FILE_P3_FIG_SUPP)
readr::write_csv(table_p3_supp, FILE_P3_SUPP_TABLE)
readr::write_csv(table_p3_rx, FILE_P3_RX_TABLE)

results_p3 <- list(
  metadata = list(
    script = "04_pillar3_anchor.R",
    seed = SEED_MAIN,
    config_file = path_for_facts(CONFIG_FILE),
    inputs = c(
      path_for_facts(FILE_PREPROCESSED_RDS),
      path_for_facts(FILE_BRANCH_AGES)
    ),
    ph908_gating_mode = PH908_GATING_MODE,
    age_bin_n = AGE_BIN_N,
    size_bin_n = SIZE_BIN_N,
    n_perm_main = N_PERM_MAIN,
    n_perm_relic_age = N_PERM_RELIC_AGE,
    primary_anchor_id = PRIMARY_ANCHOR_ID,
    anchor_top_k_values = ANCHOR_TOP_K_VALUES,
    anchor_cum_values = ANCHOR_CUM_VALUES,
    git_commit = git_commit
  ),
  filter_diag = filter_diag,
  age_diag = age_diag,
  branch_table = branch_tbl_main,
  ph908_relic_rank = ph908_relic_rank,
  anchor_definitions = anchor_def_tbl,
  panel_counts = panel_counts,
  strata_diag = strata_diag,
  observed_metrics = observed_tbl,
  packet_composition = packet_comp_tbl,
  null_main = null_tbl,
  formal_tests = formal_tests,
  final_summary = final_summary_tbl,
  sensitivity_summary = sensitivity_tbl,
  relic_table = relic_tbl,
  relic_age_counts = relic_age_counts_tbl,
  relic_age_observed = relic_age_obs_tbl,
  relic_age_null = relic_age_null_tbl,
  relic_age_tests = relic_age_tests_tbl,
  leave_one_region_out = table_p3_rx,
  figuredata_main = figuredata_p3_main,
  figuredata_sens = figuredata_p3_sens,
  figuredata_supp = figuredata_p3_supp,
  table_supp = table_p3_supp
)

saveRDS(results_p3, FILE_P3_RESULTS)

cat("============================================================\n")
cat("PILLAR 3 COMPLETE\n")
cat("============================================================\n")
cat("Results object: ", FILE_P3_RESULTS, "\n", sep = "")
cat("Facts file: ", FILE_P3_FACTS, "\n", sep = "")
cat("Figuredata main: ", FILE_P3_FIG_MAIN, "\n", sep = "")
cat("Figuredata sens: ", FILE_P3_FIG_SENS, "\n", sep = "")
cat("Figuredata supp: ", FILE_P3_FIG_SUPP, "\n", sep = "")
cat("Supplement table: ", FILE_P3_SUPP_TABLE, "\n", sep = "")
cat("Robustness table: ", FILE_P3_RX_TABLE, "\n", sep = "")
cat("============================================================\n")