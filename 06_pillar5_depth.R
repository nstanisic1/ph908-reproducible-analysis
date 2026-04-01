# 06_pillar5_depth_coherence.R
# Depth-structured regional coherence under matched branch-link rules.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

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

SCRIPT_NAME <- "06_pillar5_depth_coherence.R"
DATA_FILE <- FILE_PREPROCESSED_RDS
AGE_FILE <- FILE_BRANCH_AGES

FILE_P5_RESULTS <- file.path(DIR_RESULTS, "results_p5.rds")
FILE_P5_FIG_MAIN <- file.path(DIR_FIGUREDATA, "figuredata_p5_main.csv")
FILE_P5_FIG_INSET <- file.path(DIR_FIGUREDATA, "figuredata_p5_inset.csv")
FILE_P5_FIG_NULL <- file.path(DIR_FIGUREDATA, "figuredata_p5_null.csv")
FILE_P5_SUPP_TABLE <- file.path(DIR_TABLES, "table_p5_supp.csv")
FILE_P5_FACTS <- file.path(DIR_FACTS, "facts_p5.txt")
FILE_P5_WARNINGS <- file.path(DIR_FACTS, "warnings_p5.txt")

FILE_P5_FILTER_DIAG <- file.path(DIR_TABLES, "table_p5_filter_diagnostics.csv")
FILE_P5_AGE_DIAG <- file.path(DIR_TABLES, "table_p5_age_join_diagnostics.csv")
FILE_P5_BRANCH <- file.path(DIR_TABLES, "table_p5_branch_table.csv")
FILE_P5_EDGES <- file.path(DIR_TABLES, "table_p5_edges_observed.csv")
FILE_P5_PROFILE <- file.path(DIR_TABLES, "table_p5_depth_profile_observed.csv")
FILE_P5_PROFILE_SUM <- file.path(DIR_TABLES, "table_p5_profile_summary.csv")
FILE_P5_TESTS <- file.path(DIR_TABLES, "table_p5_formal_tests.csv")
FILE_P5_NULL <- file.path(DIR_TABLES, "table_p5_permutation_null.csv")
FILE_P5_FINAL <- file.path(DIR_TABLES, "table_p5_final_summary.csv")
FILE_P5_NET_DIAG <- file.path(DIR_TABLES, "table_p5_network_diagnostics.csv")
FILE_P5_STRATA <- file.path(DIR_TABLES, "table_p5_permutation_strata_diagnostics.csv")

RNG_SEED <- SEED_MAIN
set.seed(RNG_SEED)

PH908_GATING_MODE <- "primary"
MAX_AGE_GAP_WINDOW <- 1200
LAMBDA_AGE_COST <- 1.0
N_PERM <- 5000L
RELIC_MAX_N <- 3L
FOUNDER_MIN_N <- 15L
AGE_BIN_N <- 4L
SIZE_BIN_N <- 4L
DEPTH_BIN_N <- 4L
MIN_SHARED_DEPTH_TIERS <- 2L
MIN_EDGES_PER_TIER <- 5L

stop_msg <- function(...) stop(paste0(...), call. = FALSE)

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

haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371.0
  to_rad <- function(x) x * pi / 180
  lat1 <- to_rad(lat1)
  lon1 <- to_rad(lon1)
  lat2 <- to_rad(lat2)
  lon2 <- to_rad(lon2)
  dlat <- lat2 - lat1
  dlon <- lon2 - lon1
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

permute_within_strata <- function(x, strata) {
  out <- x
  u <- unique(strata)
  for (s in u) {
    idx <- which(strata == s)
    if (length(idx) > 1) {
      out[idx] <- sample(out[idx], length(idx), replace = FALSE)
    }
  }
  out
}

make_quantile_bins <- function(x, n_bins = 4L, prefix = "bin") {
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

one_sided_p_greater <- function(null_vals, obs) {
  null_vals <- null_vals[is.finite(null_vals)]
  if (!is.finite(obs) || length(null_vals) == 0) return(NA_real_)
  (1 + sum(null_vals >= obs)) / (length(null_vals) + 1)
}

two_sided_p <- function(null_vals, obs) {
  null_vals <- null_vals[is.finite(null_vals)]
  if (!is.finite(obs) || length(null_vals) == 0) return(NA_real_)
  (1 + sum(abs(null_vals) >= abs(obs))) / (length(null_vals) + 1)
}

ci95 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 20) return(c(NA_real_, NA_real_))
  as.numeric(stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
}

safe_lm_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2 || length(unique(x)) < 2) return(NA_real_)
  coef(stats::lm(y ~ x))[2] %>% as.numeric()
}

write_warning_log <- function(warning_vec, path) {
  if (length(warning_vec) == 0) {
    writeLines("No warnings generated.", path)
  } else {
    writeLines(unique(warning_vec), path)
  }
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

warning_log <- character()

withCallingHandlers({
  
  if (!file.exists(DATA_FILE)) stop_msg("Missing data file: ", DATA_FILE)
  if (!file.exists(AGE_FILE)) stop_msg("Missing age file: ", AGE_FILE)
  
  cl <- readRDS(DATA_FILE)
  
  if (!"is_balkan_country" %in% names(cl)) {
    if (!"country_code" %in% names(cl)) {
      stop_msg("Missing both is_balkan_country and country_code in dataset")
    }
    cl <- cl %>%
      dplyr::mutate(is_balkan_country = toupper(country_code) %in% BALKAN_ISO2)
    message("Reconstructed is_balkan_country from country_code for Pillar 5 filtering.")
  }
  
  required_cols <- c(
    "matched", "exclude_geo_primary", "is_balkan_country",
    "terminal_snp", "Region", "lat", "long",
    "major_hg", "tmrca", "is_ph908_primary"
  )
  
  missing_required <- setdiff(required_cols, names(cl))
  if (length(missing_required) > 0) {
    stop_msg("Missing required columns in dataset: ", paste(missing_required, collapse = ", "))
  }
  
  age_tbl_raw <- readr::read_csv(AGE_FILE, show_col_types = FALSE)
  if (!all(c("terminal_snp", "tmrca_ybp") %in% names(age_tbl_raw))) {
    stop_msg("branch_ages.csv must contain terminal_snp and tmrca_ybp")
  }
  
  age_tbl <- age_tbl_raw %>%
    dplyr::mutate(
      terminal_snp_norm = normalize_snp(terminal_snp),
      tmrca_ybp = as.numeric(tmrca_ybp)
    ) %>%
    dplyr::filter(is.finite(tmrca_ybp), tmrca_ybp > 0) %>%
    dplyr::group_by(terminal_snp_norm) %>%
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
      nrow(cl),
      sum(cl$matched %in% TRUE, na.rm = TRUE),
      sum(cl$exclude_geo_primary %in% FALSE, na.rm = TRUE),
      sum(cl$is_balkan_country %in% TRUE, na.rm = TRUE),
      sum(!is.na(cl$terminal_snp) & cl$terminal_snp != ""),
      sum(!is.na(cl$Region) & cl$Region != ""),
      sum(is.finite(cl$lat) & is.finite(cl$long))
    )
  )
  readr::write_csv(filter_diag, FILE_P5_FILTER_DIAG)
  
  base <- cl %>%
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
    ) %>%
    dplyr::mutate(
      terminal_snp_norm = normalize_snp(terminal_snp),
      hg_group = dplyr::case_when(
        is_ph908_primary %in% TRUE ~ "PH908",
        major_hg == "R1a" ~ "R1a",
        TRUE ~ "Other"
      )
    ) %>%
    dplyr::filter(hg_group %in% c("PH908", "R1a"))
  
  if (nrow(base) == 0) stop_msg("No usable rows after filtering.")
  
  branch_tbl <- base %>%
    dplyr::group_by(hg_group, terminal_snp_norm) %>%
    dplyr::summarise(
      branch_size = dplyr::n(),
      lat_c = mean(lat, na.rm = TRUE),
      lon_c = mean(long, na.rm = TRUE),
      dominant_region = dominant_region_fun(Region),
      dataset_tmrca = safe_mean(tmrca),
      .groups = "drop"
    ) %>%
    dplyr::left_join(age_tbl, by = "terminal_snp_norm") %>%
    dplyr::mutate(
      tmrca_final = dplyr::coalesce(tmrca_ybp, dataset_tmrca),
      age_source = dplyr::case_when(
        is.finite(tmrca_ybp) ~ "branch_ages_csv",
        is.finite(dataset_tmrca) ~ "dataset_tmrca",
        TRUE ~ "missing"
      ),
      class = classify_branch(branch_size)
    ) %>%
    dplyr::filter(is.finite(tmrca_final))
  
  if (nrow(branch_tbl) == 0) stop_msg("No branches with usable TMRCA after merging age information.")
  
  readr::write_csv(branch_tbl, FILE_P5_BRANCH)
  
  age_diag <- branch_tbl %>%
    dplyr::group_by(hg_group) %>%
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
  readr::write_csv(age_diag, FILE_P5_AGE_DIAG)
  
  choose_parent_index <- function(current_row, older_df, rule, max_age_gap_window = 1200, lambda_age_cost = 1.0) {
    if (nrow(older_df) == 0) return(NA_integer_)
    dist_vec <- haversine_km(
      lat1 = current_row$lat_c[[1]],
      lon1 = current_row$lon_c[[1]],
      lat2 = older_df$lat_c,
      lon2 = older_df$lon_c
    )
    age_gap_vec <- older_df$tmrca_final - current_row$tmrca_final[[1]]
    valid <- is.finite(dist_vec) & is.finite(age_gap_vec) & age_gap_vec >= 0
    if (!any(valid)) return(NA_integer_)
    dist_vec <- dist_vec[valid]
    age_gap_vec <- age_gap_vec[valid]
    if (rule == "nearest_older") {
      return(which(valid)[which.min(dist_vec)])
    }
    if (rule == "nearest_older_window") {
      in_window <- age_gap_vec <= max_age_gap_window
      if (any(in_window)) {
        return(which(valid)[in_window][which.min(dist_vec[in_window])])
      } else {
        return(which(valid)[which.min(dist_vec)])
      }
    }
    if (rule == "distance_age_cost") {
      dist_scale <- stats::median(dist_vec, na.rm = TRUE)
      age_scale <- stats::median(age_gap_vec, na.rm = TRUE)
      if (!is.finite(dist_scale) || dist_scale <= 0) dist_scale <- 1
      if (!is.finite(age_scale) || age_scale <= 0) age_scale <- 1
      score <- (dist_vec / dist_scale) + lambda_age_cost * (age_gap_vec / age_scale)
      return(which(valid)[which.min(score)])
    }
    stop_msg("Unknown rule: ", rule)
  }
  
  build_network_by_rule <- function(df_group, rule, max_age_gap_window = 1200, lambda_age_cost = 1.0) {
    df_group <- df_group %>%
      dplyr::filter(
        is.finite(tmrca_final),
        is.finite(lat_c),
        is.finite(lon_c),
        !is.na(dominant_region),
        dominant_region != ""
      ) %>%
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
        parent_lon = numeric(),
        parent_lat = numeric(),
        child_lon = numeric(),
        child_lat = numeric(),
        edge_km = numeric(),
        parent_region = character(),
        child_region = character(),
        same_region = logical(),
        parent_class = character(),
        child_class = character(),
        rule = character(),
        hg_group = character()
      ))
    }
    
    edges <- vector("list", nrow(df_group) - 1)
    
    for (i in seq.int(2, nrow(df_group))) {
      current <- df_group[i, , drop = FALSE]
      older <- df_group[seq_len(i - 1), , drop = FALSE]
      
      j <- choose_parent_index(
        current_row = current,
        older_df = older,
        rule = rule,
        max_age_gap_window = max_age_gap_window,
        lambda_age_cost = lambda_age_cost
      )
      
      if (!is.na(j)) {
        edge_km <- haversine_km(
          lat1 = current$lat_c[[1]],
          lon1 = current$lon_c[[1]],
          lat2 = older$lat_c[[j]],
          lon2 = older$lon_c[[j]]
        )
        
        edges[[i - 1]] <- tibble::tibble(
          parent_snp = older$terminal_snp_norm[[j]],
          child_snp = current$terminal_snp_norm[[1]],
          parent_tmrca = older$tmrca_final[[j]],
          child_tmrca = current$tmrca_final[[1]],
          age_gap = older$tmrca_final[[j]] - current$tmrca_final[[1]],
          parent_lon = older$lon_c[[j]],
          parent_lat = older$lat_c[[j]],
          child_lon = current$lon_c[[1]],
          child_lat = current$lat_c[[1]],
          edge_km = edge_km,
          parent_region = older$dominant_region[[j]],
          child_region = current$dominant_region[[1]],
          same_region = older$dominant_region[[j]] == current$dominant_region[[1]],
          parent_class = older$class[[j]],
          child_class = current$class[[1]],
          rule = rule,
          hg_group = current$hg_group[[1]]
        )
      }
    }
    
    dplyr::bind_rows(edges)
  }
  
  build_all_rules_observed <- function(branch_df) {
    rules <- c("nearest_older", "nearest_older_window", "distance_age_cost")
    out <- vector("list", length(rules) * 2)
    k <- 1
    for (rr in rules) {
      for (gg in c("PH908", "R1a")) {
        out[[k]] <- build_network_by_rule(
          df_group = branch_df %>% dplyr::filter(hg_group == gg),
          rule = rr,
          max_age_gap_window = MAX_AGE_GAP_WINDOW,
          lambda_age_cost = LAMBDA_AGE_COST
        )
        k <- k + 1
      }
    }
    dplyr::bind_rows(out)
  }
  
  edges_obs <- build_all_rules_observed(branch_tbl)
  if (nrow(edges_obs) == 0) stop_msg("No observed edges were generated.")
  
  parent_age_breaks <- unique(stats::quantile(
    edges_obs$parent_tmrca[is.finite(edges_obs$parent_tmrca)],
    probs = seq(0, 1, length.out = DEPTH_BIN_N + 1),
    na.rm = TRUE
  ))
  
  actual_depth_bins <- length(parent_age_breaks) - 1
  if (actual_depth_bins < 3) stop_msg("Too few unique parent_tmrca breaks to define depth tiers.")
  
  depth_labels <- paste0("Tier_", seq_len(actual_depth_bins))
  
  edges_obs <- edges_obs %>%
    dplyr::mutate(
      depth_tier = cut(
        parent_tmrca,
        breaks = parent_age_breaks,
        include.lowest = TRUE,
        labels = depth_labels,
        ordered_result = TRUE
      ),
      depth_index = as.integer(depth_tier)
    )
  
  readr::write_csv(edges_obs, FILE_P5_EDGES)
  
  profile_obs <- edges_obs %>%
    dplyr::filter(!is.na(depth_tier)) %>%
    dplyr::group_by(rule, hg_group, depth_tier, depth_index) %>%
    dplyr::summarise(
      n_edges = dplyr::n(),
      same_region_rate = mean(same_region, na.rm = TRUE),
      cross_region_rate = mean(!same_region, na.rm = TRUE),
      median_edge_km = stats::median(edge_km, na.rm = TRUE),
      median_age_gap = stats::median(age_gap, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(rule, hg_group, depth_index)
  
  readr::write_csv(profile_obs, FILE_P5_PROFILE)
  
  network_diag <- edges_obs %>%
    dplyr::group_by(rule, hg_group) %>%
    dplyr::summarise(
      n_edges_output = dplyr::n(),
      median_edge_km = stats::median(edge_km, na.rm = TRUE),
      median_age_gap = stats::median(age_gap, na.rm = TRUE),
      median_parent_tmrca = stats::median(parent_tmrca, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      branch_tbl %>% dplyr::count(hg_group, name = "n_branches_input"),
      by = "hg_group"
    ) %>%
    dplyr::mutate(
      n_unmatched_branches = n_branches_input - 1L - n_edges_output
    ) %>%
    dplyr::arrange(rule, hg_group)
  
  readr::write_csv(network_diag, FILE_P5_NET_DIAG)
  
  compute_profile_stats <- function(profile_df) {
    wide <- profile_df %>%
      dplyr::select(rule, hg_group, depth_index, n_edges, same_region_rate) %>%
      tidyr::pivot_wider(
        names_from = hg_group,
        values_from = c(n_edges, same_region_rate)
      )
    
    if (!"same_region_rate_PH908" %in% names(wide)) wide$same_region_rate_PH908 <- NA_real_
    if (!"same_region_rate_R1a"   %in% names(wide)) wide$same_region_rate_R1a   <- NA_real_
    if (!"n_edges_PH908"          %in% names(wide)) wide$n_edges_PH908          <- NA_real_
    if (!"n_edges_R1a"            %in% names(wide)) wide$n_edges_R1a            <- NA_real_
    if (!"depth_index"            %in% names(wide)) wide$depth_index            <- NA_real_
    
    wide <- wide %>%
      dplyr::mutate(
        tier_shared = is.finite(same_region_rate_PH908) &
          is.finite(same_region_rate_R1a) &
          is.finite(n_edges_PH908) &
          is.finite(n_edges_R1a)
      ) %>%
      dplyr::filter(tier_shared)
    
    if (nrow(wide) == 0) {
      return(tibble::tibble(
        rule = unique(profile_df$rule)[1],
        n_shared_tiers = 0L,
        n_shared_tiers_meeting_min_edges = 0L,
        profile_auc_diff = NA_real_,
        profile_mean_diff = NA_real_,
        deep_weighted_diff = NA_real_,
        slope_diff = NA_real_,
        deep_minus_shallow_diff = NA_real_,
        all_shared_tiers_PH908_gt_R1a = NA,
        all_shared_tiers_min_edges_ok = NA
      ))
    }
    
    wide <- wide %>%
      dplyr::mutate(
        tier_diff = same_region_rate_PH908 - same_region_rate_R1a,
        min_edges_ok = n_edges_PH908 >= MIN_EDGES_PER_TIER & n_edges_R1a >= MIN_EDGES_PER_TIER,
        depth_weight = depth_index
      )
    
    shared_ok <- wide %>% dplyr::filter(min_edges_ok)
    analysis_df <- if (nrow(shared_ok) >= MIN_SHARED_DEPTH_TIERS) shared_ok else shared_ok[0, , drop = FALSE]
    
    profile_auc_val <- NA_real_
    if (nrow(analysis_df) >= 2) {
      x <- analysis_df$depth_index
      y <- analysis_df$tier_diff
      ord <- order(x)
      x <- x[ord]
      y <- y[ord]
      profile_auc_val <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
    } else if (nrow(analysis_df) == 1) {
      profile_auc_val <- analysis_df$tier_diff[[1]]
    }
    
    ph908_slope <- safe_lm_slope(analysis_df$depth_index, analysis_df$same_region_rate_PH908)
    r1a_slope <- safe_lm_slope(analysis_df$depth_index, analysis_df$same_region_rate_R1a)
    
    deep_minus_shallow <- NA_real_
    if (nrow(analysis_df) >= 2) {
      shallow_row <- analysis_df %>% dplyr::slice_min(depth_index, n = 1, with_ties = FALSE)
      deep_row <- analysis_df %>% dplyr::slice_max(depth_index, n = 1, with_ties = FALSE)
      deep_minus_shallow <- deep_row$tier_diff[[1]] - shallow_row$tier_diff[[1]]
    }
    
    tibble::tibble(
      rule = unique(profile_df$rule)[1],
      n_shared_tiers = nrow(wide),
      n_shared_tiers_meeting_min_edges = nrow(shared_ok),
      profile_auc_diff = profile_auc_val,
      profile_mean_diff = if (nrow(analysis_df) > 0) mean(analysis_df$tier_diff, na.rm = TRUE) else NA_real_,
      deep_weighted_diff = if (nrow(analysis_df) > 0 && sum(analysis_df$depth_weight) > 0) {
        sum(analysis_df$tier_diff * analysis_df$depth_weight, na.rm = TRUE) /
          sum(analysis_df$depth_weight, na.rm = TRUE)
      } else {
        NA_real_
      },
      slope_diff = if (is.finite(ph908_slope) && is.finite(r1a_slope)) ph908_slope - r1a_slope else NA_real_,
      deep_minus_shallow_diff = deep_minus_shallow,
      all_shared_tiers_PH908_gt_R1a = if (nrow(analysis_df) > 0) all(analysis_df$tier_diff > 0, na.rm = TRUE) else NA,
      all_shared_tiers_min_edges_ok = if (nrow(wide) > 0) all(wide$min_edges_ok, na.rm = TRUE) else NA
    )
  }
  
  profile_summary <- dplyr::bind_rows(lapply(split(profile_obs, profile_obs$rule), compute_profile_stats))
  readr::write_csv(profile_summary, FILE_P5_PROFILE_SUM)
  
  branch_perm_base <- branch_tbl %>%
    dplyr::filter(hg_group %in% c("PH908", "R1a")) %>%
    dplyr::mutate(
      age_bin = make_quantile_bins(tmrca_final, n_bins = AGE_BIN_N, prefix = "age"),
      size_bin = make_quantile_bins(log(branch_size + 1), n_bins = SIZE_BIN_N, prefix = "size"),
      strata = interaction(age_bin, size_bin, class, drop = TRUE)
    )
  
  strata_diag <- branch_perm_base %>%
    dplyr::count(strata, hg_group) %>%
    dplyr::arrange(strata, hg_group)
  
  readr::write_csv(strata_diag, FILE_P5_STRATA)
  
  null_list <- vector("list", N_PERM)
  
  for (b in seq_len(N_PERM)) {
    perm_group <- permute_within_strata(branch_perm_base$hg_group, branch_perm_base$strata)
    
    perm_df <- branch_perm_base %>%
      dplyr::mutate(hg_group = perm_group) %>%
      dplyr::select(
        terminal_snp_norm, branch_size,
        lat_c, lon_c, dominant_region,
        tmrca_final, class, hg_group
      )
    
    perm_edges <- build_all_rules_observed(perm_df)
    
    if (nrow(perm_edges) == 0) {
      null_list[[b]] <- tidyr::expand_grid(
        perm_id = b,
        rule = c("nearest_older", "nearest_older_window", "distance_age_cost"),
        stat_name = c(
          "profile_auc_diff",
          "profile_mean_diff",
          "deep_weighted_diff",
          "slope_diff",
          "deep_minus_shallow_diff"
        )
      ) %>%
        dplyr::mutate(stat_value = NA_real_)
      next
    }
    
    perm_edges <- perm_edges %>%
      dplyr::mutate(
        depth_tier = cut(
          parent_tmrca,
          breaks = parent_age_breaks,
          include.lowest = TRUE,
          labels = depth_labels,
          ordered_result = TRUE
        ),
        depth_index = as.integer(depth_tier)
      )
    
    perm_profile <- perm_edges %>%
      dplyr::filter(!is.na(depth_tier)) %>%
      dplyr::group_by(rule, hg_group, depth_tier, depth_index) %>%
      dplyr::summarise(
        n_edges = dplyr::n(),
        same_region_rate = mean(same_region, na.rm = TRUE),
        .groups = "drop"
      )
    
    perm_summary <- dplyr::bind_rows(lapply(split(perm_profile, perm_profile$rule), compute_profile_stats))
    
    perm_long <- perm_summary %>%
      dplyr::select(
        rule,
        profile_auc_diff,
        profile_mean_diff,
        deep_weighted_diff,
        slope_diff,
        deep_minus_shallow_diff
      ) %>%
      tidyr::pivot_longer(
        cols = c(
          profile_auc_diff,
          profile_mean_diff,
          deep_weighted_diff,
          slope_diff,
          deep_minus_shallow_diff
        ),
        names_to = "stat_name",
        values_to = "stat_value"
      ) %>%
      dplyr::mutate(perm_id = b) %>%
      dplyr::select(perm_id, rule, stat_name, stat_value)
    
    null_list[[b]] <- perm_long
  }
  
  null_tbl <- dplyr::bind_rows(null_list)
  readr::write_csv(null_tbl, FILE_P5_NULL)
  
  stat_names <- c(
    "profile_auc_diff",
    "profile_mean_diff",
    "deep_weighted_diff",
    "slope_diff",
    "deep_minus_shallow_diff"
  )
  
  formal_tests <- dplyr::bind_rows(lapply(
    split(profile_summary, profile_summary$rule),
    function(obs_rule_df) {
      rule_name <- unique(obs_rule_df$rule)[1]
      dplyr::bind_rows(lapply(stat_names, function(stat_i) {
        obs_val <- obs_rule_df[[stat_i]][1]
        null_vals <- null_tbl %>%
          dplyr::filter(rule == rule_name, stat_name == stat_i) %>%
          dplyr::pull(stat_value)
        ci <- ci95(null_vals)
        tibble::tibble(
          rule = rule_name,
          stat_name = stat_i,
          observed_value = obs_val,
          perm_mean = mean(null_vals, na.rm = TRUE),
          perm_sd = stats::sd(null_vals, na.rm = TRUE),
          ci_lower = ci[1],
          ci_upper = ci[2],
          p_two_sided = two_sided_p(null_vals, obs_val),
          p_one_sided_PH908_more_coherent = one_sided_p_greater(null_vals, obs_val)
        )
      }))
    }
  ))
  
  readr::write_csv(formal_tests, FILE_P5_TESTS)
  
  final_summary_tbl <- profile_summary %>%
    dplyr::left_join(
      formal_tests %>%
        dplyr::filter(stat_name == "profile_auc_diff") %>%
        dplyr::select(
          rule,
          profile_auc_diff_p_one_sided = p_one_sided_PH908_more_coherent,
          profile_auc_diff_p_two_sided = p_two_sided
        ),
      by = "rule"
    ) %>%
    dplyr::left_join(
      formal_tests %>%
        dplyr::filter(stat_name == "profile_mean_diff") %>%
        dplyr::select(
          rule,
          profile_mean_diff_p_one_sided = p_one_sided_PH908_more_coherent
        ),
      by = "rule"
    ) %>%
    dplyr::left_join(
      formal_tests %>%
        dplyr::filter(stat_name == "deep_weighted_diff") %>%
        dplyr::select(
          rule,
          deep_weighted_diff_p_one_sided = p_one_sided_PH908_more_coherent
        ),
      by = "rule"
    ) %>%
    dplyr::left_join(
      network_diag %>%
        dplyr::select(rule, hg_group, n_edges_output) %>%
        tidyr::pivot_wider(
          names_from = hg_group,
          values_from = n_edges_output,
          names_prefix = "n_edges_"
        ),
      by = "rule"
    ) %>%
    dplyr::mutate(
      rule_label = dplyr::case_when(
        rule == "nearest_older" ~ "Nearest older branch",
        rule == "nearest_older_window" ~ "Nearest older within age window",
        rule == "distance_age_cost" ~ "Distance + age cost",
        TRUE ~ rule
      ),
      support_label = dplyr::case_when(
        is.finite(profile_auc_diff) &
          profile_auc_diff > 0 &
          is.finite(profile_auc_diff_p_one_sided) &
          profile_auc_diff_p_one_sided < 0.05 ~ "Supported",
        TRUE ~ "Not supported"
      )
    ) %>%
    dplyr::select(
      rule,
      rule_label,
      n_edges_PH908 = n_edges_PH908,
      n_edges_R1a = n_edges_R1a,
      n_shared_tiers,
      n_shared_tiers_meeting_min_edges,
      profile_auc_diff,
      profile_auc_diff_p_one_sided,
      profile_auc_diff_p_two_sided,
      profile_mean_diff,
      profile_mean_diff_p_one_sided,
      deep_weighted_diff,
      deep_weighted_diff_p_one_sided,
      slope_diff,
      deep_minus_shallow_diff,
      all_shared_tiers_PH908_gt_R1a,
      all_shared_tiers_min_edges_ok,
      support_label
    ) %>%
    dplyr::arrange(factor(rule, levels = c(
      "nearest_older",
      "nearest_older_window",
      "distance_age_cost"
    )))
  
  readr::write_csv(final_summary_tbl, FILE_P5_FINAL)
  
  figuredata_p5_main <- profile_obs %>%
    dplyr::mutate(
      rule_label = dplyr::case_when(
        rule == "nearest_older" ~ "Nearest older branch",
        rule == "nearest_older_window" ~ "Nearest older within age window",
        rule == "distance_age_cost" ~ "Distance + age cost",
        TRUE ~ rule
      ),
      hg_group = factor(hg_group, levels = c("PH908", "R1a"))
    ) %>%
    dplyr::select(
      rule, rule_label, hg_group, depth_tier, depth_index,
      n_edges, same_region_rate, cross_region_rate,
      median_edge_km, median_age_gap
    )
  
  figuredata_p5_inset <- final_summary_tbl
  
  figuredata_p5_null <- null_tbl %>%
    dplyr::filter(stat_name == "profile_auc_diff") %>%
    dplyr::left_join(
      formal_tests %>%
        dplyr::filter(stat_name == "profile_auc_diff") %>%
        dplyr::select(rule, observed_value),
      by = "rule"
    ) %>%
    dplyr::mutate(
      rule_label = dplyr::case_when(
        rule == "nearest_older" ~ "Nearest older branch",
        rule == "nearest_older_window" ~ "Nearest older within age window",
        rule == "distance_age_cost" ~ "Distance + age cost",
        TRUE ~ rule
      )
    )
  
  table_p5_supp <- dplyr::bind_rows(
    age_diag %>%
      tidyr::pivot_longer(
        cols = c(n_branches, total_mass, n_relic, n_middle, n_founder, age_from_age_file, age_from_dataset, pct_age_from_age_file, pct_age_from_dataset, median_branch_size, max_branch_size),
        names_to = "metric",
        values_to = "value"
      ) %>%
      tidyr::pivot_wider(names_from = hg_group, values_from = value),
    network_diag %>%
      dplyr::transmute(
        metric = paste0(rule, "__", hg_group, "__network"),
        PH908 = dplyr::if_else(hg_group == "PH908", median_edge_km, NA_real_),
        R1a = dplyr::if_else(hg_group == "R1a", median_edge_km, NA_real_)
      ),
    final_summary_tbl %>%
      dplyr::transmute(
        metric = paste0(rule, "__", c("profile_auc_diff", "profile_mean_diff", "deep_weighted_diff")),
        PH908 = NA_real_,
        R1a = NA_real_
      ) %>% dplyr::slice(0)
  )
  
  readr::write_csv(figuredata_p5_main, FILE_P5_FIG_MAIN)
  readr::write_csv(figuredata_p5_inset, FILE_P5_FIG_INSET)
  readr::write_csv(figuredata_p5_null, FILE_P5_FIG_NULL)
  readr::write_csv(table_p5_supp, FILE_P5_SUPP_TABLE)
  
  row_nearest <- final_summary_tbl %>% dplyr::filter(rule == "nearest_older")
  row_window <- final_summary_tbl %>% dplyr::filter(rule == "nearest_older_window")
  row_cost <- final_summary_tbl %>% dplyr::filter(rule == "distance_age_cost")
  
  net_nearest_ph <- network_diag %>% dplyr::filter(rule == "nearest_older", hg_group == "PH908")
  net_nearest_r1a <- network_diag %>% dplyr::filter(rule == "nearest_older", hg_group == "R1a")
  net_window_ph <- network_diag %>% dplyr::filter(rule == "nearest_older_window", hg_group == "PH908")
  net_window_r1a <- network_diag %>% dplyr::filter(rule == "nearest_older_window", hg_group == "R1a")
  net_cost_ph <- network_diag %>% dplyr::filter(rule == "distance_age_cost", hg_group == "PH908")
  net_cost_r1a <- network_diag %>% dplyr::filter(rule == "distance_age_cost", hg_group == "R1a")
  
  assert_equal(age_diag$n_branches[age_diag$hg_group == "PH908"][[1]], 79, "PH908 branch count drift in Pillar 5.")
  assert_equal(age_diag$n_branches[age_diag$hg_group == "R1a"][[1]], 81, "R1a branch count drift in Pillar 5.")
  
  assert_equal(net_nearest_ph$n_edges_output[[1]], 78, "Nearest older PH908 edge count drift.")
  assert_equal(net_nearest_r1a$n_edges_output[[1]], 80, "Nearest older R1a edge count drift.")
  assert_equal(net_window_ph$n_edges_output[[1]], 78, "Nearest older window PH908 edge count drift.")
  assert_equal(net_window_r1a$n_edges_output[[1]], 80, "Nearest older window R1a edge count drift.")
  assert_equal(net_cost_ph$n_edges_output[[1]], 78, "Distance-age-cost PH908 edge count drift.")
  assert_equal(net_cost_r1a$n_edges_output[[1]], 80, "Distance-age-cost R1a edge count drift.")
  
  assert_close(row_nearest$profile_auc_diff[[1]], -0.004, tol = 0.01, "Nearest older profile_auc_diff drift.")
  assert_close(row_window$profile_auc_diff[[1]], 0.092, tol = 0.01, "Nearest older window profile_auc_diff drift.")
  assert_close(row_cost$profile_auc_diff[[1]], 0.178, tol = 0.01, "Distance-age-cost profile_auc_diff drift.")
  
  assert_close(row_window$profile_mean_diff[[1]], 0.092, tol = 0.01, "Nearest older window profile_mean_diff drift.")
  assert_close(row_window$deep_weighted_diff[[1]], 0.037, tol = 0.015, "Nearest older window deep_weighted_diff drift.")
  assert_close(row_cost$profile_mean_diff[[1]], 0.178, tol = 0.01, "Distance-age-cost profile_mean_diff drift.")
  assert_close(row_cost$deep_weighted_diff[[1]], 0.187, tol = 0.015, "Distance-age-cost deep_weighted_diff drift.")
  assert_close(row_nearest$deep_weighted_diff[[1]], -0.069, tol = 0.02, "Nearest older deep_weighted_diff drift.")
  
  if (N_PERM >= 3000L) {
    assert_close(row_window$profile_auc_diff_p_one_sided[[1]], 0.020, tol = 0.02, "Nearest older window profile_auc_diff p-value drift.")
    assert_close(row_cost$profile_auc_diff_p_one_sided[[1]], 0.004, tol = 0.01, "Distance-age-cost profile_auc_diff p-value drift.")
  }
  
  git_commit <- safe_git_commit()
  
  results_sentence <- paste0(
    "Under nearest_older, the profile-wide contrast was effectively null (profile_auc_diff = ",
    sprintf("%.3f", row_nearest$profile_auc_diff[[1]]),
    "). Under nearest_older_window, the contrast became positive (",
    sprintf("%.3f", row_window$profile_auc_diff[[1]]),
    "), and under distance_age_cost it was strongest (",
    sprintf("%.3f", row_cost$profile_auc_diff[[1]]),
    ")."
  )
  
  caption_sentence <- paste0(
    "Same-region linkage rate is summarized across parent-age depth tiers under nearest_older, nearest_older_window, and distance_age_cost. The main rule summary should emphasize that the two age-aware rules support the PH908 depth-coherence advantage, whereas nearest_older does not."
  )
  
  facts_lines <- c(
    paste0("script: ", SCRIPT_NAME),
    paste0("timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("seed: ", RNG_SEED),
    paste0("config_file: ", path_for_facts(CONFIG_FILE)),
    "project_root: .",
    paste0("git_commit: ", git_commit),
    "inputs:",
    paste0("  - ", path_for_facts(DATA_FILE)),
    paste0("  - ", path_for_facts(AGE_FILE)),
    paste0("ph908_gating_mode: ", PH908_GATING_MODE),
    paste0("n_perm: ", N_PERM),
    paste0("age_bin_n: ", AGE_BIN_N),
    paste0("size_bin_n: ", SIZE_BIN_N),
    paste0("depth_bin_n_requested: ", DEPTH_BIN_N),
    paste0("depth_bin_n_achieved: ", actual_depth_bins),
    paste0("min_shared_depth_tiers: ", MIN_SHARED_DEPTH_TIERS),
    paste0("min_edges_per_tier: ", MIN_EDGES_PER_TIER),
    "",
    "branch_structure:",
    paste0("  PH908_branches: ", age_diag$n_branches[age_diag$hg_group == "PH908"][[1]]),
    paste0("  PH908_total_mass: ", age_diag$total_mass[age_diag$hg_group == "PH908"][[1]]),
    paste0("  R1a_branches: ", age_diag$n_branches[age_diag$hg_group == "R1a"][[1]]),
    paste0("  R1a_total_mass: ", age_diag$total_mass[age_diag$hg_group == "R1a"][[1]]),
    "",
    "network_edges:",
    paste0("  nearest_older: ", net_nearest_ph$n_edges_output[[1]], "/", net_nearest_r1a$n_edges_output[[1]]),
    paste0("  nearest_older_window: ", net_window_ph$n_edges_output[[1]], "/", net_window_r1a$n_edges_output[[1]]),
    paste0("  distance_age_cost: ", net_cost_ph$n_edges_output[[1]], "/", net_cost_r1a$n_edges_output[[1]]),
    "",
    "observed_profile_summary:",
    paste0("  nearest_older_profile_auc_diff: ", sprintf("%.3f", row_nearest$profile_auc_diff[[1]])),
    paste0("  nearest_older_window_profile_auc_diff: ", sprintf("%.3f", row_window$profile_auc_diff[[1]])),
    paste0("  distance_age_cost_profile_auc_diff: ", sprintf("%.3f", row_cost$profile_auc_diff[[1]])),
    paste0("  nearest_older_window_profile_mean_diff: ", sprintf("%.3f", row_window$profile_mean_diff[[1]])),
    paste0("  nearest_older_window_deep_weighted_diff: ", sprintf("%.3f", row_window$deep_weighted_diff[[1]])),
    paste0("  distance_age_cost_profile_mean_diff: ", sprintf("%.3f", row_cost$profile_mean_diff[[1]])),
    paste0("  distance_age_cost_deep_weighted_diff: ", sprintf("%.3f", row_cost$deep_weighted_diff[[1]])),
    paste0("  nearest_older_deep_weighted_diff: ", sprintf("%.3f", row_nearest$deep_weighted_diff[[1]])),
    "",
    "rule_pattern: nearest_older is non-supporting; nearest_older_window supports a positive PH908 coherence contrast; distance_age_cost is strongest.",
    "rule_interpretation_summary:",
    "  nearest_older: non-supporting",
    "  nearest_older_window: supportive",
    "  distance_age_cost: strongest support",
    paste0("results_sentence: ", results_sentence),
    paste0("caption_sentence: ", caption_sentence),
    "scope_note: This is a comparative inferred branch-link depth-profile analysis rather than a literal parent-child phylogeny or migration-route reconstruction.",
    "note: permutation-based consistency checks were run at the configured depth.",
    "output_files:",
    paste0("  - ", path_for_facts(FILE_P5_RESULTS)),
    paste0("  - ", path_for_facts(FILE_P5_FIG_MAIN)),
    paste0("  - ", path_for_facts(FILE_P5_FIG_INSET)),
    paste0("  - ", path_for_facts(FILE_P5_FIG_NULL)),
    paste0("  - ", path_for_facts(FILE_P5_SUPP_TABLE)),
    paste0("  - ", path_for_facts(FILE_P5_FACTS)),
    paste0("  - ", path_for_facts(FILE_P5_WARNINGS)),
    paste0("  - ", path_for_facts(FILE_P5_FILTER_DIAG)),
    paste0("  - ", path_for_facts(FILE_P5_AGE_DIAG)),
    paste0("  - ", path_for_facts(FILE_P5_BRANCH)),
    paste0("  - ", path_for_facts(FILE_P5_EDGES)),
    paste0("  - ", path_for_facts(FILE_P5_PROFILE)),
    paste0("  - ", path_for_facts(FILE_P5_PROFILE_SUM)),
    paste0("  - ", path_for_facts(FILE_P5_TESTS)),
    paste0("  - ", path_for_facts(FILE_P5_NULL)),
    paste0("  - ", path_for_facts(FILE_P5_FINAL)),
    paste0("  - ", path_for_facts(FILE_P5_NET_DIAG)),
    paste0("  - ", path_for_facts(FILE_P5_STRATA))
  )
  
  writeLines(facts_lines, FILE_P5_FACTS)
  
  results_p5 <- list(
    metadata = list(
      script = SCRIPT_NAME,
      seed = RNG_SEED,
      config_file = path_for_facts(CONFIG_FILE),
      inputs = c(path_for_facts(DATA_FILE), path_for_facts(AGE_FILE)),
      ph908_gating_mode = PH908_GATING_MODE,
      max_age_gap_window = MAX_AGE_GAP_WINDOW,
      lambda_age_cost = LAMBDA_AGE_COST,
      n_perm = N_PERM,
      age_bin_n = AGE_BIN_N,
      size_bin_n = SIZE_BIN_N,
      depth_bin_n_requested = DEPTH_BIN_N,
      depth_bin_n_achieved = actual_depth_bins,
      min_shared_depth_tiers = MIN_SHARED_DEPTH_TIERS,
      min_edges_per_tier = MIN_EDGES_PER_TIER,
      git_commit = git_commit
    ),
    filter_diag = filter_diag,
    age_diag = age_diag,
    branch_tbl = branch_tbl,
    edges_obs = edges_obs,
    profile_obs = profile_obs,
    profile_summary = profile_summary,
    formal_tests = formal_tests,
    final_summary_tbl = final_summary_tbl,
    network_diag = network_diag,
    strata_diag = strata_diag,
    null_tbl = null_tbl,
    figuredata_main = figuredata_p5_main,
    figuredata_inset = figuredata_p5_inset,
    figuredata_null = figuredata_p5_null,
    table_p5_supp = table_p5_supp
  )
  
  saveRDS(results_p5, FILE_P5_RESULTS)
  
}, warning = function(w) {
  warning_log <<- c(warning_log, conditionMessage(w))
  invokeRestart("muffleWarning")
})

write_warning_log(warning_log, FILE_P5_WARNINGS)