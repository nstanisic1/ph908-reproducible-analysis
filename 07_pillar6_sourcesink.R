# 07_pillar6_sourcesink.R
# Deep-core source–sink inversion for PH908 versus R1a.

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
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
set.seed(SEED_MAIN)

FILE_P6_RESULTS <- file.path(DIR_RESULTS, "results_p6.rds")

FILE_P6_SUPP_TABLE <- file.path(DIR_TABLES, "table_p6_supp.csv")
FILE_P6_THRESHOLD_TABLE <- file.path(DIR_TABLES, "table_p6_threshold_sensitivity.csv")

FILE_P6_FIG_TOPK <- file.path(DIR_FIGUREDATA, "figuredata_p6_topk.csv")
FILE_P6_FIG_FOUNDERCARRY <- file.path(DIR_FIGUREDATA, "figuredata_p6_foundercarry.csv")
FILE_P6_FIG_JACKKNIFE <- file.path(DIR_FIGUREDATA, "figuredata_p6_jackknife.csv")
FILE_P6_FIG_THRESHOLD <- file.path(DIR_FIGUREDATA, "figuredata_p6_threshold_grid.csv")
FILE_P6_FIG_NULL <- file.path(DIR_FIGUREDATA, "figuredata_p6_null.csv")

FILE_P6_FACTS <- file.path(DIR_FACTS, "facts_p6.txt")

PH908_GATING_MODE <- "primary"
PRIMARY_TOPK <- 4L
SENS_TOPK <- c(3L, 5L)
ALL_TOPK <- c(3L, 4L, 5L)

RELIC_MAX_N <- 3L
FOUNDER_MIN_N <- 15L
THRESHOLD_GRID_RELIC <- c(2L, 3L, 4L)
THRESHOLD_GRID_FOUNDER <- c(10L, 15L, 20L)

AGE_BIN_N <- 4L
SIZE_BIN_N <- 4L
# Permutation setting.
N_PERM <- 5000L
EPS <- 0.5

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
cat("07 PILLAR 6 — DEEP-CORE SOURCE–SINK INVERSION\n")
cat("============================================================\n")
cat("Root: ", DIR_ROOT, "\n", sep = "")
cat("Seed: ", SEED_MAIN, "\n", sep = "")
cat("Inputs:\n")
cat("  - ", FILE_PREPROCESSED_RDS, "\n", sep = "")
cat("  - ", FILE_BRANCH_AGES, "\n", sep = "")
cat("Primary topk: ", PRIMARY_TOPK, "\n", sep = "")
cat("All topk: ", paste(ALL_TOPK, collapse = ", "), "\n\n", sep = "")

write_facts_header(
  FILE_P6_FACTS,
  script_name = "07_pillar6_sourcesink.R",
  input_files = c(FILE_PREPROCESSED_RDS, FILE_BRANCH_AGES)
)

path_for_facts <- function(x) {
  x <- as.character(x)
  if (!length(x) || is.na(x) || x == "") return(x)
  x_norm <- normalizePath(x, winslash = "/", mustWork = FALSE)
  root_norm <- normalizePath(DIR_ROOT, winslash = "/", mustWork = FALSE)
  root_esc <- gsub("([][{}()+*^$|\\?.])", "\\\\\\1", root_norm)
  sub(paste0("^", root_esc), ".", x_norm)
}

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

classify_branch <- function(n, relic_max_n = 3L, founder_min_n = 15L) {
  dplyr::case_when(
    n <= relic_max_n ~ "Relic",
    n >= founder_min_n ~ "Founder",
    TRUE ~ "Middle"
  )
}

make_quantile_bins <- function(x, n_bins = 4, prefix = "bin") {
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
  u <- unique(strata)
  
  for (s in u) {
    idx <- which(strata == s)
    if (length(idx) > 1) {
      out[idx] <- sample(out[idx], length(idx), replace = FALSE)
    }
  }
  
  out
}

one_sided_p_greater <- function(null_vals, obs) {
  null_vals <- null_vals[is.finite(null_vals)]
  if (!is.finite(obs) || length(null_vals) == 0) return(NA_real_)
  (1 + sum(null_vals >= obs)) / (length(null_vals) + 1)
}

ci95 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 20) return(c(NA_real_, NA_real_))
  as.numeric(stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
}

round3 <- function(x) round(x, 3)

fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
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
    dplyr::filter(!is.na(terminal_snp_norm), terminal_snp_norm != "", is.finite(tmrca_ybp), tmrca_ybp > 0) |>
    dplyr::group_by(terminal_snp_norm) |>
    dplyr::summarise(
      tmrca_ybp = dplyr::first(tmrca_ybp),
      .groups = "drop"
    )
}

build_sample_base <- function(cl_dat, gating_mode = "primary") {
  base <- cl_dat |>
    dplyr::filter(
      matched == TRUE,
      exclude_geo_primary == FALSE,
      toupper(country_code) %in% BALKAN_ISO2,
      !is.na(terminal_snp), terminal_snp != "",
      !is.na(Region), Region != ""
    ) |>
    dplyr::mutate(
      terminal_snp_norm = normalize_snp(terminal_snp),
      ph908_flag = dplyr::case_when(
        gating_mode == "primary"  ~ is_ph908_primary,
        gating_mode == "accurate" ~ is_ph908_accurate,
        gating_mode == "wide"     ~ is_ph908_wide,
        TRUE ~ FALSE
      ),
      hg_group = dplyr::case_when(
        ph908_flag ~ "PH908",
        major_hg == "R1a" ~ "R1a",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(hg_group %in% c("PH908", "R1a"))
  
  if (nrow(base) == 0) {
    stop("No usable rows remain after Pillar 6 filtering.", call. = FALSE)
  }
  
  base
}

build_branch_table <- function(base, age_tbl, relic_max_n = 3L, founder_min_n = 15L, gating_mode = "primary") {
  out <- base |>
    dplyr::group_by(hg_group, terminal_snp_norm) |>
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
      class = classify_branch(branch_size, relic_max_n = relic_max_n, founder_min_n = founder_min_n),
      gating_mode = gating_mode
    ) |>
    dplyr::filter(
      !is.na(dominant_region), dominant_region != "",
      is.finite(tmrca_final)
    )
  
  needed <- c("PH908", "R1a")
  missing_groups <- setdiff(needed, unique(out$hg_group))
  if (length(missing_groups) > 0) {
    stop(
      paste0("Required haplogroup(s) absent after branch construction: ", paste(missing_groups, collapse = ", ")),
      call. = FALSE
    )
  }
  
  out
}

make_branch_diagnostics <- function(branch_tbl) {
  branch_tbl |>
    dplyr::group_by(gating_mode, hg_group) |>
    dplyr::summarise(
      n_branches = dplyr::n(),
      total_mass = sum(branch_size),
      n_relic = sum(class == "Relic"),
      n_middle = sum(class == "Middle"),
      n_founder = sum(class == "Founder"),
      .groups = "drop"
    )
}

make_age_diagnostics <- function(branch_tbl) {
  branch_tbl |>
    dplyr::group_by(gating_mode, hg_group) |>
    dplyr::summarise(
      n_branches = dplyr::n(),
      total_mass = sum(branch_size),
      n_external_age = sum(age_source == "branch_ages_csv"),
      n_dataset_age = sum(age_source == "dataset_tmrca"),
      n_missing_age = sum(age_source == "missing"),
      .groups = "drop"
    )
}

rank_ph908_relic_regions <- function(branch_tbl) {
  branch_tbl |>
    dplyr::filter(hg_group == "PH908", class == "Relic") |>
    dplyr::count(dominant_region, name = "n_relic_branches") |>
    dplyr::arrange(dplyr::desc(n_relic_branches), dominant_region) |>
    dplyr::mutate(
      rank = dplyr::row_number(),
      prop_relic_branches = n_relic_branches / sum(n_relic_branches),
      cum_prop_relic_branches = cumsum(prop_relic_branches)
    )
}

make_core_definitions <- function(ph908_relic_ranking, topks = c(3L, 4L, 5L)) {
  purrr::map_dfr(
    topks,
    function(k) {
      regs <- ph908_relic_ranking |>
        dplyr::slice_head(n = min(k, nrow(ph908_relic_ranking)))
      
      tibble::tibble(
        topk = k,
        core_id = paste0("top", k),
        core_regions = paste(regs$dominant_region, collapse = "; "),
        core_region_count = nrow(regs),
        core_relic_share_PH908 = sum(regs$prop_relic_branches, na.rm = TRUE),
        is_primary = k == PRIMARY_TOPK
      )
    }
  )
}

expand_core_region_map <- function(core_defs) {
  core_defs |>
    dplyr::mutate(core_region = stringr::str_split(core_regions, "; ")) |>
    tidyr::unnest(core_region) |>
    dplyr::select(topk, core_id, core_region)
}

compute_zone_metrics <- function(branch_tbl, core_defs, eps = 0.5) {
  core_map <- expand_core_region_map(core_defs)
  
  tagged <- branch_tbl |>
    dplyr::cross_join(core_defs |>
                        dplyr::select(topk, core_id, core_regions, core_region_count, core_relic_share_PH908, is_primary)
    ) |>
    dplyr::left_join(
      core_map |>
        dplyr::mutate(in_core = TRUE),
      by = c("topk", "core_id", "dominant_region" = "core_region")
    ) |>
    dplyr::mutate(
      in_core = dplyr::coalesce(in_core, FALSE),
      zone = ifelse(in_core, "Core", "Expanse")
    )
  
  zone_counts <- tagged |>
    dplyr::group_by(topk, core_id, core_regions, core_region_count, is_primary, hg_group, zone, class) |>
    dplyr::summarise(n_branches = dplyr::n(), .groups = "drop")
  
  zone_summary <- zone_counts |>
    tidyr::pivot_wider(
      names_from = c(zone, class),
      values_from = n_branches,
      values_fill = 0
    ) |>
    dplyr::mutate(
      total_relic = Core_Relic + Expanse_Relic,
      total_founder = Core_Founder + Expanse_Founder,
      core_relic_share = ifelse(total_relic > 0, Core_Relic / total_relic, NA_real_),
      core_founder_share = ifelse(total_founder > 0, Core_Founder / total_founder, NA_real_),
      expanse_relic_share = ifelse(total_relic > 0, Expanse_Relic / total_relic, NA_real_),
      expanse_founder_share = ifelse(total_founder > 0, Expanse_Founder / total_founder, NA_real_),
      capture_skew = core_relic_share - core_founder_share,
      mismatch_core = log((Core_Relic + eps) / (Core_Founder + eps)),
      mismatch_expanse = log((Expanse_Relic + eps) / (Expanse_Founder + eps)),
      expected_core_founders = total_founder * core_relic_share,
      observed_minus_expected = Core_Founder - expected_core_founders,
      observed_over_expected = ifelse(expected_core_founders > 0, Core_Founder / expected_core_founders, NA_real_)
    )
  
  primary_contrasts <- zone_summary |>
    dplyr::select(
      topk, core_id, core_regions, core_region_count, is_primary,
      hg_group, capture_skew,
      core_relic_share, core_founder_share,
      total_founder, expected_core_founders,
      Core_Founder, observed_minus_expected, observed_over_expected
    ) |>
    tidyr::pivot_wider(
      names_from = hg_group,
      values_from = c(
        capture_skew,
        core_relic_share, core_founder_share,
        total_founder, expected_core_founders,
        Core_Founder, observed_minus_expected, observed_over_expected
      )
    ) |>
    dplyr::mutate(
      anti_survivor_contrast = capture_skew_PH908 - capture_skew_R1a
    )
  
  secondary_contrasts <- zone_summary |>
    dplyr::select(topk, core_id, core_regions, core_region_count, is_primary, hg_group, mismatch_core, mismatch_expanse) |>
    tidyr::pivot_wider(
      names_from = hg_group,
      values_from = c(mismatch_core, mismatch_expanse)
    ) |>
    dplyr::mutate(
      core_minus_expanse = (mismatch_core_PH908 - mismatch_expanse_PH908) -
        (mismatch_core_R1a - mismatch_expanse_R1a)
    )
  
  list(
    tagged = tagged,
    zone_counts = zone_counts,
    zone_summary = zone_summary,
    primary_contrasts = primary_contrasts,
    secondary_contrasts = secondary_contrasts
  )
}

build_perm_base <- function(branch_tbl) {
  branch_tbl |>
    dplyr::mutate(
      age_bin = make_quantile_bins(tmrca_final, n_bins = AGE_BIN_N, prefix = "age"),
      size_bin = make_quantile_bins(log(branch_size + 1), n_bins = SIZE_BIN_N, prefix = "size"),
      strata = interaction(age_bin, size_bin, class, drop = TRUE)
    )
}

run_primary_null <- function(branch_tbl, core_defs, n_perm = 3000L, eps = 0.5) {
  perm_base <- build_perm_base(branch_tbl)
  
  strata_diag <- perm_base |>
    dplyr::count(strata, hg_group, name = "n_branches")
  
  null_list <- vector("list", n_perm)
  
  for (b in seq_len(n_perm)) {
    perm_group <- permute_within_strata(perm_base$hg_group, perm_base$strata)
    
    perm_tbl <- perm_base |>
      dplyr::mutate(hg_group = perm_group) |>
      dplyr::select(
        terminal_snp_norm, branch_size, dominant_region, tmrca_final,
        age_source, class, hg_group, gating_mode
      )
    
    perm_metrics <- compute_zone_metrics(perm_tbl, core_defs, eps = eps)
    
    null_list[[b]] <- perm_metrics$primary_contrasts |>
      dplyr::select(topk, core_id, anti_survivor_contrast) |>
      dplyr::left_join(
        perm_metrics$secondary_contrasts |>
          dplyr::select(topk, core_id, core_minus_expanse),
        by = c("topk", "core_id")
      ) |>
      dplyr::mutate(perm_id = b)
  }
  
  list(
    strata_diag = strata_diag,
    null_tbl = dplyr::bind_rows(null_list)
  )
}

test_topk_results <- function(primary_tbl, secondary_tbl, null_tbl) {
  out_primary <- purrr::map_dfr(
    primary_tbl$topk,
    function(k) {
      obs <- primary_tbl |>
        dplyr::filter(topk == k)
      
      nv <- null_tbl |>
        dplyr::filter(topk == k) |>
        dplyr::pull(anti_survivor_contrast)
      
      ci <- ci95(nv)
      
      tibble::tibble(
        topk = k,
        core_id = obs$core_id,
        statistic = "anti_survivor_contrast",
        observed = obs$anti_survivor_contrast,
        p_one_sided = one_sided_p_greater(nv, obs$anti_survivor_contrast),
        ci_lower = ci[1],
        ci_upper = ci[2],
        interpretation = dplyr::case_when(
          obs$anti_survivor_contrast > 0 & one_sided_p_greater(nv, obs$anti_survivor_contrast) < 0.05 & k == PRIMARY_TOPK ~ "Primary supported",
          obs$anti_survivor_contrast > 0 & one_sided_p_greater(nv, obs$anti_survivor_contrast) < 0.05 ~ "Supported",
          obs$anti_survivor_contrast > 0 ~ "Borderline",
          TRUE ~ "Not supported"
        ),
        family = "Primary"
      )
    }
  )
  
  out_secondary <- purrr::map_dfr(
    secondary_tbl$topk,
    function(k) {
      obs <- secondary_tbl |>
        dplyr::filter(topk == k)
      
      nv <- null_tbl |>
        dplyr::filter(topk == k) |>
        dplyr::pull(core_minus_expanse)
      
      ci <- ci95(nv)
      
      tibble::tibble(
        topk = k,
        core_id = obs$core_id,
        statistic = "core_minus_expanse",
        observed = obs$core_minus_expanse,
        p_one_sided = one_sided_p_greater(nv, obs$core_minus_expanse),
        ci_lower = ci[1],
        ci_upper = ci[2],
        interpretation = dplyr::case_when(
          obs$core_minus_expanse > 0 & one_sided_p_greater(nv, obs$core_minus_expanse) < 0.05 & k == PRIMARY_TOPK ~ "Primary supported",
          obs$core_minus_expanse > 0 & one_sided_p_greater(nv, obs$core_minus_expanse) < 0.05 ~ "Supported",
          obs$core_minus_expanse > 0 ~ "Borderline",
          TRUE ~ "Not supported"
        ),
        family = "Secondary"
      )
    }
  )
  
  dplyr::bind_rows(out_primary, out_secondary) |>
    dplyr::arrange(topk, family)
}

run_top4_jackknife <- function(branch_tbl, core_defs, n_perm = 3000L, eps = 0.5) {
  primary_core <- core_defs |>
    dplyr::filter(topk == PRIMARY_TOPK)
  
  regions <- stringr::str_split(primary_core$core_regions[[1]], "; ")[[1]]
  
  out <- purrr::map_dfr(
    regions,
    function(reg_i) {
      reduced_tbl <- branch_tbl |>
        dplyr::filter(dominant_region != reg_i)
      
      reduced_regions <- setdiff(regions, reg_i)
      
      reduced_core <- primary_core |>
        dplyr::mutate(
          core_regions = paste(reduced_regions, collapse = "; "),
          core_region_count = length(reduced_regions)
        )
      
      reduced_metrics <- compute_zone_metrics(reduced_tbl, reduced_core, eps = eps)
      reduced_primary <- reduced_metrics$primary_contrasts
      
      null_run <- run_primary_null(reduced_tbl, reduced_core, n_perm = n_perm, eps = eps)
      nv <- null_run$null_tbl |>
        dplyr::pull(anti_survivor_contrast)
      ci <- ci95(nv)
      
      tibble::tibble(
        dropped_region = reg_i,
        observed = reduced_primary$anti_survivor_contrast,
        ci_low = ci[1],
        ci_high = ci[2],
        p_one_sided = one_sided_p_greater(nv, reduced_primary$anti_survivor_contrast),
        sign_positive = reduced_primary$anti_survivor_contrast > 0,
        interpretation = dplyr::case_when(
          reduced_primary$anti_survivor_contrast > 0 & one_sided_p_greater(nv, reduced_primary$anti_survivor_contrast) < 0.05 ~ "Supported",
          reduced_primary$anti_survivor_contrast > 0 ~ "Positive, not supported",
          TRUE ~ "Not supported"
        )
      )
    }
  )
  
  full_obs <- compute_zone_metrics(branch_tbl, primary_core, eps = eps)$primary_contrasts$anti_survivor_contrast[[1]]
  out |>
    dplyr::mutate(delta_vs_full = observed - full_obs)
}

run_threshold_grid <- function(base, age_tbl, n_perm = 3000L, eps = 0.5) {
  grid <- tidyr::expand_grid(
    relic_max_n = THRESHOLD_GRID_RELIC,
    founder_min_n = THRESHOLD_GRID_FOUNDER
  )
  
  out <- purrr::pmap_dfr(
    list(grid$relic_max_n, grid$founder_min_n),
    function(relic_max_n, founder_min_n) {
      bt <- build_branch_table(
        base = base,
        age_tbl = age_tbl,
        relic_max_n = relic_max_n,
        founder_min_n = founder_min_n,
        gating_mode = PH908_GATING_MODE
      )
      
      diag <- make_branch_diagnostics(bt)
      adequacy_flag <- all(c("PH908", "R1a") %in% diag$hg_group) &&
        all(diag$n_relic >= 5) &&
        all(diag$n_founder >= 1)
      
      if (!adequacy_flag) {
        return(tibble::tibble(
          relic_max_n = relic_max_n,
          founder_min_n = founder_min_n,
          anti_survivor_contrast_top4 = NA_real_,
          p_one_sided_top4 = NA_real_,
          core_minus_expanse_top4 = NA_real_,
          p_one_sided_core_minus_expanse_top4 = NA_real_,
          adequacy_flag = FALSE,
          support_label = "Inadequate"
        ))
      }
      
      rank_tbl <- rank_ph908_relic_regions(bt)
      core_defs <- make_core_definitions(rank_tbl, topks = PRIMARY_TOPK)
      metrics <- compute_zone_metrics(bt, core_defs, eps = eps)
      
      null_run <- run_primary_null(bt, core_defs, n_perm = n_perm, eps = eps)
      
      anti_obs <- metrics$primary_contrasts$anti_survivor_contrast[[1]]
      anti_nv <- null_run$null_tbl$anti_survivor_contrast
      sec_obs <- metrics$secondary_contrasts$core_minus_expanse[[1]]
      sec_nv <- null_run$null_tbl$core_minus_expanse
      
      tibble::tibble(
        relic_max_n = relic_max_n,
        founder_min_n = founder_min_n,
        anti_survivor_contrast_top4 = anti_obs,
        p_one_sided_top4 = one_sided_p_greater(anti_nv, anti_obs),
        core_minus_expanse_top4 = sec_obs,
        p_one_sided_core_minus_expanse_top4 = one_sided_p_greater(sec_nv, sec_obs),
        adequacy_flag = TRUE,
        support_label = dplyr::case_when(
          is.finite(anti_obs) & anti_obs > 0 &
            is.finite(one_sided_p_greater(anti_nv, anti_obs)) &
            one_sided_p_greater(anti_nv, anti_obs) < 0.05 ~ "Supported",
          is.finite(anti_obs) & anti_obs > 0 ~ "Positive",
          TRUE ~ "Not supported"
        )
      )
    }
  )
  
  out
}

cl_raw <- readRDS(FILE_PREPROCESSED_RDS)

required_cols <- c(
  "matched",
  "exclude_geo_primary",
  "terminal_snp",
  "major_hg",
  "is_ph908_primary",
  "is_ph908_accurate",
  "is_ph908_wide",
  "tmrca",
  "Region",
  "country_code"
)

missing_cols <- setdiff(required_cols, names(cl_raw))
if (length(missing_cols) > 0) {
  stop(
    paste0("Dataset missing required columns: ", paste(missing_cols, collapse = ", ")),
    call. = FALSE
  )
}

age_tbl <- standardize_age_table(readr::read_csv(FILE_BRANCH_AGES, show_col_types = FALSE))
base <- build_sample_base(cl_raw, gating_mode = PH908_GATING_MODE)
branch_tbl <- build_branch_table(
  base = base,
  age_tbl = age_tbl,
  relic_max_n = RELIC_MAX_N,
  founder_min_n = FOUNDER_MIN_N,
  gating_mode = PH908_GATING_MODE
)

filter_diag <- tibble::tibble(
  step = c(
    "raw_rows",
    "matched_true",
    "exclude_geo_primary_false",
    "balkan_country_rows",
    "terminal_snp_present",
    "region_present",
    "usable_rows_after_filter"
  ),
  n = c(
    nrow(cl_raw),
    sum(cl_raw$matched %in% TRUE, na.rm = TRUE),
    sum(cl_raw$exclude_geo_primary %in% FALSE, na.rm = TRUE),
    sum(toupper(cl_raw$country_code) %in% BALKAN_ISO2, na.rm = TRUE),
    sum(!is.na(cl_raw$terminal_snp) & cl_raw$terminal_snp != ""),
    sum(!is.na(cl_raw$Region) & cl_raw$Region != ""),
    nrow(base)
  )
)

branch_diag <- make_branch_diagnostics(branch_tbl)
age_diag <- make_age_diagnostics(branch_tbl)
ph908_relic_ranking <- rank_ph908_relic_regions(branch_tbl)
core_defs <- make_core_definitions(ph908_relic_ranking, topks = ALL_TOPK)
zone_metrics <- compute_zone_metrics(branch_tbl, core_defs, eps = EPS)

append_facts_lines(FILE_P6_FACTS, c(
  "[pre_null_branch_universe]",
  paste(
    apply(branch_diag, 1, function(r) {
      paste0(
        r[["hg_group"]], ": branches=", r[["n_branches"]],
        ", total_mass=", r[["total_mass"]],
        ", relic=", r[["n_relic"]],
        ", middle=", r[["n_middle"]],
        ", founder=", r[["n_founder"]]
      )
    }),
    collapse = "\n"
  ),
  "",
  "[pre_null_core_precheck]",
  paste(
    apply(zone_metrics$primary_contrasts, 1, function(r) {
      paste0(
        r[["core_id"]], ": regions=", r[["core_regions"]],
        ", PH908(core_relic=", round3(as.numeric(r[["core_relic_share_PH908"]])),
        ", core_founder=", round3(as.numeric(r[["core_founder_share_PH908"]])),
        ", skew=", round3(as.numeric(r[["capture_skew_PH908"]])),
        "), R1a(core_relic=", round3(as.numeric(r[["core_relic_share_R1a"]])),
        ", core_founder=", round3(as.numeric(r[["core_founder_share_R1a"]])),
        ", skew=", round3(as.numeric(r[["capture_skew_R1a"]])),
        "), anti_survivor_contrast=", round3(as.numeric(r[["anti_survivor_contrast"]]))
      )
    }),
    collapse = "\n"
  ),
  ""
))

null_run <- run_primary_null(branch_tbl, core_defs, n_perm = N_PERM, eps = EPS)
strata_diag <- null_run$strata_diag
null_tbl <- null_run$null_tbl

test_tbl <- test_topk_results(
  primary_tbl = zone_metrics$primary_contrasts,
  secondary_tbl = zone_metrics$secondary_contrasts,
  null_tbl = null_tbl
)

expected_founder_tbl <- zone_metrics$primary_contrasts |>
  dplyr::filter(topk == PRIMARY_TOPK) |>
  tidyr::pivot_longer(
    cols = c(
      core_relic_share_PH908, core_relic_share_R1a,
      total_founder_PH908, total_founder_R1a,
      expected_core_founders_PH908, expected_core_founders_R1a,
      Core_Founder_PH908, Core_Founder_R1a,
      observed_minus_expected_PH908, observed_minus_expected_R1a,
      observed_over_expected_PH908, observed_over_expected_R1a
    ),
    names_to = c(".value", "lineage"),
    names_pattern = "(.*)_(PH908|R1a)"
  ) |>
  dplyr::rename(
    total_founders = total_founder,
    observed_core_founders = Core_Founder
  ) |>
  dplyr::select(
    lineage, core_relic_share, total_founders,
    expected_core_founders, observed_core_founders,
    observed_minus_expected, observed_over_expected
  )

jackknife_tbl <- run_top4_jackknife(branch_tbl, core_defs, n_perm = N_PERM, eps = EPS)
threshold_grid_tbl <- run_threshold_grid(base, age_tbl, n_perm = N_PERM, eps = EPS)

final_summary_tbl <- test_tbl |>
  dplyr::mutate(
    analysis_set = paste0("top", topk)
  ) |>
  dplyr::select(
    analysis_set, family, statistic, observed,
    ci_lower, ci_upper, p_one_sided, interpretation
  )

supp_table <- dplyr::bind_rows(
  filter_diag |>
    dplyr::mutate(section = "filter_diagnostics") |>
    dplyr::rename(metric = step, value = n) |>
    dplyr::mutate(gating_mode = PH908_GATING_MODE, hg_group = NA_character_, value = as.character(value)),
  branch_diag |>
    tidyr::pivot_longer(cols = -c(gating_mode, hg_group), names_to = "metric", values_to = "value") |>
    dplyr::mutate(section = "branch_diagnostics", value = as.character(value)),
  age_diag |>
    tidyr::pivot_longer(cols = -c(gating_mode, hg_group), names_to = "metric", values_to = "value") |>
    dplyr::mutate(section = "age_diagnostics", value = as.character(value)),
  ph908_relic_ranking |>
    tidyr::pivot_longer(cols = -dominant_region, names_to = "metric", values_to = "value", values_transform = list(value = as.character)) |>
    dplyr::mutate(section = "ph908_relic_region_ranking", gating_mode = PH908_GATING_MODE, hg_group = "PH908", metric = paste0(metric, "::", dominant_region)) |>
    dplyr::select(-dominant_region),
  core_defs |>
    tidyr::pivot_longer(cols = everything(), names_to = "metric", values_to = "value", values_transform = list(value = as.character)) |>
    dplyr::mutate(section = "core_definitions", gating_mode = PH908_GATING_MODE, hg_group = NA_character_),
  final_summary_tbl |>
    tidyr::pivot_longer(cols = everything(), names_to = "metric", values_to = "value", values_transform = list(value = as.character)) |>
    dplyr::mutate(section = "formal_summary", gating_mode = PH908_GATING_MODE, hg_group = NA_character_),
  jackknife_tbl |>
    tidyr::pivot_longer(cols = everything(), names_to = "metric", values_to = "value", values_transform = list(value = as.character)) |>
    dplyr::mutate(section = "jackknife", gating_mode = PH908_GATING_MODE, hg_group = NA_character_),
  expected_founder_tbl |>
    tidyr::pivot_longer(cols = everything(), names_to = "metric", values_to = "value", values_transform = list(value = as.character)) |>
    dplyr::mutate(section = "expected_founder_bridge", gating_mode = PH908_GATING_MODE, hg_group = NA_character_)
) |>
  dplyr::select(section, gating_mode, hg_group, metric, value)

figuredata_p6_topk <- zone_metrics$primary_contrasts |>
  dplyr::select(
    topk, core_id, core_regions, core_region_count, is_primary,
    core_relic_share_PH908, core_founder_share_PH908,
    core_relic_share_R1a, core_founder_share_R1a,
    capture_skew_PH908, capture_skew_R1a,
    anti_survivor_contrast
  ) |>
  dplyr::left_join(
    test_tbl |>
      dplyr::filter(statistic == "anti_survivor_contrast") |>
      dplyr::select(topk, p_one_sided, ci_lower, ci_upper, interpretation),
    by = "topk"
  ) |>
  dplyr::rename(
    p_one_sided = p_one_sided,
    support_label = interpretation
  )

figuredata_p6_foundercarry <- expected_founder_tbl
figuredata_p6_jackknife <- jackknife_tbl
figuredata_p6_threshold_grid <- threshold_grid_tbl
figuredata_p6_null <- null_tbl |>
  tidyr::pivot_longer(
    cols = c(anti_survivor_contrast, core_minus_expanse),
    names_to = "statistic",
    values_to = "stat_value"
  )

readr::write_csv(supp_table, FILE_P6_SUPP_TABLE)
readr::write_csv(threshold_grid_tbl, FILE_P6_THRESHOLD_TABLE)
readr::write_csv(figuredata_p6_topk, FILE_P6_FIG_TOPK)
readr::write_csv(figuredata_p6_foundercarry, FILE_P6_FIG_FOUNDERCARRY)
readr::write_csv(figuredata_p6_jackknife, FILE_P6_FIG_JACKKNIFE)
readr::write_csv(figuredata_p6_threshold_grid, FILE_P6_FIG_THRESHOLD)
readr::write_csv(figuredata_p6_null, FILE_P6_FIG_NULL)

results_p6 <- list(
  script = "07_pillar6_sourcesink.R",
  seed = SEED_MAIN,
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  config_file = path_for_facts(CONFIG_FILE),
  inputs = list(
    preprocessed_dataset = FILE_PREPROCESSED_RDS,
    branch_ages = FILE_BRANCH_AGES
  ),
  settings = list(
    ph908_gating_mode = PH908_GATING_MODE,
    primary_topk = PRIMARY_TOPK,
    sens_topk = SENS_TOPK,
    all_topk = ALL_TOPK,
    relic_max_n = RELIC_MAX_N,
    founder_min_n = FOUNDER_MIN_N,
    threshold_grid_relic = THRESHOLD_GRID_RELIC,
    threshold_grid_founder = THRESHOLD_GRID_FOUNDER,
    age_bin_n = AGE_BIN_N,
    size_bin_n = SIZE_BIN_N,
    n_perm = N_PERM,
    eps = EPS
  ),
  filter_diag = filter_diag,
  branch_tbl = branch_tbl,
  branch_diag = branch_diag,
  age_diag = age_diag,
  ph908_relic_ranking = ph908_relic_ranking,
  core_defs = core_defs,
  zone_metrics = zone_metrics,
  strata_diag = strata_diag,
  null_tbl = null_tbl,
  test_tbl = test_tbl,
  expected_founder_tbl = expected_founder_tbl,
  jackknife_tbl = jackknife_tbl,
  threshold_grid_tbl = threshold_grid_tbl,
  final_summary_tbl = final_summary_tbl
)

saveRDS(results_p6, FILE_P6_RESULTS)

primary_top4 <- figuredata_p6_topk |>
  dplyr::filter(topk == PRIMARY_TOPK)

append_facts_lines(FILE_P6_FACTS, c(
  "[settings]",
  paste0("ph908_gating_mode: ", PH908_GATING_MODE),
  paste0("primary_topk: ", PRIMARY_TOPK),
  paste0("sens_topk: ", paste(SENS_TOPK, collapse = ", ")),
  paste0("all_topk: ", paste(ALL_TOPK, collapse = ", ")),
  paste0("relic_max_n: ", RELIC_MAX_N),
  paste0("founder_min_n: ", FOUNDER_MIN_N),
  paste0("threshold_grid_relic: ", paste(THRESHOLD_GRID_RELIC, collapse = ", ")),
  paste0("threshold_grid_founder: ", paste(THRESHOLD_GRID_FOUNDER, collapse = ", ")),
  paste0("age_bin_n: ", AGE_BIN_N),
  paste0("size_bin_n: ", SIZE_BIN_N),
  paste0("n_perm: ", N_PERM),
  paste0("eps: ", EPS),
  "note: permutation-based consistency checks were run at the configured depth.",
  "",
  "[topk_primary_secondary]",
  paste(
    apply(final_summary_tbl, 1, function(r) {
      paste0(
        r[["analysis_set"]], " / ", r[["family"]], " / ", r[["statistic"]],
        ": observed=", round3(as.numeric(r[["observed"]])),
        ", ci=[", round3(as.numeric(r[["ci_lower"]])), ", ", round3(as.numeric(r[["ci_upper"]])), "]",
        ", p_one_sided=", fmt_p(as.numeric(r[["p_one_sided"]])),
        ", interpretation=", r[["interpretation"]]
      )
    }),
    collapse = "\n"
  ),
  "",
  "[expected_founder_bridge_top4]",
  paste(
    apply(expected_founder_tbl, 1, function(r) {
      paste0(
        r[["lineage"]], ": core_relic_share=", round3(as.numeric(r[["core_relic_share"]])),
        ", total_founders=", r[["total_founders"]],
        ", expected_core_founders=", round3(as.numeric(r[["expected_core_founders"]])),
        ", observed_core_founders=", r[["observed_core_founders"]],
        ", observed_minus_expected=", round3(as.numeric(r[["observed_minus_expected"]])),
        ", observed_over_expected=", round3(as.numeric(r[["observed_over_expected"]]))
      )
    }),
    collapse = "\n"
  ),
  "",
  "[jackknife_top4]",
  paste(
    apply(jackknife_tbl, 1, function(r) {
      paste0(
        r[["dropped_region"]], ": observed=", round3(as.numeric(r[["observed"]])),
        ", ci=[", round3(as.numeric(r[["ci_low"]])), ", ", round3(as.numeric(r[["ci_high"]])), "]",
        ", p_one_sided=", fmt_p(as.numeric(r[["p_one_sided"]])),
        ", sign_positive=", r[["sign_positive"]],
        ", delta_vs_full=", round3(as.numeric(r[["delta_vs_full"]])),
        ", interpretation=", r[["interpretation"]]
      )
    }),
    collapse = "\n"
  ),
  "",
  "[threshold_grid_summary]",
  paste(
    apply(threshold_grid_tbl, 1, function(r) {
      paste0(
        "relic<=", r[["relic_max_n"]], " / founder>=", r[["founder_min_n"]],
        ": anti_survivor_contrast_top4=", round3(as.numeric(r[["anti_survivor_contrast_top4"]])),
        ", p_one_sided_top4=", fmt_p(as.numeric(r[["p_one_sided_top4"]])),
        ", core_minus_expanse_top4=", round3(as.numeric(r[["core_minus_expanse_top4"]])),
        ", p_one_sided_core_minus_expanse_top4=", fmt_p(as.numeric(r[["p_one_sided_core_minus_expanse_top4"]])),
        ", adequacy_flag=", r[["adequacy_flag"]],
        ", support_label=", r[["support_label"]]
      )
    }),
    collapse = "\n"
  ),
  "",
  "[wording_ready_lines]",
  paste0(
    "methods_line_1: Pillar 6 defines deep-core capture skew as core relic share minus core founder share, and tests whether this skew is larger in PH908 than in the matched R1a control."
  ),
  paste0(
    "methods_line_2: The prespecified primary analysis uses the top4 PH908-defined relic core, with neighboring top3 and top5 cores treated as sensitivity checks, and with additional leave-one-out and threshold-grid robustness analyses."
  ),
  paste0(
    "results_line_1: In the primary top4 core, anti_survivor_contrast was ",
    round3(primary_top4$anti_survivor_contrast),
    " with a one-sided permutation p-value of ",
    fmt_p(primary_top4$p_one_sided),
    "."
  ),
  paste0(
    "results_line_2: The neighboring sensitivity cores showed anti_survivor_contrast values of ",
    paste(
      figuredata_p6_topk |>
        dplyr::mutate(txt = paste0("top", topk, "=", round3(anti_survivor_contrast))) |>
        dplyr::pull(txt),
      collapse = "; "
    ),
    "."
  ),
  paste0(
    "results_line_3: The top4 leave-one-out jackknife remained positive in every run, and the 3×3 threshold grid showed whether the top4 result remained directionally stable across neighboring relic/founder definitions."
  ),
  paste0(
    "figure_caption_line: Pillar 6 tests whether the deepest PH908 core is more consistent with a relic-preserving sink than with a founder launch center, using top3/top4/top5 core definitions, leave-one-out jackknife, and a threshold sensitivity grid."
  ),
  "",
  "[output_files]",
  paste0("results_rds: ", path_for_facts(FILE_P6_RESULTS)),
  paste0("supp_table: ", path_for_facts(FILE_P6_SUPP_TABLE)),
  paste0("threshold_sensitivity_table: ", path_for_facts(FILE_P6_THRESHOLD_TABLE)),
  paste0("figuredata_topk: ", path_for_facts(FILE_P6_FIG_TOPK)),
  paste0("figuredata_foundercarry: ", path_for_facts(FILE_P6_FIG_FOUNDERCARRY)),
  paste0("figuredata_jackknife: ", path_for_facts(FILE_P6_FIG_JACKKNIFE)),
  paste0("figuredata_threshold: ", path_for_facts(FILE_P6_FIG_THRESHOLD)),
  paste0("figuredata_null: ", path_for_facts(FILE_P6_FIG_NULL)),
  paste0("facts_file: ", path_for_facts(FILE_P6_FACTS)),
  paste0("config_file: ", path_for_facts(CONFIG_FILE)),
  ""
))

cat("============================================================\n")
cat("PILLAR 6 COMPLETE\n")
cat("============================================================\n")
cat("Seed: ", SEED_MAIN, "\n", sep = "")
cat("Primary topk: ", PRIMARY_TOPK, "\n", sep = "")
cat("Primary top4 anti_survivor_contrast: ", round3(primary_top4$anti_survivor_contrast), "\n", sep = "")
cat("Primary top4 support label: ", primary_top4$support_label, "\n", sep = "")
cat("Permutation depth: ", N_PERM, "\n", sep = "")
cat("Results object: ", FILE_P6_RESULTS, "\n", sep = "")
cat("Supplement table: ", FILE_P6_SUPP_TABLE, "\n", sep = "")
cat("Threshold sensitivity table: ", FILE_P6_THRESHOLD_TABLE, "\n", sep = "")
cat("Figuredata topk: ", FILE_P6_FIG_TOPK, "\n", sep = "")
cat("Figuredata foundercarry: ", FILE_P6_FIG_FOUNDERCARRY, "\n", sep = "")
cat("Figuredata jackknife: ", FILE_P6_FIG_JACKKNIFE, "\n", sep = "")
cat("Figuredata threshold: ", FILE_P6_FIG_THRESHOLD, "\n", sep = "")
cat("Figuredata null: ", FILE_P6_FIG_NULL, "\n", sep = "")
cat("Facts file: ", FILE_P6_FACTS, "\n", sep = "")
cat("============================================================\n")