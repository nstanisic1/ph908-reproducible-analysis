# 02_pillar1_burden.R
# Phylogenetic burden-capacity bound for PH908 versus R1a.

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

FILE_P1_RESULTS <- file.path(DIR_RESULTS, "results_p1.rds")
FILE_P1_FIG_MAIN <- file.path(DIR_FIGUREDATA, "figuredata_p1_main.csv")
FILE_P1_FIG_S2 <- file.path(DIR_FIGUREDATA, "figuredata_p1_s2.csv")
FILE_P1_FIG_S3 <- file.path(DIR_FIGUREDATA, "figuredata_p1_s3.csv")
FILE_P1_SUPP_TABLE <- file.path(DIR_TABLES, "table_p1_supp.csv")
FILE_P1_RX_TABLE <- file.path(DIR_TABLES, "table_p1_rx.csv")
FILE_P1_FACTS <- file.path(DIR_FACTS, "facts_p1.txt")

RELIC_MAX_N <- 3L
FOUNDER_MIN_N <- 15L
YOUNG_BRANCH_CUTOFF_YBP <- 2000
SURVIVAL_FLOORS <- c(1L, 2L, 3L)
GATING_MODE_MAIN <- "primary"

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
cat("02 PILLAR 1 — PHYLOGENETIC BURDEN-CAPACITY BOUND\n")
cat("============================================================\n")
cat("Root: ", DIR_ROOT, "\n", sep = "")
cat("Seed: ", SEED_MAIN, "\n", sep = "")
cat("Inputs:\n")
cat("  - ", FILE_PREPROCESSED_RDS, "\n", sep = "")
cat("  - ", FILE_BRANCH_AGES, "\n\n", sep = "")

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

writeLines(
  c(
    "script: 02_pillar1_burden.R",
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
  FILE_P1_FACTS,
  useBytes = TRUE
)

normalize_snp <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("^I-|^R-", "", x)
  x <- gsub("\\*$", "", x)
  x
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_prop <- function(num, den) {
  ifelse(den > 0, num / den, NA_real_)
}

top_share <- function(x, k) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  sum(sort(x, decreasing = TRUE)[seq_len(min(k, length(x)))]) / sum(x)
}

classify_branch <- function(n) {
  dplyr::case_when(
    n <= RELIC_MAX_N ~ "Relic",
    n >= FOUNDER_MIN_N ~ "Founder",
    TRUE ~ "Middle"
  )
}

choose_candidate <- function(branch_size, class_now, tmrca_final, model) {
  dplyr::case_when(
    model == "any_branch" ~ TRUE,
    model == "non_relic_only" ~ branch_size > RELIC_MAX_N,
    model == "non_singleton" ~ branch_size >= 2,
    model == "young_non_relic" ~ branch_size > RELIC_MAX_N &
      is.finite(tmrca_final) &
      tmrca_final <= YOUNG_BRANCH_CUTOFF_YBP,
    TRUE ~ FALSE
  )
}

round1 <- function(x) round(x, 1)

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

assert_one_dec <- function(actual, expected, msg, tol = 0.35) {
  if (length(actual) != 1 || is.na(actual)) {
    stop(sprintf("%s Expected %.1f, got NA/invalid.", msg, expected), call. = FALSE)
  }
  if (abs(actual - expected) > tol) {
    stop(sprintf("%s Expected %.1f, got %.1f.", msg, expected, round(actual, 1)), call. = FALSE)
  }
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

build_branch_table <- function(cl_dat, age_tbl, gating_mode = "primary") {
  base <- cl_dat |>
    dplyr::filter(
      matched == TRUE,
      exclude_geo_primary == FALSE,
      is_balkan_country == TRUE
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
    dplyr::filter(!is.na(hg_group))
  
  if (nrow(base) == 0) {
    stop("No rows remain after Pillar 1 filtering.", call. = FALSE)
  }
  
  branch_tbl <- base |>
    dplyr::group_by(hg_group, terminal_snp_norm) |>
    dplyr::summarise(
      branch_size = dplyr::n(),
      dataset_tmrca = safe_mean(tmrca),
      dominant_region = names(sort(table(Region), decreasing = TRUE))[1],
      .groups = "drop"
    ) |>
    dplyr::left_join(age_tbl, by = "terminal_snp_norm") |>
    dplyr::mutate(
      tmrca_final = dplyr::coalesce(tmrca_ybp, dataset_tmrca),
      age_source_final = dplyr::case_when(
        !is.na(tmrca_ybp) ~ "external_age",
        is.na(tmrca_ybp) & !is.na(dataset_tmrca) ~ "dataset_age",
        TRUE ~ "missing"
      ),
      class_now = classify_branch(branch_size),
      gating_mode = gating_mode
    ) |>
    dplyr::arrange(hg_group, dplyr::desc(branch_size), terminal_snp_norm)
  
  required_hgs <- c("PH908", "R1a")
  missing_hgs <- setdiff(required_hgs, unique(branch_tbl$hg_group))
  if (length(missing_hgs) > 0) {
    stop(
      paste0("Required haplogroup(s) absent after filtering: ", paste(missing_hgs, collapse = ", ")),
      call. = FALSE
    )
  }
  
  branch_tbl
}

make_branch_diagnostics <- function(branch_tbl) {
  branch_tbl |>
    dplyr::group_by(gating_mode, hg_group) |>
    dplyr::summarise(
      n_branches = dplyr::n(),
      total_branch_mass = sum(branch_size),
      mean_branch_size = mean(branch_size),
      median_branch_size = stats::median(branch_size),
      min_branch_size = min(branch_size),
      max_branch_size = max(branch_size),
      n_relic = sum(class_now == "Relic"),
      n_middle = sum(class_now == "Middle"),
      n_founder = sum(class_now == "Founder"),
      prop_relic = mean(class_now == "Relic"),
      prop_middle = mean(class_now == "Middle"),
      prop_founder = mean(class_now == "Founder"),
      .groups = "drop"
    )
}

make_age_coverage_diagnostics <- function(branch_tbl) {
  branch_tbl |>
    dplyr::group_by(gating_mode, hg_group) |>
    dplyr::summarise(
      n_branches_total = dplyr::n(),
      total_branch_mass = sum(branch_size),
      n_external_age = sum(age_source_final == "external_age"),
      n_dataset_age = sum(age_source_final == "dataset_age"),
      n_missing_age = sum(age_source_final == "missing"),
      mass_external_age = sum(branch_size[age_source_final == "external_age"]),
      mass_dataset_age = sum(branch_size[age_source_final == "dataset_age"]),
      mass_missing_age = sum(branch_size[age_source_final == "missing"]),
      prop_branches_missing_age = mean(age_source_final == "missing"),
      prop_mass_missing_age = safe_prop(sum(branch_size[age_source_final == "missing"]), sum(branch_size)),
      n_branches_age_eligible_for_young_model = sum(is.finite(tmrca_final)),
      mass_age_eligible_for_young_model = sum(branch_size[is.finite(tmrca_final)]),
      prop_branches_age_eligible = mean(is.finite(tmrca_final)),
      prop_mass_age_eligible = safe_prop(sum(branch_size[is.finite(tmrca_final)]), sum(branch_size)),
      .groups = "drop"
    )
}

make_branch_concentration_audit <- function(branch_tbl) {
  branch_tbl |>
    dplyr::group_by(gating_mode, hg_group) |>
    dplyr::summarise(
      n_branches = dplyr::n(),
      total_branch_mass = sum(branch_size),
      largest_branch = max(branch_size),
      second_largest_branch = sort(branch_size, decreasing = TRUE)[min(2, length(branch_size))],
      fifth_largest_branch = sort(branch_size, decreasing = TRUE)[min(5, length(branch_size))],
      share_top1 = top_share(branch_size, 1),
      share_top2 = top_share(branch_size, 2),
      share_top5 = top_share(branch_size, 5),
      .groups = "drop"
    )
}

run_capacity <- function(branch_tbl, scenario_tbl, carrier_models) {
  results <- list()
  k <- 0
  
  for (ms in SURVIVAL_FLOORS) {
    for (hg in sort(unique(branch_tbl$hg_group))) {
      df_hg <- branch_tbl |>
        dplyr::filter(hg_group == hg)
      total_mass <- sum(df_hg$branch_size)
      
      scen <- scenario_tbl |>
        dplyr::filter(hg_group == hg)
      
      for (s in seq_len(nrow(scen))) {
        amp <- scen$amplification_factor[s]
        target_start_mass <- total_mass / amp
        required_reduction <- total_mass - target_start_mass
        
        for (model in carrier_models) {
          candidate_flag <- choose_candidate(
            df_hg$branch_size,
            df_hg$class_now,
            df_hg$tmrca_final,
            model
          )
          
          reducible_capacity <- sum(pmax(0, df_hg$branch_size[candidate_flag] - ms))
          capacity_margin <- reducible_capacity - required_reduction
          
          k <- k + 1
          results[[k]] <- tibble::tibble(
            gating_mode = unique(df_hg$gating_mode),
            hg_group = hg,
            scenario_family = scen$scenario_family[s],
            scenario_id = scen$scenario_id[s],
            carrier_model = model,
            min_start_size = ms,
            total_mass = total_mass,
            amplification_factor = amp,
            target_start_mass = target_start_mass,
            required_reduction = required_reduction,
            n_candidate_branches = sum(candidate_flag),
            max_reduction_capacity = reducible_capacity,
            capacity_margin = capacity_margin,
            feasible = capacity_margin >= 0
          )
        }
      }
    }
  }
  
  dplyr::bind_rows(results)
}

make_supp_table <- function(branch_diag, age_diag, branch_audit) {
  dplyr::bind_rows(
    branch_diag |>
      dplyr::mutate(section = "branch_diagnostics") |>
      tidyr::pivot_longer(cols = -c(section, gating_mode, hg_group), names_to = "metric", values_to = "value"),
    age_diag |>
      dplyr::mutate(section = "age_coverage") |>
      tidyr::pivot_longer(cols = -c(section, gating_mode, hg_group), names_to = "metric", values_to = "value"),
    branch_audit |>
      dplyr::mutate(section = "branch_concentration") |>
      tidyr::pivot_longer(cols = -c(section, gating_mode, hg_group), names_to = "metric", values_to = "value")
  ) |>
    dplyr::arrange(section, gating_mode, hg_group, metric)
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
  
  stop(sprintf("Unknown perturbation: %s", perturbation), call. = FALSE)
}

summarise_main <- function(main_results) {
  main_results |>
    dplyr::group_by(gating_mode, hg_group, min_start_size) |>
    dplyr::summarise(
      worst_margin = min(capacity_margin),
      best_margin = max(capacity_margin),
      median_margin = median(capacity_margin),
      all_carrier_models_feasible = all(feasible),
      any_carrier_model_feasible = any(feasible),
      n_feasible_models = sum(feasible),
      n_total_models = dplyr::n(),
      .groups = "drop"
    )
}

summarise_grid <- function(grid_results) {
  grid_results |>
    dplyr::group_by(gating_mode, hg_group, min_start_size) |>
    dplyr::summarise(
      worst_margin = min(capacity_margin),
      best_margin = max(capacity_margin),
      median_margin = median(capacity_margin),
      prop_grid_feasible = mean(feasible),
      all_grid_points_feasible = all(feasible),
      n_grid_points = dplyr::n(),
      n_feasible_grid_points = sum(feasible),
      .groups = "drop"
    )
}

carrier_models <- c(
  "any_branch",
  "non_relic_only",
  "non_singleton",
  "young_non_relic"
)

cl_raw <- readRDS(FILE_PREPROCESSED_RDS)

if (!"is_balkan_country" %in% names(cl_raw)) {
  if (!"country_code" %in% names(cl_raw)) {
    stop(
      "Dataset is missing both is_balkan_country and country_code, so the Pillar 1 Balkan-country filter cannot be reconstructed.",
      call. = FALSE
    )
  }
  cl_raw <- cl_raw |>
    dplyr::mutate(
      is_balkan_country = toupper(country_code) %in% BALKAN_ISO2
    )
  message("Reconstructed is_balkan_country from country_code for Pillar 1 filtering.")
}

required_cols <- c(
  "matched", "exclude_geo_primary", "is_balkan_country",
  "terminal_snp", "major_hg",
  "is_ph908_primary", "is_ph908_accurate", "is_ph908_wide",
  "tmrca", "Region"
)
missing_cols <- setdiff(required_cols, names(cl_raw))
if (length(missing_cols) > 0) {
  stop(
    paste0("Dataset missing required columns after Pillar 1 input preparation: ", paste(missing_cols, collapse = ", ")),
    call. = FALSE
  )
}

age_tbl <- standardize_age_table(readr::read_csv(FILE_BRANCH_AGES, show_col_types = FALSE))

main_scenario <- tibble::tribble(
  ~scenario_family,   ~scenario_id,        ~hg_group, ~modern_freq, ~Pi, ~M,  ~Ps,
  "strict_exogenous", "strict_exogenous", "PH908",     0.65,         0,   0.5, 0.2,
  "strict_exogenous", "strict_exogenous", "R1a",       0.15,         0,   0.5, 0.2
) |>
  dplyr::mutate(
    immediate_post_admix = ((1 - M) * Pi) + (M * Ps),
    amplification_factor = modern_freq / immediate_post_admix
  )

grid_scenarios <- tidyr::crossing(
  M = c(0.45, 0.50, 0.55),
  Ps = c(0.18, 0.20, 0.22, 0.25),
  hg_group = c("PH908", "R1a")
) |>
  dplyr::mutate(
    scenario_family = "strict_exogenous_grid",
    Pi = 0,
    modern_freq = ifelse(hg_group == "PH908", 0.65, 0.15),
    immediate_post_admix = M * Ps,
    amplification_factor = modern_freq / immediate_post_admix,
    scenario_id = paste0("grid_M", format(M, nsmall = 2), "_Ps", format(Ps, nsmall = 2))
  ) |>
  dplyr::select(
    scenario_family, scenario_id, hg_group, modern_freq, Pi, M, Ps,
    immediate_post_admix, amplification_factor
  )

branch_tbl_main <- build_branch_table(cl_raw, age_tbl, gating_mode = GATING_MODE_MAIN)
branch_diag <- make_branch_diagnostics(branch_tbl_main)
age_diag <- make_age_coverage_diagnostics(branch_tbl_main)
branch_audit <- make_branch_concentration_audit(branch_tbl_main)

main_results <- run_capacity(branch_tbl_main, main_scenario, carrier_models) |>
  dplyr::arrange(hg_group, carrier_model, min_start_size)

main_summary <- summarise_main(main_results)

grid_results <- run_capacity(branch_tbl_main, grid_scenarios, carrier_models) |>
  dplyr::left_join(
    grid_scenarios |>
      dplyr::select(scenario_id, hg_group, M, Ps, Pi, modern_freq),
    by = c("scenario_id", "hg_group")
  ) |>
  dplyr::arrange(hg_group, scenario_id, carrier_model, min_start_size)

grid_summary <- summarise_grid(grid_results)

manuscript_tbl <- main_summary |>
  dplyr::left_join(
    grid_summary |>
      dplyr::select(
        gating_mode, hg_group, min_start_size,
        grid_worst_margin = worst_margin,
        grid_best_margin = best_margin,
        prop_grid_feasible
      ),
    by = c("gating_mode", "hg_group", "min_start_size")
  ) |>
  dplyr::left_join(
    branch_diag |>
      dplyr::select(
        gating_mode, hg_group,
        n_branches, total_branch_mass,
        n_relic, n_middle, n_founder,
        prop_relic, prop_middle, prop_founder
      ),
    by = c("gating_mode", "hg_group")
  ) |>
  dplyr::arrange(hg_group, min_start_size)

# Consistency checks
ph_diag <- branch_diag |> dplyr::filter(hg_group == "PH908")
r1a_diag <- branch_diag |> dplyr::filter(hg_group == "R1a")

assert_equal(ph_diag$n_branches, 79L, "PH908 branch count drift detected.")
assert_equal(ph_diag$total_branch_mass, 623L, "PH908 total branch mass drift detected.")
assert_equal(ph_diag$n_relic, 50L, "PH908 relic count drift detected.")
assert_equal(ph_diag$n_middle, 20L, "PH908 middle count drift detected.")
assert_equal(ph_diag$n_founder, 9L, "PH908 founder count drift detected.")

assert_equal(r1a_diag$n_branches, 82L, "R1a branch count drift detected.")
assert_equal(r1a_diag$total_branch_mass, 730L, "R1a total branch mass drift detected.")
assert_equal(r1a_diag$n_relic, 39L, "R1a relic count drift detected.")
assert_equal(r1a_diag$n_middle, 29L, "R1a middle count drift detected.")
assert_equal(r1a_diag$n_founder, 14L, "R1a founder count drift detected.")

get_main_row <- function(hg, ms) {
  main_summary |> dplyr::filter(hg_group == hg, min_start_size == ms)
}

get_grid_row <- function(hg, ms) {
  grid_summary |> dplyr::filter(hg_group == hg, min_start_size == ms)
}

for (spec in list(
  list(hg = "PH908", ms = 1L, worst = -15.2, best = 16.8, feasible = 2L),
  list(hg = "PH908", ms = 2L, worst = -44.2, best = -32.2, feasible = 0L),
  list(hg = "PH908", ms = 3L, worst = -73.2, best = -73.2, feasible = 0L),
  list(hg = "R1a", ms = 1L, worst = 117.0, best = 405.0, feasible = 4L),
  list(hg = "R1a", ms = 2L, worst = 85.7, best = 347.0, feasible = 4L),
  list(hg = "R1a", ms = 3L, worst = 54.7, best = 301.0, feasible = 4L)
)) {
  row <- get_main_row(spec$hg, spec$ms)
  assert_one_dec(row$worst_margin, spec$worst, sprintf("Main worst margin drift for %s floor %s.", spec$hg, spec$ms))
  assert_one_dec(row$best_margin, spec$best, sprintf("Main best margin drift for %s floor %s.", spec$hg, spec$ms))
  assert_equal(row$n_feasible_models, spec$feasible, sprintf("Main feasible-model count drift for %s floor %s.", spec$hg, spec$ms))
}

for (spec in list(
  list(hg = "PH908", ms = 1L, worst = -33.4, best = 52.8, feasible = 0.5833333),
  list(hg = "PH908", ms = 2L, worst = -62.4, best = 3.8, feasible = 0.0416667),
  list(hg = "PH908", ms = 3L, worst = -91.4, best = -37.2, feasible = 0.0),
  list(hg = "R1a", ms = 1L, worst = 24.2, best = 587.0, feasible = 1.0),
  list(hg = "R1a", ms = 2L, worst = -6.8, best = 529.0, feasible = 0.9791667),
  list(hg = "R1a", ms = 3L, worst = -37.8, best = 483.0, feasible = 0.9791667)
)) {
  row <- get_grid_row(spec$hg, spec$ms)
  assert_one_dec(row$worst_margin, spec$worst, sprintf("Grid worst margin drift for %s floor %s.", spec$hg, spec$ms))
  assert_one_dec(row$best_margin, spec$best, sprintf("Grid best margin drift for %s floor %s.", spec$hg, spec$ms))
  assert_one_dec(100 * row$prop_grid_feasible, 100 * spec$feasible, sprintf("Grid feasible proportion drift for %s floor %s.", spec$hg, spec$ms))
}

# Figuredata
figuredata_p1_main <- manuscript_tbl |>
  dplyr::transmute(
    gating_mode,
    hg_group,
    min_start_size,
    n_branches,
    total_branch_mass,
    n_relic,
    n_middle,
    n_founder,
    prop_relic,
    prop_middle,
    prop_founder,
    worst_margin,
    best_margin,
    median_margin,
    n_feasible_models,
    n_total_models,
    grid_worst_margin,
    grid_best_margin,
    prop_grid_feasible
  )

figuredata_p1_s2 <- grid_results |>
  dplyr::transmute(
    gating_mode,
    hg_group,
    scenario_id,
    carrier_model,
    min_start_size,
    M,
    Ps,
    Pi,
    modern_freq,
    amplification_factor,
    capacity_margin,
    feasible
  )

figuredata_p1_s3 <- branch_tbl_main |>
  dplyr::transmute(
    gating_mode,
    hg_group,
    terminal_snp_norm,
    branch_size,
    class_now,
    tmrca_final,
    age_source_final,
    dominant_region
  )

table_p1_supp <- make_supp_table(branch_diag, age_diag, branch_audit)

# Robustness
robustness_specs <- tibble::tribble(
  ~perturbation,
  "baseline",
  "drop_largest_PH908",
  "drop_largest_R1a",
  "drop_top2_cumulative"
)

robustness_results <- purrr::map_dfr(robustness_specs$perturbation, function(perturbation) {
  branch_tbl_rx <- make_robustness_branch_table(branch_tbl_main, perturbation)
  main_results_rx <- run_capacity(branch_tbl_rx, main_scenario, carrier_models)
  main_summary_rx <- summarise_main(main_results_rx)
  rx_counts <- branch_tbl_rx |>
    dplyr::group_by(hg_group) |>
    dplyr::summarise(
      n_branches_after = dplyr::n(),
      total_mass_after = sum(branch_size),
      .groups = "drop"
    )
  
  main_summary_rx |>
    dplyr::left_join(rx_counts, by = "hg_group") |>
    dplyr::mutate(perturbation = perturbation)
}) |>
  dplyr::relocate(perturbation, .before = gating_mode) |>
  dplyr::arrange(perturbation, hg_group, min_start_size)

baseline_sign <- main_summary |>
  dplyr::select(hg_group, min_start_size, baseline_worst_margin = worst_margin)

robustness_results <- robustness_results |>
  dplyr::left_join(baseline_sign, by = c("hg_group", "min_start_size")) |>
  dplyr::mutate(
    sign_flip_vs_baseline = sign(worst_margin) != sign(baseline_worst_margin)
  )

readr::write_csv(figuredata_p1_main, FILE_P1_FIG_MAIN)
readr::write_csv(figuredata_p1_s2, FILE_P1_FIG_S2)
readr::write_csv(figuredata_p1_s3, FILE_P1_FIG_S3)
readr::write_csv(table_p1_supp, FILE_P1_SUPP_TABLE)
readr::write_csv(robustness_results, FILE_P1_RX_TABLE)

# Exports
results_p1 <- list(
  metadata = list(
    script = "02_pillar1_burden.R",
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    seed = SEED_MAIN,
    config_file = path_for_facts(CONFIG_FILE),
    inputs = c(
      path_for_facts(FILE_PREPROCESSED_RDS),
      path_for_facts(FILE_BRANCH_AGES)
    ),
    gating_mode_main = GATING_MODE_MAIN,
    thresholds = list(
      relic_max_n = RELIC_MAX_N,
      founder_min_n = FOUNDER_MIN_N,
      young_branch_cutoff_ybp = YOUNG_BRANCH_CUTOFF_YBP,
      survival_floors = SURVIVAL_FLOORS
    ),
    main_scenario = main_scenario,
    grid_values = list(M = c(0.45, 0.50, 0.55), Ps = c(0.18, 0.20, 0.22, 0.25)),
    git_commit = git_commit
  ),
  branch_tbl_main = branch_tbl_main,
  branch_diagnostics = branch_diag,
  age_coverage_diagnostics = age_diag,
  branch_concentration_audit = branch_audit,
  main_results = main_results,
  main_summary = main_summary,
  grid_results = grid_results,
  grid_summary = grid_summary,
  manuscript_table = manuscript_tbl,
  figuredata_main = figuredata_p1_main,
  figuredata_s2 = figuredata_p1_s2,
  figuredata_s3 = figuredata_p1_s3,
  table_p1_supp = table_p1_supp,
  robustness_results = robustness_results,
  freeze_checks = list(
    ph908_branches = 79L,
    ph908_mass = 623L,
    r1a_branches = 82L,
    r1a_mass = 730L,
    ph908_classes = c(relic = 50L, middle = 20L, founder = 9L),
    r1a_classes = c(relic = 39L, middle = 29L, founder = 14L)
  )
)

saveRDS(results_p1, FILE_P1_RESULTS)

# Facts
sign_flip_any <- any(
  robustness_results$perturbation != "baseline" &
    robustness_results$sign_flip_vs_baseline,
  na.rm = TRUE
)

facts_lines <- c(
  readLines(FILE_P1_FACTS, warn = FALSE),
  "thresholds_used:",
  paste0("  relic_max_n: ", RELIC_MAX_N),
  paste0("  founder_min_n: ", FOUNDER_MIN_N),
  paste0("  young_branch_cutoff_ybp: ", YOUNG_BRANCH_CUTOFF_YBP),
  paste0("  survival_floors: ", paste(SURVIVAL_FLOORS, collapse = ", ")),
  paste0("  gating_mode_main: ", GATING_MODE_MAIN),
  paste0("  main_scenario: Pi=0, M=0.5, Ps=0.2; modern PH908=0.65, R1a=0.15"),
  paste0("  grid_M: ", paste(c(0.45, 0.50, 0.55), collapse = ", ")),
  paste0("  grid_Ps: ", paste(c(0.18, 0.20, 0.22, 0.25), collapse = ", ")),
  "",
  "branch_counts_and_mass:",
  paste0("  PH908_branches: ", ph_diag$n_branches),
  paste0("  PH908_total_mass: ", ph_diag$total_branch_mass),
  paste0("  R1a_branches: ", r1a_diag$n_branches),
  paste0("  R1a_total_mass: ", r1a_diag$total_branch_mass),
  paste0("  PH908_classes: relic=", ph_diag$n_relic, ", middle=", ph_diag$n_middle, ", founder=", ph_diag$n_founder),
  paste0("  R1a_classes: relic=", r1a_diag$n_relic, ", middle=", r1a_diag$n_middle, ", founder=", r1a_diag$n_founder),
  "",
  "main_summary_table1:",
  vapply(seq_len(nrow(manuscript_tbl)), function(i) {
    row <- manuscript_tbl[i, ]
    paste0(
      "  ", row$hg_group[[1]], " floor_", row$min_start_size[[1]],
      ": main_worst=", round1(as.numeric(row$worst_margin[[1]])),
      ", main_best=", round1(as.numeric(row$best_margin[[1]])),
      ", feasible_models=", row$n_feasible_models[[1]], "/", row$n_total_models[[1]],
      ", grid_worst=", round1(as.numeric(row$grid_worst_margin[[1]])),
      ", grid_best=", round1(as.numeric(row$grid_best_margin[[1]])),
      ", grid_feasible_pct=", round1(100 * as.numeric(row$prop_grid_feasible[[1]]))
    )
  }, character(1)),
  "",
  "wording_ready_lines:",
  paste0(
    "  PH908 is represented by ", ph_diag$n_branches, " terminal branches (total mass ", ph_diag$total_branch_mass,
    "), versus ", r1a_diag$n_branches, " branches (total mass ", r1a_diag$total_branch_mass,
    ") for R1a under the primary Balkan-filtered comparison set."
  ),
  paste0(
    "  Under the strict exogenous main scenario (Pi=0, M=0.5, Ps=0.2), PH908 is feasible in ",
    get_main_row("PH908", 1L)$n_feasible_models, "/", get_main_row("PH908", 1L)$n_total_models,
    " carrier models at floor 1, but in ", get_main_row("PH908", 2L)$n_feasible_models, "/",
    get_main_row("PH908", 2L)$n_total_models, " at floor 2 and ", get_main_row("PH908", 3L)$n_feasible_models,
    "/", get_main_row("PH908", 3L)$n_total_models, " at floor 3."
  ),
  paste0(
    "  R1a remains feasible in ", get_main_row("R1a", 1L)$n_feasible_models, "/", get_main_row("R1a", 1L)$n_total_models,
    ", ", get_main_row("R1a", 2L)$n_feasible_models, "/", get_main_row("R1a", 2L)$n_total_models,
    ", and ", get_main_row("R1a", 3L)$n_feasible_models, "/", get_main_row("R1a", 3L)$n_total_models,
    " carrier models across floors 1-3, respectively."
  ),
  "",
  "sensitivity_summary:",
  paste0(
    "  Influential-unit robustness scenarios run: ",
    paste(robustness_specs$perturbation[robustness_specs$perturbation != "baseline"], collapse = ", "),
    "."
  ),
  paste0("  Any sign flip versus baseline worst-margin direction: ", ifelse(sign_flip_any, "YES", "NO")),
  paste0("sign_flip_anywhere: ", ifelse(sign_flip_any, "yes", "no")),
  "",
  "output_files:",
  paste0("  - ", path_for_facts(FILE_P1_RESULTS)),
  paste0("  - ", path_for_facts(FILE_P1_FIG_MAIN)),
  paste0("  - ", path_for_facts(FILE_P1_FIG_S2)),
  paste0("  - ", path_for_facts(FILE_P1_FIG_S3)),
  paste0("  - ", path_for_facts(FILE_P1_SUPP_TABLE)),
  paste0("  - ", path_for_facts(FILE_P1_RX_TABLE))
)

writeLines(facts_lines, FILE_P1_FACTS, useBytes = TRUE)

cat("============================================================\n")
cat("PILLAR 1 COMPLETE\n")
cat("============================================================\n")
cat("Results object: ", FILE_P1_RESULTS, "\n", sep = "")
cat("Facts file: ", FILE_P1_FACTS, "\n", sep = "")
cat("Figuredata main: ", FILE_P1_FIG_MAIN, "\n", sep = "")
cat("Figuredata S2: ", FILE_P1_FIG_S2, "\n", sep = "")
cat("Figuredata S3: ", FILE_P1_FIG_S3, "\n", sep = "")
cat("Supplement table: ", FILE_P1_SUPP_TABLE, "\n", sep = "")
cat("Robustness table: ", FILE_P1_RX_TABLE, "\n", sep = "")
cat("============================================================\n")