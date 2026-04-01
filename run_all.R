# run_all.R
# Master runner for the PH908 reproducible analysis pipeline.

# -------------------------------------------------------------------
# Runner settings
# -------------------------------------------------------------------
CONTINUE_ON_ERROR <- FALSE
STRICT_OUTPUT_CHECKS <- TRUE

# -------------------------------------------------------------------
# Helper: locate config
# -------------------------------------------------------------------
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

# -------------------------------------------------------------------
# Boot
# -------------------------------------------------------------------
CONFIG_FILE <- locate_config_file()
source(CONFIG_FILE)
setwd(DIR_ROOT)

DIR_LOGS <- file.path(DIR_ROOT, "logs")
dir.create(DIR_LOGS, showWarnings = FALSE, recursive = TRUE)

run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
LOG_FILE <- file.path(DIR_LOGS, paste0("pipeline_run_", run_stamp, ".log"))
SUMMARY_CSV <- file.path(DIR_LOGS, paste0("pipeline_summary_", run_stamp, ".csv"))
WARNINGS_TXT <- file.path(DIR_LOGS, paste0("pipeline_warnings_", run_stamp, ".txt"))

# -------------------------------------------------------------------
# Logging
# -------------------------------------------------------------------
log_line <- function(..., sep = "") {
  msg <- paste0(...)
  timestamped <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(timestamped, "\n", sep = "")
  cat(timestamped, "\n", sep = "", file = LOG_FILE, append = TRUE)
}

write_block <- function(title_char = "=", width = 60) {
  line <- paste(rep(title_char, width), collapse = "")
  log_line(line)
}

# -------------------------------------------------------------------
# Assertions
# -------------------------------------------------------------------
assert_files_exist <- function(paths, label) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(
      paste0(
        label, " missing:\n  - ",
        paste(normalizePath(missing, winslash = "/", mustWork = FALSE), collapse = "\n  - ")
      ),
      call. = FALSE
    )
  }
}

# -------------------------------------------------------------------
# Required static inputs
# -------------------------------------------------------------------
INPUT_FILES <- c(
  file.path(DIR_ROOT, "config_analysis.R"),
  file.path(DIR_ROOT, "plot_style.R"),
  file.path(DIR_ROOT, "run_all.R"),
  file.path(DIR_ROOT, "01_preprocessing.R"),
  file.path(DIR_ROOT, "02_pillar1_burden.R"),
  file.path(DIR_ROOT, "03_pillar2_topology.R"),
  file.path(DIR_ROOT, "04_pillar3_anchor.R"),
  file.path(DIR_ROOT, "05_pillar4_concentration.R"),
  file.path(DIR_ROOT, "06_pillar5_depth.R"),
  file.path(DIR_ROOT, "07_pillar6_sourcesink.R"),
  file.path(DIR_ROOT, "08_presentation.R"),
  file.path(DIR_ROOT, "dnk_source_snapshot.csv"),
  file.path(DIR_ROOT, "dnk_source_table.csv"),
  file.path(DIR_ROOT, "branch_ages.csv"),
  file.path(DIR_ROOT, "geocoded_locations_cache.rds"),
  file.path(DIR_ROOT, "ph908_branch_map.tsv"),
  file.path(DIR_ROOT, "ph908_node_age_tiers.tsv")
)

# -------------------------------------------------------------------
# Ensure output directories exist
# -------------------------------------------------------------------
dir.create(DIR_RESULTS, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_TABLES, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIGUREDATA, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FACTS, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIGURES, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_SUPP, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DIR_SUPP, "tables_preview_html"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DIR_SUPP, "tables_preview_csv"), showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------
# Step registry
# -------------------------------------------------------------------
PIPELINE <- list(
  list(
    id = "01",
    label = "Preprocessing",
    script = file.path(DIR_ROOT, "01_preprocessing.R"),
    outputs = c(
      file.path(DIR_RESULTS, "results_preprocessing.rds"),
      file.path(DIR_RESULTS, "preprocessed_dataset.rds"),
      file.path(DIR_FACTS, "facts_preprocessing.txt"),
      file.path(DIR_TABLES, "table_preprocessing_supp.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_preprocessing_flow.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_preprocessing_geocodes.csv")
    )
  ),
  list(
    id = "02",
    label = "Pillar 1 burden-capacity bound",
    script = file.path(DIR_ROOT, "02_pillar1_burden.R"),
    outputs = c(
      file.path(DIR_RESULTS, "results_p1.rds"),
      file.path(DIR_FACTS, "facts_p1.txt"),
      file.path(DIR_TABLES, "table_p1_supp.csv"),
      file.path(DIR_TABLES, "table_p1_rx.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p1_main.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p1_s2.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p1_s3.csv")
    )
  ),
  list(
    id = "03",
    label = "Pillar 2 topology",
    script = file.path(DIR_ROOT, "03_pillar2_topology.R"),
    outputs = c(
      file.path(DIR_RESULTS, "results_p2.rds"),
      file.path(DIR_FACTS, "facts_p2.txt"),
      file.path(DIR_TABLES, "table_p2_supp.csv"),
      file.path(DIR_TABLES, "table_p2_rx.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p2_main.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p2_null.csv")
    )
  ),
  list(
    id = "04",
    label = "Pillar 3 anchor-packet occupancy",
    script = file.path(DIR_ROOT, "04_pillar3_anchor.R"),
    outputs = c(
      file.path(DIR_RESULTS, "results_p3.rds"),
      file.path(DIR_FACTS, "facts_p3.txt"),
      file.path(DIR_TABLES, "table_p3_supp.csv"),
      file.path(DIR_TABLES, "table_p3_rx.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p3_main.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p3_sens.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p3_supp.csv")
    )
  ),
  list(
    id = "05",
    label = "Pillar 4 concentration",
    script = file.path(DIR_ROOT, "05_pillar4_concentration.R"),
    outputs = c(
      file.path(DIR_RESULTS, "results_p4.rds"),
      file.path(DIR_FACTS, "facts_p4.txt"),
      file.path(DIR_TABLES, "table_p4_supp.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p4_main.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p4_sens.csv")
    )
  ),
  list(
    id = "06",
    label = "Pillar 5 depth coherence",
    script = file.path(DIR_ROOT, "06_pillar5_depth.R"),
    outputs = c(
      file.path(DIR_RESULTS, "results_p5.rds"),
      file.path(DIR_FACTS, "facts_p5.txt"),
      file.path(DIR_TABLES, "table_p5_supp.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p5_main.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p5_inset.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p5_null.csv")
    )
  ),
  list(
    id = "07",
    label = "Pillar 6 source-sink inversion",
    script = file.path(DIR_ROOT, "07_pillar6_sourcesink.R"),
    outputs = c(
      file.path(DIR_RESULTS, "results_p6.rds"),
      file.path(DIR_FACTS, "facts_p6.txt"),
      file.path(DIR_TABLES, "table_p6_supp.csv"),
      file.path(DIR_TABLES, "table_p6_threshold_sensitivity.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p6_topk.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p6_foundercarry.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p6_jackknife.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p6_threshold_grid.csv"),
      file.path(DIR_FIGUREDATA, "figuredata_p6_null.csv")
    )
  ),
  list(
    id = "08",
    label = "Presentation and supplement rendering",
    script = file.path(DIR_ROOT, "08_presentation.R"),
    outputs = c(
      file.path(DIR_FIGURES, "Figure_Burden_Main.pdf"),
      file.path(DIR_FIGURES, "Figure_Topology_Main.pdf"),
      file.path(DIR_FIGURES, "Figure_Anchor_Main.pdf"),
      file.path(DIR_FIGURES, "Figure_Concentration_Main.pdf"),
      file.path(DIR_FIGURES, "Figure_Depth_Main.pdf"),
      file.path(DIR_FIGURES, "Figure_SourceSink_Main.pdf"),
      file.path(DIR_SUPP, "tables_preview_html", "Table_S1_preprocessing_summary_metrics.html"),
      file.path(DIR_SUPP, "tables_preview_html", "Table_S12_influential_unit_topology.html")
    )
  )
)

# -------------------------------------------------------------------
# Summary object
# -------------------------------------------------------------------
step_results <- vector("list", length(PIPELINE))
all_warning_lines <- character(0)

# -------------------------------------------------------------------
# Step runner
# -------------------------------------------------------------------
run_step <- function(step, strict_output_checks = TRUE) {
  script_path <- step$script
  expected_outputs <- step$outputs
  step_id <- step$id
  step_label <- step$label
  step_name <- basename(script_path)
  
  write_block("=")
  log_line("STEP ", step_id, " — ", step_label)
  log_line("Script: ", step_name)
  log_line("Path: ", script_path)
  write_block("-")
  
  if (!file.exists(script_path)) {
    stop(paste0("Script not found: ", script_path), call. = FALSE)
  }
  
  start_time <- Sys.time()
  warning_bucket <- character(0)
  
  # Duplicate regular output to the log file while keeping console visible.
  log_con <- file(LOG_FILE, open = "at")
  sink(log_con, type = "output", split = TRUE)
  
  on.exit({
    while (sink.number(type = "output") > 0) sink(type = "output")
    try(close(log_con), silent = TRUE)
  }, add = TRUE)
  
  eval_env <- new.env(parent = globalenv())
  
  err <- tryCatch(
    withCallingHandlers(
      {
        sys.source(script_path, envir = eval_env, chdir = TRUE)
        NULL
      },
      warning = function(w) {
        warning_bucket <<- c(warning_bucket, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  
  elapsed_sec <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
  
  output_check_ok <- TRUE
  output_check_msg <- "OK"
  
  if (!inherits(err, "error") && strict_output_checks && length(expected_outputs) > 0) {
    missing_outputs <- expected_outputs[!file.exists(expected_outputs)]
    if (length(missing_outputs) > 0) {
      output_check_ok <- FALSE
      output_check_msg <- paste(
        c("Missing expected outputs:", paste0("  - ", normalizePath(missing_outputs, winslash = "/", mustWork = FALSE))),
        collapse = "\n"
      )
    }
  }
  
  status <- if (inherits(err, "error")) "FAILED" else if (!output_check_ok) "FAILED" else "OK"
  
  log_line("")
  log_line("Step status: ", status)
  log_line("Elapsed seconds: ", elapsed_sec)
  
  if (length(unique(warning_bucket)) > 0) {
    log_line("Warnings captured: ", length(unique(warning_bucket)))
    for (w in unique(warning_bucket)) {
      log_line("  WARNING: ", w)
    }
  } else {
    log_line("Warnings captured: 0")
  }
  
  if (!output_check_ok) {
    log_line(output_check_msg)
  }
  
  if (inherits(err, "error")) {
    log_line("ERROR: ", conditionMessage(err))
  }
  
  write_block("=")
  log_line("")
  
  list(
    step_id = step_id,
    label = step_label,
    script = step_name,
    status = status,
    elapsed_sec = elapsed_sec,
    warnings_n = length(unique(warning_bucket)),
    warnings_text = paste(unique(warning_bucket), collapse = " || "),
    error_text = if (inherits(err, "error")) conditionMessage(err) else "",
    output_check_ok = output_check_ok,
    output_check_msg = output_check_msg
  )
}

# -------------------------------------------------------------------
# Header
# -------------------------------------------------------------------
writeLines(character(0), LOG_FILE)

write_block("=")
log_line("PH908 REPRODUCIBLE ANALYSIS PIPELINE")
write_block("=")
log_line("Root: ", DIR_ROOT)
log_line("Config: ", CONFIG_FILE)
log_line("Seed: ", SEED_MAIN)
log_line("Continue on error: ", CONTINUE_ON_ERROR)
log_line("Strict output checks: ", STRICT_OUTPUT_CHECKS)
log_line("Run started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log_line("Log file: ", LOG_FILE)
log_line("Summary CSV: ", SUMMARY_CSV)
log_line("Warnings TXT: ", WARNINGS_TXT)
write_block("=")
log_line("")

assert_files_exist(INPUT_FILES, "Required pipeline inputs")
log_line("Static input check: OK")
log_line("")

# -------------------------------------------------------------------
# Execute pipeline
# -------------------------------------------------------------------
pipeline_start <- Sys.time()

for (i in seq_along(PIPELINE)) {
  res <- run_step(PIPELINE[[i]], strict_output_checks = STRICT_OUTPUT_CHECKS)
  step_results[[i]] <- res
  
  if (nzchar(res$warnings_text)) {
    all_warning_lines <- c(
      all_warning_lines,
      paste0("[", res$step_id, " ", res$script, "] ", gsub(" \\|\\| ", "\n[continued] ", res$warnings_text))
    )
  }
  
  if (identical(res$status, "FAILED") && !CONTINUE_ON_ERROR) {
    break
  }
}

total_elapsed_min <- round(as.numeric(difftime(Sys.time(), pipeline_start, units = "mins")), 2)

# -------------------------------------------------------------------
# Finalize summary
# -------------------------------------------------------------------
summary_df <- do.call(
  rbind,
  lapply(step_results[!vapply(step_results, is.null, logical(1))], as.data.frame)
)

if (is.null(summary_df) || nrow(summary_df) == 0) {
  summary_df <- data.frame(
    step_id = character(),
    label = character(),
    script = character(),
    status = character(),
    elapsed_sec = numeric(),
    warnings_n = integer(),
    warnings_text = character(),
    error_text = character(),
    output_check_ok = logical(),
    output_check_msg = character(),
    stringsAsFactors = FALSE
  )
}

utils::write.csv(summary_df, SUMMARY_CSV, row.names = FALSE, na = "")

if (length(all_warning_lines) == 0) {
  writeLines("No warnings captured during this run.", WARNINGS_TXT, useBytes = TRUE)
} else {
  writeLines(unique(all_warning_lines), WARNINGS_TXT, useBytes = TRUE)
}

n_ok <- sum(summary_df$status == "OK", na.rm = TRUE)
n_failed <- sum(summary_df$status == "FAILED", na.rm = TRUE)

write_block("=")
log_line("PIPELINE SUMMARY")
write_block("=")
log_line("Run finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log_line("Elapsed minutes: ", total_elapsed_min)
log_line("Successful steps: ", n_ok)
log_line("Failed steps: ", n_failed)
log_line("Summary CSV: ", SUMMARY_CSV)
log_line("Warnings TXT: ", WARNINGS_TXT)
log_line("Results dir: ", DIR_RESULTS)
log_line("Tables dir: ", DIR_TABLES)
log_line("Figuredata dir: ", DIR_FIGUREDATA)
log_line("Figures dir: ", DIR_FIGURES)
log_line("Facts dir: ", DIR_FACTS)
log_line("Supplement dir: ", DIR_SUPP)
write_block("=")

if (n_failed > 0) {
  failed_rows <- summary_df[summary_df$status == "FAILED", , drop = FALSE]
  for (j in seq_len(nrow(failed_rows))) {
    log_line("FAILED STEP ", failed_rows$step_id[j], " — ", failed_rows$script[j])
    if (nzchar(failed_rows$error_text[j])) {
      log_line("  Error: ", failed_rows$error_text[j])
    }
    if (!isTRUE(failed_rows$output_check_ok[j])) {
      log_line("  Output check: ", failed_rows$output_check_msg[j])
    }
  }
  stop("Pipeline finished with one or more failed steps. See log and summary files.", call. = FALSE)
} else {
  log_line("Pipeline completed successfully.")
}