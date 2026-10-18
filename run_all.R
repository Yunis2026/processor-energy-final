# Reproduce all analyses in a clean temporary project before updating outputs.
# Run with: Rscript run_all.R
# Requires base R only.
# Existing source files and paper_data.csv are not changed.

run_all <- function() {

  # Resolve the original project directory.
  file_arg <- grep(
    "^--file=",
    commandArgs(FALSE),
    value = TRUE
  )

  project_dir <- if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    getwd()
  }

  previous_dir <- getwd()
  on.exit(setwd(previous_dir), add = TRUE)

  scripts <- c(
    "pde_model.R",
    "sensitivity_analysis.R",
    "strategy_map.R",
    "input_robustness.R"
  )

  inputs <- c(
    "paper_data.csv",
    scripts
  )

  missing_files <- inputs[
    !file.exists(file.path(project_dir, inputs))
  ]

  if (length(missing_files)) {
    stop(
      "Missing project files: ",
      paste(missing_files, collapse = ", ")
    )
  }

  # Use a clean workspace so old results cannot count as new outputs.
  stage_dir <- tempfile(
    "processor_energy_reproduce_"
  )

  dir.create(stage_dir)

  copied <- file.copy(
    file.path(project_dir, inputs),
    stage_dir
  )

  if (!all(copied)) {
    stop(
      "Could not copy inputs into the temporary project."
    )
  }

  setwd(stage_dir)

  dir.create(
    file.path("outputs", "logs"),
    recursive = TRUE
  )

  # Prefer Rscript.
  # A direct R executable supports minimal Unix installations.
  rscript <- Sys.which("Rscript")

  r_binary <- file.path(
    R.home("bin"),
    "exec",
    "R"
  )

  if (
    !nzchar(rscript) &&
    !file.exists(r_binary)
  ) {
    stop(
      "Could not locate an Rscript or R executable."
    )
  }

  expected <- list(
    c(
      "pde/summary.csv",
      "pde/convergence.csv",
      "pde/energy_comparison.png"
    ),

    c(
      "sensitivity/oat_results.csv",
      "sensitivity/factor_screening.csv",
      "sensitivity/oat_saving_curves.png"
    ),

    c(
      "strategy_map/strategy_map.csv",
      "strategy_map/break_even.csv",
      "strategy_map/strategy_map.png"
    ),

    c(
      "input_robustness/current_oat.csv",
      "input_robustness/current_break_even.csv",
      "input_robustness/joint_corner_summary.csv",
      "input_robustness/current_perturbation.png"
    )
  )

  run_rows <- list()

  # Run each analysis in a separate R process.
  for (i in seq_along(scripts)) {

    script <- scripts[i]

    cat(sprintf(
      "\n[%d/%d] Running %s...\n",
      i,
      length(scripts),
      script
    ))

    log_file <- file.path(
      "outputs",
      "logs",
      paste0(script, ".log")
    )

    start_time <- proc.time()[["elapsed"]]

    if (nzchar(rscript)) {

      status <- system2(
        rscript,

        args = c(
          "--vanilla",
          shQuote(script)
        ),

        stdout = log_file,
        stderr = log_file
      )

    } else {

      status <- system2(
        r_binary,

        args = c(
          "--vanilla",
          "--slave",
          paste0("--file=", shQuote(script))
        ),

        stdout = log_file,
        stderr = log_file
      )
    }

    log_lines <- readLines(
      log_file,
      warn = FALSE
    )

    cat(
      paste(
        tail(log_lines, 8),
        collapse = "\n"
      ),
      "\n"
    )

    if (!identical(
      as.integer(status),
      0L
    )) {
      stop(
        "Analysis failed. Diagnostic files remain in: ",
        stage_dir
      )
    }

    generated <- file.path(
      "outputs",
      expected[[i]]
    )

    if (
      !all(file.exists(generated)) ||
      any(file.info(generated)$size <= 0)
    ) {
      stop(
        "Expected outputs were not generated. ",
        "Check the final function call in ",
        script,
        ". Temporary project: ",
        stage_dir
      )
    }

    run_rows[[i]] <- data.frame(
      script = script,

      exit_status = status,

      elapsed_seconds =
        proc.time()[["elapsed"]] - start_time,

      expected_outputs_present = TRUE
    )
  }

  # Read the newly generated tables.
  read_result <- function(path) {
    read.csv(
      file.path("outputs", path),
      stringsAsFactors = FALSE
    )
  }

  summary_table <- read_result(
    "pde/summary.csv"
  )

  oat <- read_result(
    "sensitivity/oat_results.csv"
  )

  map <- read_result(
    "strategy_map/strategy_map.csv"
  )

  current_oat <- read_result(
    "input_robustness/current_oat.csv"
  )

  corners <- read_result(
    "input_robustness/joint_corners.csv"
  )

  # Check expected scenario counts and successful status.
  stopifnot(
    nrow(summary_table) == 4L,
    nrow(oat) == 28L,
    nrow(map) == 91L,
    nrow(current_oat) == 14L,
    nrow(corners) == 12L,

    all(oat$status == "ok"),
    all(map$status == "ok")
  )

  # Check that all reported energies are finite and nonnegative.
  for (
    table in list(
      summary_table,
      oat,
      map,
      current_oat,
      corners
    )
  ) {

    energy_columns <- intersect(
      c(
        "total_energy_J",
        "high_energy_J",
        "low_energy_J"
      ),

      names(table)
    )

    energies <- unlist(
      table[energy_columns],
      use.names = FALSE
    )

    stopifnot(
      all(is.finite(energies)),
      all(energies >= 0)
    )
  }

  # All scripts must reproduce the same coupled baseline.
  reference_high <- summary_table$total_energy_J[
    summary_table$mode == "high" &
      summary_table$feedback
  ]

  reference_low <- summary_table$total_energy_J[
    summary_table$mode == "low" &
      summary_table$feedback
  ]

  stopifnot(
    length(reference_high) == 1L,
    length(reference_low) == 1L
  )

  baseline_files <- c(
    "sensitivity/baseline_check.csv",
    "strategy_map/baseline.csv",
    "input_robustness/baseline.csv"
  )

  for (path in baseline_files) {

    check <- read_result(path)

    if ("feedback" %in% names(check)) {
      check <- check[check$feedback, ]
    }

    stopifnot(
      nrow(check) == 1L,

      abs(
        check$high_energy_J -
          reference_high
      ) < 1e-8,

      abs(
        check$low_energy_J -
          reference_low
      ) < 1e-8
    )
  }

  # Check energy agreement at successfully calculated crossings.
  crossing_files <- c(
    "strategy_map/break_even.csv",
    "input_robustness/current_break_even.csv"
  )

  for (path in crossing_files) {

    roots <- read_result(path)

    found <- roots$status == "crossing_found"

    stopifnot(
      all(
        is.finite(
          roots$root_energy_difference_J[found]
        )
      ),

      all(
        abs(
          roots$root_energy_difference_J[found]
        ) <= 1e-6
      )
    )
  }

  write.csv(
    do.call(rbind, run_rows),

    file.path(
      "outputs",
      "run_manifest.csv"
    ),

    row.names = FALSE
  )

  writeLines(
    capture.output(sessionInfo()),

    file.path(
      "outputs",
      "run_session_info.txt"
    )
  )

  # Update project outputs only after all analyses and checks pass.
  generated_files <- list.files(
    "outputs",

    recursive = TRUE,
    all.files = TRUE,
    full.names = FALSE
  )

  generated_files <- generated_files[
    !file.info(
      file.path(
        "outputs",
        generated_files
      )
    )$isdir
  ]

  for (relative_path in generated_files) {

    destination <- file.path(
      project_dir,
      "outputs",
      relative_path
    )

    dir.create(
      dirname(destination),

      recursive = TRUE,
      showWarnings = FALSE
    )

    copied <- file.copy(
      file.path(
        "outputs",
        relative_path
      ),

      destination,

      overwrite = TRUE
    )

    if (!copied) {
      stop(
        "Could not update output: ",
        destination,
        ". Regenerated files remain in: ",
        stage_dir
      )
    }
  }

  cat(
    "\nAll four analyses and output checks passed.\n"
  )

  cat(
    "Regenerated results saved to: ",
    file.path(project_dir, "outputs"),
    "\n",
    sep = ""
  )

  cat(
    paste(
      "This verifies reproducibility and numerical consistency,",
      "not hardware accuracy.\n"
    )
  )

  invisible(
    do.call(rbind, run_rows)
  )
}

run_all()