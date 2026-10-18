# Four-factor one-at-a-time sensitivity analysis.
# This script uses the existing PDE solver without modifying pde_model.R.
# Run with: Rscript sensitivity_analysis.R
# No external R packages are required.

run_sensitivity_analysis <- function() {

  # ------------------------------------------------------------
  # 1. Resolve the project directory
  # ------------------------------------------------------------

  file_arg <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )

  script_dir <- if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    getwd()
  }

  previous_dir <- getwd()
  on.exit(setwd(previous_dir), add = TRUE)
  setwd(script_dir)

  required_files <- c("pde_model.R", "paper_data.csv")

  for (path in required_files) {
    if (!file.exists(path)) {
      stop("Required file not found: ", path)
    }
  }

  # ------------------------------------------------------------
  # 2. Load the existing model definitions
  # ------------------------------------------------------------

  # A private environment isolates the model parameters.
  model <- new.env(parent = globalenv())

  model$pde_definitions_only <- TRUE

  # Read initialization and solver definitions only.
  # Stop after simulate_mode is defined, before the original driver.
  expressions <- parse(
    file = "pde_model.R",
    keep.source = FALSE
  )

  for (expression in expressions) {

    eval(expression, envir = model)

    if (exists(
      "simulate_mode",
      envir = model,
      inherits = FALSE
    )) {
      break
    }
  }

  required_definitions <- c(
    "simulate_mode",
    "alpha_feedback",
    "beta_per_K",
    "R_th_K_per_W",
    "idle_reference_power_W"
  )

  definitions_available <- vapply(
    required_definitions,
    exists,
    logical(1),
    envir = model,
    inherits = FALSE
  )

  if (!all(definitions_available)) {
    stop("Required model definitions are missing.")
  }

  if (exists(
    "scenario_results",
    envir = model,
    inherits = FALSE
  )) {
    stop(
      "The model driver ran unexpectedly; check the model file structure."
    )
  }

  # ------------------------------------------------------------
  # 3. Create the output directory
  # ------------------------------------------------------------

  output_dir <- file.path("outputs", "sensitivity")

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  writeLines(
    capture.output(sessionInfo()),
    file.path(output_dir, "session_info.txt")
  )

  # ------------------------------------------------------------
  # 4. Define the four sensitivity factors
  # ------------------------------------------------------------

  # These are exploratory ranges, not measured confidence intervals.
  factors <- data.frame(
    factor = c(
      "idle_power_mW",
      "thermal_resistance_K_per_W",
      "feedback_alpha",
      "feedback_beta_per_K"
    ),

    label = c(
      "Idle power (mW)",
      "Thermal resistance (K/W)",
      "Feedback alpha",
      "Feedback beta (1/K)"
    ),

    baseline = c(
      model$idle_reference_power_W * 1000,
      model$R_th_K_per_W,
      model$alpha_feedback,
      model$beta_per_K
    ),

    lower = c(
      0,
      30,
      0,
      0.025
    ),

    upper = c(
      0.30,
      120,
      0.40,
      0.10
    ),

    status = "exploratory assumption; not a confidence interval",

    stringsAsFactors = FALSE
  )

  baseline_parameters <- setNames(
    factors$baseline,
    factors$factor
  )

  nx <- 21L
  dt <- 1

  # ------------------------------------------------------------
  # 5. Update model parameters consistently
  # ------------------------------------------------------------

  set_parameters <- function(parameters) {

    model$idle_reference_power_W <-
      parameters[["idle_power_mW"]] / 1000

    model$R_th_K_per_W <-
      parameters[["thermal_resistance_K_per_W"]]

    # Recompute distributed cooling whenever thermal resistance changes.
    model$h_v_W_per_m3K <-
      1 / (model$R_th_K_per_W * model$volume_m3)

    model$alpha_feedback <-
      parameters[["feedback_alpha"]]

    model$beta_per_K <-
      parameters[["feedback_beta_per_K"]]
  }

  # Reuse identical scenarios instead of running them repeatedly.
  cache <- new.env(parent = emptyenv())

  # ------------------------------------------------------------
  # 6. Evaluate both operating strategies
  # ------------------------------------------------------------

  evaluate_pair <- function(
    parameters,
    feedback = TRUE,
    grid_cells = nx,
    time_step = dt
  ) {

    key <- paste(
      c(
        sprintf("%.16g", parameters),
        feedback,
        grid_cells,
        time_step
      ),
      collapse = "|"
    )

    if (exists(key, envir = cache, inherits = FALSE)) {
      return(cache[[key]])
    }

    # Reset all four parameters before every new scenario.
    set_parameters(parameters)

    result <- tryCatch({

      high <- model$simulate_mode(
        "high",
        nx = grid_cells,
        dt_nominal_s = time_step,
        feedback = feedback
      )

      low <- model$simulate_mode(
        "low",
        nx = grid_cells,
        dt_nominal_s = time_step,
        feedback = feedback
      )

      residual <- max(abs(c(
        high$energy_balance_residual_J,
        low$energy_balance_residual_J
      )))

      if (!is.finite(residual) || residual > 1e-7) {
        stop("Energy-balance residual exceeds 1e-7 J.")
      }

      data.frame(
        status = "ok",
        error = "",

        nx = grid_cells,
        dt_s = time_step,
        feedback = feedback,

        high_energy_J = high$total_energy_J,
        low_energy_J = low$total_energy_J,

        saving_percent =
          100 *
          (low$total_energy_J - high$total_energy_J) /
          low$total_energy_J,

        energy_difference_high_minus_low_J =
          high$total_energy_J - low$total_energy_J,

        high_peak_temperature_C =
          high$peak_temperature_C,

        low_peak_temperature_C =
          low$peak_temperature_C,

        high_carbon_mg_CO2e =
          high$carbon_mg_CO2e,

        low_carbon_mg_CO2e =
          low$carbon_mg_CO2e,

        max_energy_balance_residual_J = residual,

        stringsAsFactors = FALSE
      )

    }, error = function(e) {

      data.frame(
        status = "failed",
        error = conditionMessage(e),

        nx = grid_cells,
        dt_s = time_step,
        feedback = feedback,

        high_energy_J = NA_real_,
        low_energy_J = NA_real_,
        saving_percent = NA_real_,

        energy_difference_high_minus_low_J = NA_real_,

        high_peak_temperature_C = NA_real_,
        low_peak_temperature_C = NA_real_,

        high_carbon_mg_CO2e = NA_real_,
        low_carbon_mg_CO2e = NA_real_,

        max_energy_balance_residual_J = NA_real_,

        stringsAsFactors = FALSE
      )
    })

    cache[[key]] <- result

    result
  }

  # ------------------------------------------------------------
  # 7. Check the no-feedback baseline
  # ------------------------------------------------------------

  cat("Checking the no-feedback energy baseline...\n")

  no_feedback <- evaluate_pair(
    baseline_parameters,
    feedback = FALSE
  )

  if (no_feedback$status != "ok") {
    stop(no_feedback$error)
  }

  active_time_s <-
    model$task_count * model$paper_data$task_time_s

  analytic_energy_J <-
    model$paper_data$power_W * active_time_s +
    model$idle_reference_power_W *
    (model$observation_time_s - active_time_s)

  stopifnot(
    abs(
      no_feedback$high_energy_J -
      analytic_energy_J[model$high_row]
    ) < 1e-8,

    abs(
      no_feedback$low_energy_J -
      analytic_energy_J[model$low_row]
    ) < 1e-8
  )

  # ------------------------------------------------------------
  # 8. Calculate the coupled baseline
  # ------------------------------------------------------------

  cat("Computing the coupled baseline...\n")

  baseline <- evaluate_pair(baseline_parameters)

  if (baseline$status != "ok") {
    stop(baseline$error)
  }

  write.csv(
    rbind(no_feedback, baseline),
    file.path(output_dir, "baseline_check.csv"),
    row.names = FALSE
  )

  write.csv(
    factors,
    file.path(output_dir, "factor_ranges.csv"),
    row.names = FALSE
  )

  # ------------------------------------------------------------
  # 9. Run the one-at-a-time sensitivity analysis
  # ------------------------------------------------------------

  rows <- list()

  for (i in seq_len(nrow(factors))) {

    values <- seq(
      factors$lower[i],
      factors$upper[i],
      length.out = 7
    )

    # Use the exact baseline value at the matching level.
    near_baseline <-
      abs(values - factors$baseline[i]) < 1e-12

    values[near_baseline] <- factors$baseline[i]

    for (value in values) {

      parameters <- baseline_parameters

      parameters[[factors$factor[i]]] <- value

      cat(sprintf(
        "%s = %.6g\n",
        factors$factor[i],
        value
      ))

      result <- evaluate_pair(parameters)

      rows[[length(rows) + 1L]] <- cbind(
        data.frame(
          factor = factors$factor[i],
          factor_value = value,

          is_baseline =
            abs(value - factors$baseline[i]) < 1e-12,

          stringsAsFactors = FALSE
        ),

        as.data.frame(as.list(parameters)),

        result
      )

      # Save after each case to retain progress if a later case fails.
      write.csv(
        do.call(rbind, rows),
        file.path(output_dir, "oat_results.csv"),
        row.names = FALSE
      )
    }
  }

  results <- do.call(rbind, rows)

  if (any(results$status != "ok")) {
    stop(
      paste(
        "Some cases failed.",
        "Inspect oat_results.csv;",
        "no ranking was produced."
      )
    )
  }

  # Four factors with seven levels produce 28 exported rows.
  stopifnot(nrow(results) == 28L)

  # The baseline must be identical in all four factor sweeps.
  stopifnot(
    all(
      abs(
        results$saving_percent[results$is_baseline] -
        baseline$saving_percent
      ) < 1e-9
    )
  )

  # ------------------------------------------------------------
  # 10. Verify the alpha = 0 limiting case
  # ------------------------------------------------------------

  alpha_zero <- results[
    results$factor == "feedback_alpha" &
    results$factor_value == 0,
  ]

  stopifnot(
    abs(
      alpha_zero$high_energy_J -
      no_feedback$high_energy_J
    ) < 1e-8,

    abs(
      alpha_zero$low_energy_J -
      no_feedback$low_energy_J
    ) < 1e-8
  )

  # ------------------------------------------------------------
  # 11. Summarize factor influence
  # ------------------------------------------------------------

  # The response span depends on the selected exploration range.
  # It is not a global sensitivity index.
  screening <- do.call(
    rbind,

    lapply(factors$factor, function(name) {

      selected <- results[results$factor == name, ]

      data.frame(
        factor = name,

        minimum_saving_percent =
          min(selected$saving_percent),

        maximum_saving_percent =
          max(selected$saving_percent),

        span_percentage_points =
          diff(range(selected$saving_percent)),

        both_strategies_win_at_sampled_levels =
          any(selected$saving_percent < 0) &&
          any(selected$saving_percent > 0)
      )
    })
  )

  screening <- screening[
    order(-screening$span_percentage_points),
  ]

  write.csv(
    screening,
    file.path(output_dir, "factor_screening.csv"),
    row.names = FALSE
  )

  # ------------------------------------------------------------
  # 12. Check the baseline under time-step refinement
  # ------------------------------------------------------------

  # This checks baseline integrated outputs only.
  # It does not prove convergence of every sensitivity scenario.
  cat("Checking coupled baseline at dt = 0.5 s...\n")

  refined <- evaluate_pair(
    baseline_parameters,
    time_step = 0.5
  )

  if (refined$status != "ok") {
    stop(refined$error)
  }

  time_step_checks <- rbind(baseline, refined)

  time_step_checks$saving_change_from_dt_1_percentage_points <-
    time_step_checks$saving_percent -
    baseline$saving_percent

  write.csv(
    time_step_checks,
    file.path(output_dir, "time_step_check.csv"),
    row.names = FALSE
  )

  # ------------------------------------------------------------
  # 13. Plot the four sensitivity curves
  # ------------------------------------------------------------

  png(
    file.path(output_dir, "oat_saving_curves.png"),
    width = 2000,
    height = 1500,
    res = 200
  )

  par(
    mfrow = c(2, 2),
    mar = c(4.6, 4.8, 3, 1.2),
    oma = c(1, 0, 1, 0)
  )

  # Use a common vertical scale to compare factor effects.
  response_range <- range(c(results$saving_percent, 0))

  padding <- max(
    diff(response_range) * 0.08,
    0.1
  )

  for (i in seq_len(nrow(factors))) {

    selected <- results[
      results$factor == factors$factor[i],
    ]

    plot(
      selected$factor_value,
      selected$saving_percent,

      type = "b",
      pch = 19,
      col = "#245A81",
      lwd = 2,

      xlab = factors$label[i],
      ylab = "High-mode energy saving (%)",
      main = factors$label[i],

      ylim = response_range + c(-padding, padding)
    )

    abline(
      h = 0,
      col = "gray40",
      lty = 2
    )

    abline(
      v = factors$baseline[i],
      col = "gray60",
      lty = 3
    )

    points(
      factors$baseline[i],
      baseline$saving_percent,

      pch = 21,
      bg = "#E69F00",
      col = "black",
      cex = 1.3
    )
  }

  mtext(
    paste(
      "Positive: high + idle wins;",
      "negative: low wins.",
      "Orange: baseline."
    ),
    side = 1,
    outer = TRUE,
    line = -0.2,
    cex = 0.85
  )

  dev.off()

  # ------------------------------------------------------------
  # 14. Plot the response-span comparison
  # ------------------------------------------------------------

  png(
    file.path(output_dir, "factor_screening.png"),
    width = 1800,
    height = 1000,
    res = 200
  )

  par(mar = c(5, 11, 4, 2))

  factor_labels <- factors$label[
    match(screening$factor, factors$factor)
  ]

  barplot(
    rev(screening$span_percentage_points),

    names.arg = rev(factor_labels),
    horiz = TRUE,
    las = 1,

    col = "#245A81",
    border = NA,

    xlab = "Range of energy saving (percentage points)",

    main = "Factor influence within selected exploration ranges"
  )

  dev.off()

  # ------------------------------------------------------------
  # 15. Save interpretation notes
  # ------------------------------------------------------------

  notes <- c(
    "Four-factor coupled PDE screening",

    sprintf(
      "Coupled baseline high energy: %.9f J",
      baseline$high_energy_J
    ),

    sprintf(
      "Coupled baseline low energy: %.9f J",
      baseline$low_energy_J
    ),

    sprintf(
      "Coupled baseline high-mode saving: %.6f%%",
      baseline$saving_percent
    ),

    "Seven levels per factor; 28 rows.",
    "Repeated baseline calculations are cached.",

    paste(
      "Only one factor varies at a time.",
      "Interactions are not measured here."
    ),

    paste(
      "Factor ranking depends on selected ranges.",
      "The ranges are not confidence intervals."
    ),

    paste(
      "No measurement uncertainty distributions",
      "or hardware calibration are claimed."
    ),

    paste(
      "Common fixed positive carbon intensity",
      "preserves the energy ranking."
    ),

    paste(
      "The time-step check covers the baseline.",
      "Inspect sensitive cases before final conclusions."
    ),

    paste(
      "Solver failures are recorded explicitly",
      "and are not evidence of physical thermal runaway."
    )
  )

  writeLines(
    notes,
    file.path(output_dir, "analysis_notes.txt")
  )

  # ------------------------------------------------------------
  # 16. Print the main results
  # ------------------------------------------------------------

  cat("\nSensitivity screening summary:\n")

  print(screening, row.names = FALSE)

  cat(sprintf(
    "\nBaseline high-mode saving: %.6f%%\n",
    baseline$saving_percent
  ))

  cat("Results written to outputs/sensitivity/\n")

  invisible(results)
}

run_sensitivity_analysis()