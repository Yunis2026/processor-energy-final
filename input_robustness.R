# Exploratory current-input perturbations for the coupled processor PDE model.
# Run with: Rscript input_robustness.R
# Requires pde_model.R and paper_data.csv.
# No extra R packages are needed.

run_input_robustness <- function() {

  # ------------------------------------------------------------
  # 1. Resolve the project directory
  # ------------------------------------------------------------

  file_arg <- grep(
    "^--file=",
    commandArgs(FALSE),
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

  required_files <- c(
    "pde_model.R",
    "paper_data.csv"
  )

  if (!all(file.exists(required_files))) {
    stop(
      "Keep input_robustness.R, pde_model.R and paper_data.csv together."
    )
  }

  # ------------------------------------------------------------
  # 2. Load the existing PDE solver
  # ------------------------------------------------------------

  model <- new.env(parent = globalenv())

  model$pde_definitions_only <- TRUE

  expressions <- parse(
    "pde_model.R",
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

  solver_available <- exists(
    "simulate_mode",
    envir = model,
    inherits = FALSE
  )

  driver_executed <- exists(
    "scenario_results",
    envir = model,
    inherits = FALSE
  )

  if (!solver_available || driver_executed) {
    stop(
      "Could not load the solver definitions without running the driver."
    )
  }

  # ------------------------------------------------------------
  # 3. Create the output directory
  # ------------------------------------------------------------

  output_dir <- file.path(
    "outputs",
    "input_robustness"
  )

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
  # 4. Preserve the original inputs
  # ------------------------------------------------------------

  original_data <- model$paper_data

  original_task_count <- model$task_count

  original_idle_W <- model$idle_reference_power_W

  original_resistance <- model$R_th_K_per_W

  original_alpha <- model$alpha_feedback

  original_beta <- model$beta_per_K

  nx <- 21L
  dt_s <- 1

  equality_tolerance_J <- 1e-6

  levels_percent <- c(
    -5, -3, -1, 0, 1, 3, 5
  )

  joint_magnitudes_percent <- c(
    1, 3, 5
  )

  write.csv(
    original_data,
    file.path(output_dir, "input_snapshot.csv"),
    row.names = FALSE
  )

  settings <- data.frame(
    parameter = c(
      "idle_power_W",
      "thermal_resistance_K_per_W",
      "alpha",
      "beta_per_K",
      "nx",
      "dt_s",
      "task_count",
      "observation_time_s",
      "carbon_intensity_kg_CO2e_per_kWh"
    ),

    value = c(
      original_idle_W,
      original_resistance,
      original_alpha,
      original_beta,
      nx,
      dt_s,
      original_task_count,
      model$observation_time_s,
      model$carbon_intensity_kg_CO2e_per_kWh
    )
  )

  write.csv(
    settings,
    file.path(output_dir, "fixed_settings.csv"),
    row.names = FALSE
  )

  # ------------------------------------------------------------
  # 5. Define paired simulations with perturbed currents
  # ------------------------------------------------------------

  cache <- new.env(parent = emptyenv())

  evaluations <- list()

  evaluate_pair <- function(
    high_change_percent,
    low_change_percent,
    time_step = dt_s
  ) {

    key <- paste(
      sprintf(
        "%.16g",
        c(
          high_change_percent,
          low_change_percent,
          time_step
        )
      ),
      collapse = "|"
    )

    if (exists(
      key,
      envir = cache,
      inherits = FALSE
    )) {
      return(cache[[key]])
    }

    if (
      min(high_change_percent, low_change_percent) <= -100
    ) {
      stop(
        "Current changes must be greater than -100 percent."
      )
    }

    # Reset to the original inputs before each scenario.
    # The original CSV is never overwritten.
    model$paper_data <- original_data

    model$paper_data$current_mA[model$high_row] <-
      original_data$current_mA[model$high_row] *
      (1 + high_change_percent / 100)

    model$paper_data$current_mA[model$low_row] <-
      original_data$current_mA[model$low_row] *
      (1 + low_change_percent / 100)

    # Recalculate reference power consistently from P = V * I.
    model$paper_data$current_A <-
      model$paper_data$current_mA / 1000

    model$paper_data$power_W <-
      model$paper_data$voltage_V *
      model$paper_data$current_A

    model$idle_reference_power_W <-
      original_idle_W

    model$R_th_K_per_W <-
      original_resistance

    model$h_v_W_per_m3K <-
      1 / (original_resistance * model$volume_m3)

    model$alpha_feedback <- original_alpha

    model$beta_per_K <- original_beta

    # Workload and execution times remain unchanged.
    stopifnot(
      model$task_count == original_task_count,

      identical(
        model$paper_data$task_time_s,
        original_data$task_time_s
      )
    )

    result <- tryCatch({

      high <- model$simulate_mode(
        "high",
        nx = nx,
        dt_nominal_s = time_step,
        feedback = TRUE
      )

      low <- model$simulate_mode(
        "low",
        nx = nx,
        dt_nominal_s = time_step,
        feedback = TRUE
      )

      gap <-
        high$total_energy_J -
        low$total_energy_J

      residual <- max(abs(c(
        high$energy_balance_residual_J,
        low$energy_balance_residual_J
      )))

      if (
        !is.finite(residual) ||
        residual > 1e-7
      ) {
        stop(
          "Energy-balance residual exceeds 1e-7 J."
        )
      }

      preferred_strategy <- if (
        abs(gap) <= equality_tolerance_J
      ) {
        "approximately_equal"
      } else if (gap < 0) {
        "high_plus_idle"
      } else {
        "low"
      }

      data.frame(
        high_current_change_percent =
          high_change_percent,

        low_current_change_percent =
          low_change_percent,

        high_current_mA =
          model$paper_data$current_mA[model$high_row],

        low_current_mA =
          model$paper_data$current_mA[model$low_row],

        dt_s = time_step,

        high_energy_J =
          high$total_energy_J,

        low_energy_J =
          low$total_energy_J,

        energy_difference_high_minus_low_J =
          gap,

        saving_percent =
          -100 * gap / low$total_energy_J,

        preferred_strategy =
          preferred_strategy,

        high_peak_temperature_C =
          high$peak_temperature_C,

        low_peak_temperature_C =
          low$peak_temperature_C,

        high_carbon_mg_CO2e =
          high$carbon_mg_CO2e,

        low_carbon_mg_CO2e =
          low$carbon_mg_CO2e,

        max_energy_balance_residual_J =
          residual,

        stringsAsFactors = FALSE
      )

    }, error = function(e) {

      failure <- sprintf(
        paste(
          "high change = %.9g%%,",
          "low change = %.9g%%,",
          "dt = %.9g s: %s"
        ),

        high_change_percent,
        low_change_percent,
        time_step,
        conditionMessage(e)
      )

      cat(
        failure,
        "\n",

        file = file.path(
          output_dir,
          "failure_log.txt"
        ),

        append = TRUE
      )

      stop(failure)
    })

    evaluations[[
      length(evaluations) + 1L
    ]] <<- result

    write.csv(
      do.call(rbind, evaluations),

      file.path(
        output_dir,
        "all_evaluations.csv"
      ),

      row.names = FALSE
    )

    cache[[key]] <- result

    result
  }

  # ------------------------------------------------------------
  # 6. Calculate the unperturbed baseline
  # ------------------------------------------------------------

  cat(
    "Computing the unperturbed coupled baseline...\n"
  )

  baseline <- evaluate_pair(0, 0)

  write.csv(
    baseline,
    file.path(output_dir, "baseline.csv"),
    row.names = FALSE
  )

  # ------------------------------------------------------------
  # 7. Change one current at a time
  # ------------------------------------------------------------

  oat_rows <- list()

  for (
    factor in c("high_current", "low_current")
  ) {

    for (change in levels_percent) {

      cat(sprintf(
        "%s change = %+.1f%%...\n",
        factor,
        change
      ))

      result <- if (
        factor == "high_current"
      ) {
        evaluate_pair(change, 0)
      } else {
        evaluate_pair(0, change)
      }

      oat_rows[[
        length(oat_rows) + 1L
      ]] <- cbind(
        data.frame(
          factor = factor,
          change_percent = change
        ),

        result
      )

      write.csv(
        do.call(rbind, oat_rows),

        file.path(
          output_dir,
          "current_oat.csv"
        ),

        row.names = FALSE
      )
    }
  }

  oat <- do.call(rbind, oat_rows)

  stopifnot(nrow(oat) == 14L)

  # ------------------------------------------------------------
  # 8. Evaluate joint perturbation scenarios
  # ------------------------------------------------------------

  # Evaluate four corners at each selected magnitude:
  # high -, low -
  # high -, low +
  # high +, low -
  # high +, low +
  #
  # These are deterministic stress scenarios.
  # No measurement-error distribution is assumed.
  corner_rows <- list()

  for (
    magnitude in joint_magnitudes_percent
  ) {

    for (high_sign in c(-1, 1)) {

      for (low_sign in c(-1, 1)) {

        cat(sprintf(
          "Joint case: high %+.1f%%, low %+.1f%%...\n",

          high_sign * magnitude,

          low_sign * magnitude
        ))

        result <- evaluate_pair(
          high_sign * magnitude,

          low_sign * magnitude
        )

        corner_rows[[
          length(corner_rows) + 1L
        ]] <- cbind(
          data.frame(
            perturbation_magnitude_percent =
              magnitude
          ),

          result
        )

        write.csv(
          do.call(rbind, corner_rows),

          file.path(
            output_dir,
            "joint_corners.csv"
          ),

          row.names = FALSE
        )
      }
    }
  }

  corners <- do.call(rbind, corner_rows)

  stopifnot(nrow(corners) == 12L)

  # ------------------------------------------------------------
  # 9. Summarize the tested joint scenarios
  # ------------------------------------------------------------

  corner_summary <- do.call(
    rbind,

    lapply(
      joint_magnitudes_percent,

      function(magnitude) {

        selected <- corners[
          corners$perturbation_magnitude_percent ==
            magnitude,
        ]

        data.frame(
          perturbation_magnitude_percent =
            magnitude,

          minimum_saving_among_tested_corners_percent =
            min(selected$saving_percent),

          maximum_saving_among_tested_corners_percent =
            max(selected$saving_percent),

          ranking_reverses_among_tested_corners =
            any(
              selected$energy_difference_high_minus_low_J <
                -equality_tolerance_J
            ) &&
            any(
              selected$energy_difference_high_minus_low_J >
                equality_tolerance_J
            )
        )
      }
    )
  )

  write.csv(
    corner_summary,

    file.path(
      output_dir,
      "joint_corner_summary.csv"
    ),

    row.names = FALSE
  )

  # ------------------------------------------------------------
  # 10. Find single-input break-even changes
  # ------------------------------------------------------------

  root_rows <- list()

  for (
    factor in c("high_current", "low_current")
  ) {

    cat(
      "Finding the single-input crossing for ",
      factor,
      "...\n",
      sep = ""
    )

    gap_function <- function(change) {

      result <- if (
        factor == "high_current"
      ) {
        evaluate_pair(change, 0)
      } else {
        evaluate_pair(0, change)
      }

      result$energy_difference_high_minus_low_J
    }

    lower_gap <- gap_function(-5)

    upper_gap <- gap_function(5)

    if (lower_gap == 0) {

      root <- -5

    } else if (upper_gap == 0) {

      root <- 5

    } else if (
      sign(lower_gap) != sign(upper_gap)
    ) {

      root <- uniroot(
        gap_function,

        interval = c(-5, 5),

        f.lower = lower_gap,

        f.upper = upper_gap,

        tol = 1e-6,

        maxiter = 100
      )$root

    } else {

      root <- NA_real_
    }

    root_gap <- if (
      is.finite(root)
    ) {
      gap_function(root)
    } else {
      NA_real_
    }

    if (
      is.finite(root_gap) &&
      abs(root_gap) > equality_tolerance_J
    ) {
      stop(
        "Root energy mismatch exceeds the requested tolerance."
      )
    }

    root_rows[[
      length(root_rows) + 1L
    ]] <- data.frame(
      factor = factor,

      status = if (
        is.finite(root)
      ) {
        "crossing_found"
      } else {
        "no_bracketed_crossing"
      },

      critical_current_change_percent =
        root,

      root_energy_difference_J =
        root_gap
    )
  }

  thresholds <- do.call(
    rbind,
    root_rows
  )

  write.csv(
    thresholds,

    file.path(
      output_dir,
      "current_break_even.csv"
    ),

    row.names = FALSE
  )

  # ------------------------------------------------------------
  # 11. Check crossings with a smaller time step
  # ------------------------------------------------------------

  checks <- list()

  for (
    i in seq_len(nrow(thresholds))
  ) {

    change <-
      thresholds$critical_current_change_percent[i]

    if (!is.finite(change)) {
      next
    }

    refined <- if (
      thresholds$factor[i] == "high_current"
    ) {
      evaluate_pair(
        change,
        0,
        time_step = 0.5
      )
    } else {
      evaluate_pair(
        0,
        change,
        time_step = 0.5
      )
    }

    checks[[
      length(checks) + 1L
    ]] <- cbind(
      data.frame(
        factor = thresholds$factor[i]
      ),

      refined
    )
  }

  if (length(checks)) {

    write.csv(
      do.call(rbind, checks),

      file.path(
        output_dir,
        "crossing_time_step_check.csv"
      ),

      row.names = FALSE
    )
  }

  # ------------------------------------------------------------
  # 12. Plot the single-input perturbation curves
  # ------------------------------------------------------------

  png(
    file.path(
      output_dir,
      "current_perturbation.png"
    ),

    width = 1800,
    height = 1200,
    res = 200
  )

  par(
    mar = c(5, 5, 6.5, 2)
  )

  high_curve <- oat[
    oat$factor == "high_current",
  ]

  low_curve <- oat[
    oat$factor == "low_current",
  ]

  limits <- range(
    c(oat$saving_percent, 0)
  )

  padding <- max(
    diff(limits) * 0.08,
    0.1
  )

  plot(
    high_curve$change_percent,

    high_curve$saving_percent,

    type = "b",
    pch = 19,

    col = "#245A81",
    lwd = 2,

    ylim = limits + c(-padding, padding),

    xlab = "Change in selected active-current input (%)",

    ylab = "High-mode energy saving (%)",

    main = ""
  )

  title(
    "Sensitivity to active-current inputs",
    line = 4
  )

  lines(
    low_curve$change_percent,

    low_curve$saving_percent,

    type = "b",
    pch = 17,

    col = "#B24A3A",
    lwd = 2
  )

  abline(
    h = 0,
    lty = 2,
    col = "gray40"
  )

  abline(
    v = 0,
    lty = 3,
    col = "gray60"
  )

  points(
    0,
    baseline$saving_percent,

    pch = 23,
    bg = "#F0B323",
    cex = 1.4
  )

  legend(
    "top",

    inset = c(0, -0.09),

    xpd = NA,
    horiz = TRUE,
    cex = 0.8,
    bty = "n",

    legend = c(
      "High current only",
      "Low current only",
      "Baseline"
    ),

    col = c(
      "#245A81",
      "#B24A3A",
      "black"
    ),

    pch = c(19, 17, 23),

    pt.bg = "#F0B323",

    lty = c(1, 1, NA),

    bg = "white"
  )

  dev.off()

  # ------------------------------------------------------------
  # 13. Plot the tested joint-scenario ranges
  # ------------------------------------------------------------

  png(
    file.path(
      output_dir,
      "joint_current_scenarios.png"
    ),

    width = 1800,
    height = 1200,
    res = 200
  )

  par(
    mar = c(5, 5, 5.5, 2)
  )

  limits <- range(
    c(
      corners$saving_percent,
      baseline$saving_percent,
      0
    )
  )

  padding <- max(
    diff(limits) * 0.08,
    0.1
  )

  plot(
    c(
      0,
      max(joint_magnitudes_percent)
    ),

    limits,

    type = "n",

    ylim = limits + c(-padding, padding),

    xlab = "Magnitude of joint current perturbations (%)",

    ylab = "High-mode energy saving (%)",

    main = "Range among four tested joint scenarios"
  )

  segments(
    corner_summary$perturbation_magnitude_percent,

    corner_summary$minimum_saving_among_tested_corners_percent,

    corner_summary$perturbation_magnitude_percent,

    corner_summary$maximum_saving_among_tested_corners_percent,

    col = "#245A81",

    lwd = 4
  )

  points(
    corner_summary$perturbation_magnitude_percent,

    corner_summary$minimum_saving_among_tested_corners_percent,

    pch = 25,

    bg = "#B24A3A",

    cex = 1.3
  )

  points(
    corner_summary$perturbation_magnitude_percent,

    corner_summary$maximum_saving_among_tested_corners_percent,

    pch = 24,

    bg = "#245A81",

    cex = 1.3
  )

  points(
    0,

    baseline$saving_percent,

    pch = 23,

    bg = "#F0B323",

    cex = 1.4
  )

  abline(
    h = 0,

    lty = 2,

    col = "gray40"
  )

  legend(
    "topleft",

    c(
      "Maximum tested saving",
      "Minimum tested saving",
      "Baseline"
    ),

    pch = c(24, 25, 23),

    pt.bg = c(
      "#245A81",
      "#B24A3A",
      "#F0B323"
    ),

    bg = "white"
  )

  mtext(
    "Scenario ranges are not confidence intervals.",

    side = 3,

    line = 0.3,

    cex = 0.85
  )

  dev.off()

  # ------------------------------------------------------------
  # 14. Save interpretation notes
  # ------------------------------------------------------------

  notes <- c(
    paste(
      "Exploratory input robustness analysis;",
      "all runs use temperature feedback."
    ),

    paste(
      "Current perturbations are assumed stress scenarios,",
      "not paper-reported measurement errors."
    ),

    paste(
      "Voltage, execution times, task count, idle power",
      "and thermal assumptions remain fixed."
    ),

    paste(
      "Four joint corners are tested at each magnitude;",
      "no probability distribution is assigned."
    ),

    paste(
      "Minimum/maximum values describe tested corners,",
      "not proven extrema over a continuous box."
    ),

    paste(
      "Independent and opposing input changes are scenarios,",
      "not an established error correlation model."
    ),

    paste(
      "Root finding detects a bracketed crossing,",
      "not a proof of uniqueness."
    ),

    paste(
      "Time-step checks at crossings are numerical checks,",
      "not hardware validation."
    ),

    paste(
      "A ranking reversal establishes conditional sensitivity,",
      "not the probability of either strategy winning."
    ),

    "The original paper_data.csv is not modified."
  )

  writeLines(
    notes,

    file.path(
      output_dir,
      "analysis_notes.txt"
    )
  )

  # ------------------------------------------------------------
  # 15. Print the results
  # ------------------------------------------------------------

  cat(sprintf(
    "\nBaseline high-mode saving: %.6f%%\n",

    baseline$saving_percent
  ))

  cat(
    "\nSingle-input break-even changes:\n"
  )

  print(
    thresholds,
    row.names = FALSE
  )

  cat(
    "\nJoint-corner scenario summary:\n"
  )

  print(
    corner_summary,
    row.names = FALSE
  )

  cat(
    "Results written to outputs/input_robustness/\n"
  )

  invisible(
    list(
      oat = oat,
      corners = corners,
      thresholds = thresholds
    )
  )
}

run_input_robustness()

