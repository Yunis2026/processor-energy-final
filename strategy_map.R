# Joint idle-power/thermal-resistance scenarios and coupled break-even points.
# Run with: Rscript strategy_map.R
# Requires pde_model.R and paper_data.csv.
# No additional R packages are needed.

run_strategy_map <- function() {

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
      "Keep strategy_map.R, pde_model.R and paper_data.csv together."
    )
  }

  # ------------------------------------------------------------
  # 2. Load the existing PDE solver
  # ------------------------------------------------------------

  # Load initialization and solver definitions only.
  # Do not execute the original baseline analysis driver.
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

  if (!exists(
    "simulate_mode",
    envir = model,
    inherits = FALSE
  )) {
    stop("simulate_mode was not found in pde_model.R.")
  }

  if (exists(
    "scenario_results",
    envir = model,
    inherits = FALSE
  )) {
    stop(
      "The baseline driver ran unexpectedly; check the model structure."
    )
  }

  # ------------------------------------------------------------
  # 3. Create the output directory
  # ------------------------------------------------------------

  output_dir <- file.path(
    "outputs",
    "strategy_map"
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
  # 4. Define the joint scenarios
  # ------------------------------------------------------------

  # These ranges match the previous OAT study.
  # They are exploratory assumptions, not confidence intervals.
  idle_values_mW <- seq(
    0,
    0.30,
    length.out = 13
  )

  resistance_values <- seq(
    30,
    120,
    length.out = 7
  )

  baseline_idle_mW <-
    model$idle_reference_power_W * 1000

  baseline_resistance <-
    model$R_th_K_per_W

  fixed_alpha <- model$alpha_feedback
  fixed_beta <- model$beta_per_K

  nx <- 21L
  dt_s <- 1

  # This tolerance classifies numerical near-equality only.
  # It is not a claim about physical prediction accuracy.
  equality_tolerance_J <- 1e-6

  # ------------------------------------------------------------
  # 5. Save inputs and fixed settings
  # ------------------------------------------------------------

  write.csv(
    model$paper_data,
    file.path(output_dir, "input_snapshot.csv"),
    row.names = FALSE
  )

  settings <- data.frame(
    parameter = c(
      "alpha",
      "beta_per_K",
      "nx",
      "dt_s",
      "observation_time_s",
      "task_count",
      "L_m",
      "A_cs_m2",
      "volume_m3",
      "c_v_J_per_m3K",
      "k_W_per_mK",
      "T_ambient_K",
      "T_reference_K",
      "carbon_intensity_kg_CO2e_per_kWh",
      "equality_tolerance_J"
    ),

    value = c(
      fixed_alpha,
      fixed_beta,
      nx,
      dt_s,
      model$observation_time_s,
      model$task_count,
      model$L_m,
      model$A_cs_m2,
      model$volume_m3,
      model$c_v_J_per_m3K,
      model$k_W_per_mK,
      model$T_ambient_K,
      model$T_reference_K,
      model$carbon_intensity_kg_CO2e_per_kWh,
      equality_tolerance_J
    )
  )

  write.csv(
    settings,
    file.path(output_dir, "fixed_settings.csv"),
    row.names = FALSE
  )

  # ------------------------------------------------------------
  # 6. Define the paired simulation function
  # ------------------------------------------------------------

  cache <- new.env(parent = emptyenv())

  all_evaluations <- list()

  evaluate_pair <- function(
    idle_mW,
    resistance,
    time_step = dt_s
  ) {

    key <- paste(
      sprintf(
        "%.16g",
        c(idle_mW, resistance, time_step)
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

    # Update both the thermal resistance and its derived sink.
    model$idle_reference_power_W <- idle_mW / 1000

    model$R_th_K_per_W <- resistance

    model$h_v_W_per_m3K <-
      1 / (resistance * model$volume_m3)

    model$alpha_feedback <- fixed_alpha
    model$beta_per_K <- fixed_beta

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

      if (!is.finite(residual) || residual > 1e-7) {
        stop("Energy-balance residual exceeds 1e-7 J.")
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
        idle_power_mW = idle_mW,

        thermal_resistance_K_per_W =
          resistance,

        dt_s = time_step,

        status = "ok",
        error = "",

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

      data.frame(
        idle_power_mW = idle_mW,

        thermal_resistance_K_per_W =
          resistance,

        dt_s = time_step,

        status = "failed",
        error = conditionMessage(e),

        high_energy_J = NA_real_,
        low_energy_J = NA_real_,

        energy_difference_high_minus_low_J =
          NA_real_,

        saving_percent = NA_real_,

        preferred_strategy = NA_character_,

        high_peak_temperature_C = NA_real_,
        low_peak_temperature_C = NA_real_,

        high_carbon_mg_CO2e = NA_real_,
        low_carbon_mg_CO2e = NA_real_,

        max_energy_balance_residual_J =
          NA_real_,

        stringsAsFactors = FALSE
      )
    })

    all_evaluations[[
      length(all_evaluations) + 1L
    ]] <<- result

    # Preserve every unique evaluation, including root-finding trials.
    write.csv(
      do.call(rbind, all_evaluations),
      file.path(output_dir, "all_evaluations.csv"),
      row.names = FALSE
    )

    if (result$status != "ok") {
      stop(
        "A PDE calculation failed; inspect all_evaluations.csv: ",
        result$error
      )
    }

    cache[[key]] <- result

    result
  }

  # ------------------------------------------------------------
  # 7. Reproduce the coupled baseline
  # ------------------------------------------------------------

  cat("Computing the coupled baseline...\n")

  baseline <- evaluate_pair(
    baseline_idle_mW,
    baseline_resistance
  )

  write.csv(
    baseline,
    file.path(output_dir, "baseline.csv"),
    row.names = FALSE
  )

  # ------------------------------------------------------------
  # 8. Calculate the two-dimensional scenario grid
  # ------------------------------------------------------------

  map_rows <- list()

  for (resistance in resistance_values) {

    cat(sprintf(
      "Mapping thermal resistance = %.1f K/W...\n",
      resistance
    ))

    for (idle in idle_values_mW) {

      map_rows[[
        length(map_rows) + 1L
      ]] <- evaluate_pair(
        idle,
        resistance
      )
    }

    write.csv(
      do.call(rbind, map_rows),
      file.path(output_dir, "strategy_map.csv"),
      row.names = FALSE
    )
  }

  map_table <- do.call(rbind, map_rows)

  stopifnot(nrow(map_table) == 91L)

  # ------------------------------------------------------------
  # 9. Calculate the no-feedback analytical threshold
  # ------------------------------------------------------------

  active_times <-
    model$task_count *
    model$paper_data$task_time_s

  high_time <- active_times[model$high_row]

  low_time <- active_times[model$low_row]

  uncoupled_threshold_mW <- 1000 * (
    model$paper_data$power_W[model$low_row] *
      low_time -
      model$paper_data$power_W[model$high_row] *
      high_time
  ) / (low_time - high_time)

  # ------------------------------------------------------------
  # 10. Find the coupled break-even idle power
  # ------------------------------------------------------------

  root_rows <- list()

  for (resistance in resistance_values) {

    cat(sprintf(
      "Finding break-even idle power at %.1f K/W...\n",
      resistance
    ))

    energy_gap <- function(idle) {

      evaluate_pair(
        idle,
        resistance
      )$energy_difference_high_minus_low_J
    }

    lower <- min(idle_values_mW)
    upper <- max(idle_values_mW)

    f_lower <- energy_gap(lower)
    f_upper <- energy_gap(upper)

    # Search only when the endpoints bracket a crossing.
    if (f_lower == 0) {

      root <- lower

    } else if (f_upper == 0) {

      root <- upper

    } else if (
      sign(f_lower) != sign(f_upper)
    ) {

      root <- uniroot(
        energy_gap,

        interval = c(lower, upper),

        f.lower = f_lower,
        f.upper = f_upper,

        tol = 1e-7,
        maxiter = 100
      )$root

    } else {

      root <- NA_real_
    }

    if (is.finite(root)) {

      at_root <- evaluate_pair(
        root,
        resistance
      )

      root_gap <-
        at_root$energy_difference_high_minus_low_J

      if (
        abs(root_gap) >
        equality_tolerance_J
      ) {
        stop(
          "Root energy mismatch exceeds the requested tolerance."
        )
      }

    } else {

      root_gap <- NA_real_
    }

    root_rows[[
      length(root_rows) + 1L
    ]] <- data.frame(
      thermal_resistance_K_per_W =
        resistance,

      status = if (is.finite(root)) {
        "crossing_found"
      } else {
        "no_bracketed_crossing"
      },

      coupled_threshold_mW =
        root,

      no_feedback_threshold_mW =
        uncoupled_threshold_mW,

      threshold_change_mW =
        root - uncoupled_threshold_mW,

      root_energy_difference_J =
        root_gap
    )

    write.csv(
      do.call(rbind, root_rows),
      file.path(output_dir, "break_even.csv"),
      row.names = FALSE
    )
  }

  thresholds <- do.call(rbind, root_rows)

  # ------------------------------------------------------------
  # 11. Check the boundary points with a smaller time step
  # ------------------------------------------------------------

  # Re-evaluate the calculated roots using dt = 0.5 s.
  # A tiny gap near a root is not a robust physical preference.
  refinement_rows <- list()

  for (i in seq_len(nrow(thresholds))) {

    if (!is.finite(
      thresholds$coupled_threshold_mW[i]
    )) {
      next
    }

    cat(sprintf(
      "Checking boundary at %.1f K/W with dt = 0.5 s...\n",
      thresholds$thermal_resistance_K_per_W[i]
    ))

    refined <- evaluate_pair(
      thresholds$coupled_threshold_mW[i],

      thresholds$thermal_resistance_K_per_W[i],

      time_step = 0.5
    )

    refined$coarse_root_energy_difference_J <-
      thresholds$root_energy_difference_J[i]

    refined$gap_change_from_dt_1_J <-
      refined$energy_difference_high_minus_low_J -
      thresholds$root_energy_difference_J[i]

    refinement_rows[[
      length(refinement_rows) + 1L
    ]] <- refined
  }

  if (length(refinement_rows)) {

    write.csv(
      do.call(rbind, refinement_rows),

      file.path(
        output_dir,
        "boundary_time_step_check.csv"
      ),

      row.names = FALSE
    )
  }

  # ------------------------------------------------------------
  # 12. Prepare the strategy-map data
  # ------------------------------------------------------------

  # Idle power varies fastest within each thermal-resistance group.
  saving_matrix <- matrix(
    map_table$saving_percent,

    nrow = length(idle_values_mW),

    ncol = length(resistance_values)
  )

  color_limit <- max(abs(saving_matrix))

  if (color_limit == 0) {
    color_limit <- 1
  }

  colors <- colorRampPalette(
    c(
      "#B2182B",
      "white",
      "#2166AC"
    )
  )(101)

  finite_roots <- is.finite(
    thresholds$coupled_threshold_mW
  )

  # ------------------------------------------------------------
  # 13. Plot the two-dimensional strategy map
  # ------------------------------------------------------------

  png(
    file.path(
      output_dir,
      "strategy_map.png"
    ),

    width = 2000,
    height = 1300,
    res = 200
  )

  layout(
    matrix(c(1, 2), nrow = 1),
    widths = c(1, 0.15)
  )

  par(mar = c(5, 5, 5, 1))

  image(
    idle_values_mW,
    resistance_values,
    saving_matrix,

    col = colors,

    zlim = c(
      -color_limit,
      color_limit
    ),

    xlim = range(idle_values_mW),

    ylim = range(resistance_values),

    xlab = "Idle power (mW)",

    ylab = "Thermal resistance (K/W)",

    main = paste(
      "Operating-strategy comparison",
      "with temperature feedback"
    )
  )

  if (any(finite_roots)) {

    lines(
      thresholds$coupled_threshold_mW,

      thresholds$thermal_resistance_K_per_W,

      col = "black",
      lwd = 2
    )

    points(
      thresholds$coupled_threshold_mW,

      thresholds$thermal_resistance_K_per_W,

      pch = 16,
      cex = 0.7
    )
  }

  abline(
    v = uncoupled_threshold_mW,

    lty = 2,
    lwd = 2,
    col = "gray35"
  )

  points(
    baseline_idle_mW,
    baseline_resistance,

    pch = 23,
    bg = "#F0B323",
    cex = 1.4
  )

  legend(
    "topright",

    legend = c(
      "Coupled break-even",
      "No-feedback break-even",
      "Baseline"
    ),

    col = c(
      "black",
      "gray35",
      "black"
    ),

    lty = c(1, 2, NA),

    pch = c(NA, NA, 23),

    pt.bg = "#F0B323",

    bg = "white",
    cex = 0.8
  )

  mtext(
    paste(
      "Blue: high + idle saves energy;",
      "red: low saves energy."
    ),

    side = 3,
    line = 0.4,
    cex = 0.85
  )

  # Add a shared color bar.
  par(mar = c(5, 1, 5, 5))

  color_values <- seq(
    -color_limit,
    color_limit,

    length.out = length(colors)
  )

  image(
    c(0, 1),
    color_values,

    matrix(
      rep(color_values, each = 2),
      nrow = 2
    ),

    col = colors,

    axes = FALSE,
    xlab = "",
    ylab = ""
  )

  axis(4, las = 1)

  mtext(
    "High-mode saving (%)",

    side = 4,
    line = 3
  )

  dev.off()

  # ------------------------------------------------------------
  # 14. Plot the break-even threshold against thermal resistance
  # ------------------------------------------------------------

  png(
    file.path(
      output_dir,
      "break_even_vs_thermal_resistance.png"
    ),

    width = 1800,
    height = 1200,
    res = 200
  )

  par(mar = c(5, 5, 4, 2))

  limits <- range(
    c(
      thresholds$coupled_threshold_mW,
      uncoupled_threshold_mW
    ),

    na.rm = TRUE
  )

  padding <- max(
    diff(limits) * 0.15,
    0.001
  )

  plot(
    thresholds$thermal_resistance_K_per_W,

    thresholds$coupled_threshold_mW,

    type = "n",

    xlim = range(resistance_values),

    ylim = limits + c(-padding, padding),

    xlab = "Thermal resistance (K/W)",

    ylab = "Break-even idle power (mW)",

    main = paste(
      "Effect of thermal resistance",
      "on the strategy boundary"
    )
  )

  abline(
    h = uncoupled_threshold_mW,

    col = "gray35",
    lty = 2,
    lwd = 2
  )

  if (any(finite_roots)) {

    lines(
      thresholds$thermal_resistance_K_per_W,

      thresholds$coupled_threshold_mW,

      type = "b",
      pch = 19,

      col = "#245A81",
      lwd = 2
    )
  }

  legend(
    "bottomleft",

    c(
      "Coupled model",
      "No-feedback model"
    ),

    col = c(
      "#245A81",
      "gray35"
    ),

    lty = c(1, 2),

    pch = c(19, NA),

    lwd = 2,
    bty = "n"
  )

  dev.off()

  # ------------------------------------------------------------
  # 15. Save interpretation notes
  # ------------------------------------------------------------

  notes <- c(
    paste(
      "Joint idle-power and thermal-resistance study:",
      "13 x 7 = 91 scenarios."
    ),

    paste(
      "Alpha, beta, workload, geometry, material properties",
      "and carbon intensity are fixed."
    ),

    paste(
      "Ranges are exploratory assumptions,",
      "not calibrated uncertainty intervals."
    ),

    paste(
      "Positive saving favors high + idle;",
      "negative saving favors low frequency."
    ),

    paste(
      "Map colors represent sampled scenario results;",
      "lines connect calculated boundary points."
    ),

    paste(
      "Root finding detects a bracketed crossing,",
      "not a proof of uniqueness or monotonicity."
    ),

    paste(
      "Time-step checks at boundary points",
      "do not establish full spatial convergence."
    ),

    paste(
      "Near-equality within 1e-6 J is a numerical",
      "classification, not a physical accuracy claim."
    ),

    "Numerical consistency is not hardware validation.",

    paste(
      "A common fixed positive electricity emission factor",
      "preserves energy rankings."
    )
  )

  writeLines(
    notes,
    file.path(
      output_dir,
      "analysis_notes.txt"
    )
  )

  # ------------------------------------------------------------
  # 16. Print the final summary
  # ------------------------------------------------------------

  cat(sprintf(
    "\nBaseline high-mode saving: %.6f%%\n",
    baseline$saving_percent
  ))

  cat(sprintf(
    "No-feedback break-even idle power: %.6f mW\n",
    uncoupled_threshold_mW
  ))

  print(
    thresholds,
    row.names = FALSE
  )

  cat(
    "Results written to outputs/strategy_map/\n"
  )

  invisible(
    list(
      map = map_table,
      thresholds = thresholds
    )
  )
}

run_strategy_map()