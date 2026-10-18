# One-dimensional exploratory chip heat-transfer and use-phase emissions model.

paper_data <- read.csv("paper_data.csv", stringsAsFactors = FALSE)
paper_data$current_A <- paper_data$current_mA / 1000
paper_data$task_time_s <- paper_data$task_time_ms / 1000
paper_data$power_W <- paper_data$voltage_V * paper_data$current_A

high_row <- which(paper_data$mode == "high")
low_row <- which(paper_data$mode == "low")
if (length(high_row) != 1 || length(low_row) != 1) {
  stop("paper_data.csv must contain exactly one high row and one low row.")
}

observation_time_s <- 3600
task_count <- floor(observation_time_s / paper_data$task_time_s[low_row])
carbon_intensity_kg_CO2e_per_kWh <- 0.40
idle_reference_power_W <- 0.10 / 1000

L_m <- 0.003
A_cs_m2 <- 0.003 * 0.0003
volume_m3 <- L_m * A_cs_m2
c_v_J_per_m3K <- 1.6e6
k_W_per_mK <- 100
T_ambient_K <- 298.15
T_reference_K <- 298.15
R_th_K_per_W <- 60
h_v_W_per_m3K <- 1 / (R_th_K_per_W * volume_m3)
beta_per_K <- 0.05
alpha_feedback <- 0.20
temperature_tolerance_K <- 1e-8
max_picard_iterations <- 100L
abnormal_temperature_limit_K <- 1000
abnormal_step_change_limit_K <- 100

pde_output_dir <- file.path("outputs", "pde")
if (!dir.exists(pde_output_dir) && !dir.create(pde_output_dir, recursive = TRUE)) {
  stop("Could not create outputs/pde.")
}

make_grid <- function(nx) {
  if (length(nx) != 1 || nx < 2 || nx != as.integer(nx)) {
    stop("nx must be an integer of at least 2.")
  }
  dx <- L_m / nx
  list(
    nx = as.integer(nx),
    dx = dx,
    cell_volume = rep(A_cs_m2 * dx, as.integer(nx)),
    centers = (seq_len(nx) - 0.5) * dx
  )
}

make_reference_source <- function(grid, active_power_W, state, spatial_profile) {
  if (state == "idle") {
    return(rep(idle_reference_power_W / volume_m3, grid$nx))
  }
  if (spatial_profile == "uniform") {
    return(rep(active_power_W / volume_m3, grid$nx))
  }
  cell_left_m <- (seq_len(grid$nx) - 1) * grid$dx
  cell_right_m <- seq_len(grid$nx) * grid$dx
  overlap_m <- pmax(
    0,
    pmin(cell_right_m, 2 * L_m / 3) - pmax(cell_left_m, L_m / 3)
  )
  active_region_length_m <- L_m / 3
  source <- active_power_W / (A_cs_m2 * active_region_length_m) *
    overlap_m / grid$dx
  integrated_power_W <- sum(source * grid$cell_volume)
  if (abs(integrated_power_W - active_power_W) >
      1e-12 * active_power_W) {
    stop("Active reference source is not power-normalized on this grid.")
  }
  source
}

make_tridiagonal_solver <- function(grid, dt_s) {
  conductance_W_per_K <- k_W_per_mK * A_cs_m2 / grid$dx
  degree <- c(1, rep(2, grid$nx - 2), 1)
  lower <- rep(-dt_s * conductance_W_per_K, grid$nx - 1)
  upper <- lower
  diagonal <- c_v_J_per_m3K * grid$cell_volume +
    dt_s * (degree * conductance_W_per_K + h_v_W_per_m3K * grid$cell_volume)

  c_prime <- numeric(grid$nx - 1)
  denominator <- numeric(grid$nx)
  denominator[1] <- diagonal[1]
  if (grid$nx > 1) {
    c_prime[1] <- upper[1] / denominator[1]
  }
  if (grid$nx > 2) {
    for (i in 2:(grid$nx - 1)) {
      denominator[i] <- diagonal[i] - lower[i - 1] * c_prime[i - 1]
      c_prime[i] <- upper[i] / denominator[i]
    }
  }
  if (grid$nx > 1) {
    denominator[grid$nx] <-
      diagonal[grid$nx] - lower[grid$nx - 1] * c_prime[grid$nx - 1]
  }

  function(rhs) {
    d_prime <- numeric(grid$nx)
    d_prime[1] <- rhs[1] / denominator[1]
    if (grid$nx > 1) {
      for (i in 2:grid$nx) {
        d_prime[i] <-
          (rhs[i] - lower[i - 1] * d_prime[i - 1]) / denominator[i]
      }
    }
    solution <- numeric(grid$nx)
    solution[grid$nx] <- d_prime[grid$nx]
    if (grid$nx > 1) {
      for (i in (grid$nx - 1):1) {
        solution[i] <- d_prime[i] - c_prime[i] * solution[i + 1]
      }
    }
    solution
  }
}

conduction_power_W <- function(temperature_K, grid) {
  face_conductance_W_per_K <- k_W_per_mK * A_cs_m2 / grid$dx
  face_power_W <- face_conductance_W_per_K * diff(temperature_K)
  result <- numeric(grid$nx)
  result[-grid$nx] <- result[-grid$nx] + face_power_W
  result[-1] <- result[-1] - face_power_W
  result
}

simulate_mode <- function(mode, nx = 21L, dt_nominal_s = 1,
                          feedback = FALSE, record_history = FALSE,
                          spatial_profile = "center",
                          end_time_s = observation_time_s) {
  if (!is.finite(end_time_s) || end_time_s <= 0 ||
      end_time_s > observation_time_s) {
    stop("end_time_s must be greater than zero and no later than the observation time.")
  }
  row <- if (mode == "high") high_row else low_row
  active_power_W <- paper_data$power_W[row]
  task_time_s <- paper_data$task_time_s[row]
  active_time_s <- task_count * task_time_s
  if (active_time_s > observation_time_s) {
    stop(sprintf("%s mode active time exceeds the observation period.", mode))
  }

  grid <- make_grid(nx)
  temperature_K <- rep(T_ambient_K, grid$nx)
  time_s <- 0
  generated_energy_J <- 0
  heat_lost_J <- 0
  peak_temperature_K <- max(temperature_K)
  maximum_spatial_spread_K <- max(temperature_K) - min(temperature_K)
  history_times <- if (record_history) list(0) else NULL
  history_temperatures <- if (record_history) list(temperature_K) else NULL
  conduction_integral_max_W <- 0

  while (time_s < end_time_s) {
    remaining_s <- end_time_s - time_s
    if (remaining_s <= 1e-12) {
      time_s <- observation_time_s
      break
    }
    active_remaining_s <- active_time_s - time_s
    state <- if (active_remaining_s > 0) "active" else "idle"
    dt_s <- min(dt_nominal_s, remaining_s)
    if (state == "active") {
      dt_s <- min(dt_s, active_remaining_s)
    }
    if (!is.finite(dt_s) || dt_s <= 0) {
      stop(sprintf("Invalid time step encountered in %s mode.", mode))
    }

    q_reference <- make_reference_source(
      grid, active_power_W, state, spatial_profile
    )
    reference_power_W <- sum(q_reference * grid$cell_volume)
    expected_power_W <- if (state == "active") {
      active_power_W
    } else {
      idle_reference_power_W
    }
    if (abs(reference_power_W - expected_power_W) >
        1e-12 * max(1, expected_power_W)) {
      stop(sprintf("Reference heat source is not power-normalized in %s mode.", mode))
    }

    solve_step <- make_tridiagonal_solver(grid, dt_s)
    temperature_guess_K <- temperature_K
    if (feedback) {
      mode_feedback_W_per_m3 <-
        alpha_feedback * reference_power_W / volume_m3
      converged <- FALSE
      for (iteration in seq_len(max_picard_iterations)) {
        exponent <- beta_per_K * (temperature_guess_K - T_reference_K)
        feedback_source <- mode_feedback_W_per_m3 * (exp(exponent) - 1)
        q_generated <- q_reference + feedback_source
        if (any(!is.finite(q_generated))) {
          stop(sprintf("Non-finite feedback source in %s mode.", mode))
        }
        rhs <- c_v_J_per_m3K * grid$cell_volume *
          (temperature_K - T_ambient_K) +
          dt_s * q_generated * grid$cell_volume
        temperature_new_K <- T_ambient_K + solve_step(rhs)
        if (any(!is.finite(temperature_new_K))) {
          stop(sprintf("Non-finite temperature in %s mode.", mode))
        }
        if (max(temperature_new_K) > abnormal_temperature_limit_K ||
            max(abs(temperature_new_K - temperature_K)) >
            abnormal_step_change_limit_K) {
          stop(sprintf("Abnormal temperature growth in %s mode.", mode))
        }
        if (max(abs(temperature_new_K - temperature_guess_K)) <
            temperature_tolerance_K) {
          temperature_K_new <- temperature_new_K
          q_generated_step <- q_generated
          converged <- TRUE
          break
        }
        temperature_guess_K <- temperature_new_K
      }
      if (!converged) {
        stop(sprintf(
          "Picard iteration did not converge in %s mode at t = %.9f s.",
          mode, time_s
        ))
      }
      temperature_new_K <- temperature_K_new
    } else {
      rhs <- c_v_J_per_m3K * grid$cell_volume *
        (temperature_K - T_ambient_K) +
        dt_s * q_reference * grid$cell_volume
      temperature_new_K <- T_ambient_K + solve_step(rhs)
      q_generated_step <- q_reference
      if (any(!is.finite(temperature_new_K))) {
        stop(sprintf("Non-finite temperature in %s mode.", mode))
      }
      if (max(temperature_new_K) > abnormal_temperature_limit_K ||
          max(abs(temperature_new_K - temperature_K)) >
          abnormal_step_change_limit_K) {
        stop(sprintf("Abnormal temperature growth in %s mode.", mode))
      }
    }

    conduction <- conduction_power_W(temperature_new_K, grid)
    conduction_integral_W <- sum(conduction)
    conduction_integral_max_W <-
      max(conduction_integral_max_W, abs(conduction_integral_W))
    if (abs(conduction_integral_W) >
        1e-12 * max(1, sum(abs(conduction)))) {
      stop("Finite-volume conduction does not conserve energy.")
    }

    total_power_W <- sum(q_generated_step * grid$cell_volume)
    generated_energy_J <- generated_energy_J + total_power_W * dt_s
    heat_lost_J <- heat_lost_J +
      h_v_W_per_m3K * sum((temperature_new_K - T_ambient_K) *
                            grid$cell_volume) * dt_s

    temperature_K <- temperature_new_K
    peak_temperature_K <- max(peak_temperature_K, max(temperature_K))
    maximum_spatial_spread_K <- max(
      maximum_spatial_spread_K,
      max(temperature_K) - min(temperature_K)
    )
    time_s <- time_s + dt_s
    if (abs(time_s - active_time_s) <= 1e-12) {
      time_s <- active_time_s
    }
    if (abs(time_s - end_time_s) <= 1e-10) {
      time_s <- end_time_s
    }
    if (record_history) {
      history_times[[length(history_times) + 1L]] <- time_s
      history_temperatures[[length(history_temperatures) + 1L]] <-
        temperature_K
    }
  }

  stored_heat_change_J <- c_v_J_per_m3K *
    sum((temperature_K - T_ambient_K) * grid$cell_volume)
  energy_balance_residual_J <-
    generated_energy_J - heat_lost_J - stored_heat_change_J
  result <- list(
    mode = mode,
    feedback = feedback,
    task_count = task_count,
    active_time_s = active_time_s,
    total_energy_J = generated_energy_J,
    carbon_mg_CO2e =
      generated_energy_J / 3600000 *
      carbon_intensity_kg_CO2e_per_kWh * 1e6,
    peak_temperature_C = peak_temperature_K - 273.15,
    final_mean_temperature_C = mean(temperature_K) - 273.15,
    energy_balance_residual_J = energy_balance_residual_J,
    conduction_integral_max_W = conduction_integral_max_W,
    maximum_spatial_spread_K = maximum_spatial_spread_K,
    temperature_final_K = temperature_K,
    grid = grid
  )
  if (record_history) {
    result$times_s <- unlist(history_times, use.names = FALSE)
    result$temperature_history_K <-
      do.call(cbind, history_temperatures)
  }
  result
}

mode_order <- c("high", "low")
scenario_results <- list()
for (mode in mode_order) {
  scenario_results[[paste(mode, "no_feedback", sep = "_")]] <-
    simulate_mode(mode, feedback = FALSE)
  scenario_results[[paste(mode, "feedback", sep = "_")]] <-
    simulate_mode(mode, feedback = TRUE, record_history = TRUE)
}

summary_table <- do.call(rbind, lapply(scenario_results, function(result) {
  data.frame(
    mode = result$mode,
    feedback = result$feedback,
    task_count = result$task_count,
    active_time_s = result$active_time_s,
    total_energy_J = result$total_energy_J,
    carbon_mg_CO2e = result$carbon_mg_CO2e,
    peak_temperature_C = result$peak_temperature_C,
    final_mean_temperature_C = result$final_mean_temperature_C,
    energy_balance_residual_J = result$energy_balance_residual_J
  )
}))
row.names(summary_table) <- NULL
write.csv(summary_table, file.path(pde_output_dir, "summary.csv"), row.names = FALSE)

parameters <- data.frame(
  parameter = c(
    "observation_time", "task_count", "idle_reference_power", "carbon_intensity",
    "workload",
    "chip_length", "cross_section_area", "chip_volume", "volumetric_heat_capacity",
    "thermal_conductivity", "ambient_temperature", "reference_temperature",
    "thermal_resistance", "volumetric_heat_transfer_coefficient",
    "feedback_beta", "feedback_alpha", "active_source_region",
    "idle_source_distribution", "active_source_distribution",
    "default_spatial_cells", "default_time_step", "temperature_tolerance",
    "maximum_picard_iterations", "abnormal_temperature_limit",
    "abnormal_step_change_limit", "high_reference_power", "high_task_time",
    "low_reference_power", "low_task_time", "high_active_duration",
    "low_active_duration", "model_boundary"
  ),
  value = c(
    observation_time_s, task_count, idle_reference_power_W,
    carbon_intensity_kg_CO2e_per_kWh, "SHA-256", L_m, A_cs_m2, volume_m3,
    c_v_J_per_m3K, k_W_per_mK, T_ambient_K, T_reference_K, R_th_K_per_W,
    h_v_W_per_m3K, beta_per_K, alpha_feedback,
    "Middle third of the chip length",
    "Uniform over the chip volume",
    "Uniform power density over selected middle-third cells",
    "21", "1", temperature_tolerance_K, max_picard_iterations,
    abnormal_temperature_limit_K, abnormal_step_change_limit_K,
    paper_data$power_W[high_row], paper_data$task_time_s[high_row],
    paper_data$power_W[low_row], paper_data$task_time_s[low_row],
    task_count * paper_data$task_time_s[high_row],
    task_count * paper_data$task_time_s[low_row],
    "Use-phase processor electricity only; no conversion losses, manufacturing, transport, replacement, or end-of-life"
  ),
  unit = c(
    "s", "tasks", "W", "kg CO2e/kWh", "description", "m", "m^2", "m^3",
    "J/(m^3 K)", "W/(m K)", "K", "K", "K/W", "W/(m^3 K)",
    "1/K", "1", "description", "description", "description", "cells",
    "s", "K", "iterations", "K", "K/step", "W", "s", "W", "s", "s",
    "s", "description"
  ),
  stringsAsFactors = FALSE
)
parameters$status <- "assumed"
parameters$status[parameters$parameter %in%
                    c("task_count", "high_active_duration", "low_active_duration")] <-
  "derived"
parameters$status[parameters$parameter %in% c(
  "high_reference_power", "high_task_time", "low_reference_power", "low_task_time"
)] <- "paper_data.csv"
write.csv(parameters, file.path(pde_output_dir, "parameters.csv"), row.names = FALSE)

high_feedback <- scenario_results[["high_feedback"]]
low_feedback <- scenario_results[["low_feedback"]]
short_high_feedback <- simulate_mode(
  "high", nx = 41L, dt_nominal_s = 0.025, feedback = TRUE,
  record_history = TRUE, end_time_s = 5
)
short_low_feedback <- simulate_mode(
  "low", nx = 41L, dt_nominal_s = 0.025, feedback = TRUE,
  record_history = TRUE, end_time_s = 5
)
if (short_high_feedback$times_s[1] != 0 ||
    tail(short_high_feedback$times_s, 1) != 5 ||
    short_low_feedback$times_s[1] != 0 ||
    tail(short_low_feedback$times_s, 1) != 5) {
  stop("Short-time simulations must include t = 0 and t = 5 s.")
}
history_min_C <- min(
  high_feedback$temperature_history_K,
  low_feedback$temperature_history_K
) - 273.15
history_max_C <- max(
  high_feedback$temperature_history_K,
  low_feedback$temperature_history_K
) - 273.15
plot_ylim_C <- c(history_min_C, history_max_C)
if (diff(plot_ylim_C) == 0) {
  plot_ylim_C <- plot_ylim_C + c(-1, 1) * 1e-8
}

png(file.path(pde_output_dir, "temperature_history.png"),
    width = 1600, height = 1000, res = 200)
plot(
  high_feedback$times_s,
  apply(high_feedback$temperature_history_K, 2, max) - 273.15,
  type = "l",
  col = "red",
  xlab = "Time (s)",
  ylab = "Peak temperature (degree C)",
  ylim = plot_ylim_C,
  main = "Peak temperature with temperature feedback"
)
lines(
  low_feedback$times_s,
  apply(low_feedback$temperature_history_K, 2, max) - 273.15,
  col = "blue"
)
legend("topright", c("16 MHz", "2.097 MHz"), col = c("red", "blue"),
       lty = 1, bg = "white")
dev.off()

high_switch_s <- high_feedback$active_time_s
early_windows <- list(c(0, 30), c(max(0, high_switch_s - 15), high_switch_s + 15))
png(file.path(pde_output_dir, "temperature_early.png"),
    width = 1800, height = 1000, res = 200)
par(mfrow = c(1, 2))
for (window in early_windows) {
  high_keep <- high_feedback$times_s >= window[1] &
    high_feedback$times_s <= window[2]
  low_keep <- low_feedback$times_s >= window[1] &
    low_feedback$times_s <= window[2]
  plot(
    high_feedback$times_s[high_keep],
    apply(high_feedback$temperature_history_K[, high_keep, drop = FALSE], 2, max) -
      273.15,
    type = "l",
    col = "red",
    xlab = "Time (s)",
    ylab = "Peak temperature (degree C)",
    ylim = plot_ylim_C,
    main = sprintf("%.3f to %.3f s", window[1], window[2])
  )
  lines(
    low_feedback$times_s[low_keep],
    apply(low_feedback$temperature_history_K[, low_keep, drop = FALSE], 2, max) -
      273.15,
    col = "blue"
  )
}
legend("top", inset = -0.18, xpd = TRUE,
       c("16 MHz", "2.097 MHz"), col = c("red", "blue"), lty = 1,
       horiz = TRUE, bty = "n")
par(mfrow = c(1, 1))
dev.off()

heatmap_palette <- colorRampPalette(c("navy", "cyan", "yellow", "red"))(100)
heatmap_limits_C <- range(
  high_feedback$temperature_history_K,
  low_feedback$temperature_history_K
) - 273.15
png(file.path(pde_output_dir, "temperature_heatmap.png"),
    width = 1800, height = 850, res = 200)
layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 0.22))
par(mar = c(5, 5, 4, 2))
image(
  high_feedback$times_s, high_feedback$grid$centers * 1000,
  t(high_feedback$temperature_history_K) - 273.15,
  col = heatmap_palette, zlim = heatmap_limits_C,
  xlab = "Time (s)", ylab = "Position (mm)", main = "16 MHz"
)
image(
  low_feedback$times_s, low_feedback$grid$centers * 1000,
  t(low_feedback$temperature_history_K) - 273.15,
  col = heatmap_palette, zlim = heatmap_limits_C,
  xlab = "Time (s)", ylab = "Position (mm)", main = "2.097 MHz"
)
par(mar = c(5, 1, 4, 5))
colorbar_values_C <- seq(heatmap_limits_C[1], heatmap_limits_C[2], length.out = 100)
image(
  c(0, 1), colorbar_values_C,
  matrix(rep(colorbar_values_C, each = 2), nrow = 2),
  col = heatmap_palette, axes = FALSE, xlab = "", ylab = ""
)
axis(4, at = pretty(heatmap_limits_C), las = 1, cex.axis = 0.8)
mtext("Temperature (°C)", side = 4, line = 3, cex = 0.8)
dev.off()

short_temperature_range_C <- range(
  short_high_feedback$temperature_history_K,
  short_low_feedback$temperature_history_K
) - 273.15
early_keep_high <- short_high_feedback$times_s >= 0 &
  short_high_feedback$times_s <= 5
early_keep_low <- short_low_feedback$times_s >= 0 &
  short_low_feedback$times_s <= 5
png(file.path(pde_output_dir, "temperature_first_5s.png"),
    width = 1800, height = 850, res = 200)
layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 0.22))
par(mar = c(5, 5, 4, 2))
image(
  short_high_feedback$times_s[early_keep_high],
  short_high_feedback$grid$centers * 1000,
  t(short_high_feedback$temperature_history_K[, early_keep_high, drop = FALSE]) -
    273.15,
  col = heatmap_palette, zlim = short_temperature_range_C,
  xlim = c(0, 5),
  xlab = "Time (s)", ylab = "Position (mm)", main = "16 MHz"
)
axis(1, at = seq(0, 5, by = 1))
image(
  short_low_feedback$times_s[early_keep_low],
  short_low_feedback$grid$centers * 1000,
  t(short_low_feedback$temperature_history_K[, early_keep_low, drop = FALSE]) -
    273.15,
  col = heatmap_palette, zlim = short_temperature_range_C,
  xlim = c(0, 5),
  xlab = "Time (s)", ylab = "Position (mm)", main = "2.097 MHz"
)
axis(1, at = seq(0, 5, by = 1))
par(mar = c(5, 1, 4, 5))
image(
  c(0, 1), seq(short_temperature_range_C[1], short_temperature_range_C[2],
                length.out = 100),
  matrix(rep(seq(short_temperature_range_C[1], short_temperature_range_C[2],
                 length.out = 100), each = 2), nrow = 2),
  col = heatmap_palette, axes = FALSE, xlab = "", ylab = ""
)
axis(4, at = pretty(short_temperature_range_C), las = 1, cex.axis = 0.8)
mtext("Temperature (°C)", side = 4, line = 3, cex = 0.8)
dev.off()

high_t5_index <- which.min(abs(short_high_feedback$times_s - 5))
low_t5_index <- which.min(abs(short_low_feedback$times_s - 5))
high_t5_profile_C <-
  short_high_feedback$temperature_history_K[, high_t5_index] - 273.15
low_t5_profile_C <-
  short_low_feedback$temperature_history_K[, low_t5_index] - 273.15
shared_t5_limits_C <- range(high_t5_profile_C, low_t5_profile_C)
shared_t5_padding_C <- max(diff(shared_t5_limits_C) * 0.05, 1e-7)
shared_t5_limits_C <- shared_t5_limits_C + c(-1, 1) * shared_t5_padding_C
png(file.path(pde_output_dir, "temperature_t5_profile.png"),
    width = 1400, height = 900, res = 200)
plot(
  short_high_feedback$grid$centers * 1000,
  high_t5_profile_C,
  type = "l", col = "red", lwd = 2,
  xlab = "Position (mm)", ylab = "Temperature (°C)",
  ylim = shared_t5_limits_C,
  main = "Spatial temperature profile at t = 5 s"
)
lines(
  short_low_feedback$grid$centers * 1000,
  low_t5_profile_C,
  col = "blue", lwd = 2
)
legend("topright", c("16 MHz", "2.097 MHz"),
       col = c("red", "blue"), lty = 1, lwd = 2, bg = "white")
dev.off()

png(file.path(pde_output_dir, "temperature_t5_profile_closeup.png"),
    width = 1800, height = 900, res = 200)
par(mfrow = c(1, 2))
for (profile in list(
  list(mode = "16 MHz", x = short_high_feedback$grid$centers * 1000,
       temperature_C = high_t5_profile_C),
  list(mode = "2.097 MHz", x = short_low_feedback$grid$centers * 1000,
       temperature_C = low_t5_profile_C)
)) {
  profile_range_C <- range(profile$temperature_C)
  profile_padding_C <- max(diff(profile_range_C) * 0.05, 1e-7)
  profile_limits_C <- profile_range_C + c(-1, 1) * profile_padding_C
  plot(
    profile$x, profile$temperature_C,
    type = "l", lwd = 2, col = if (profile$mode == "16 MHz") "red" else "blue",
    xlab = "Position (mm)", ylab = "Temperature (°C)",
    ylim = profile_limits_C,
    main = sprintf("%s: %.8f to %.8f °C",
                   profile$mode, profile_limits_C[1], profile_limits_C[2])
  )
}
par(mfrow = c(1, 1))
dev.off()

energy_matrix <- rbind(
  c(summary_table$total_energy_J[summary_table$mode == "high" &
                                   !summary_table$feedback],
    summary_table$total_energy_J[summary_table$mode == "low" &
                                   !summary_table$feedback]),
  c(summary_table$total_energy_J[summary_table$mode == "high" &
                                   summary_table$feedback],
    summary_table$total_energy_J[summary_table$mode == "low" &
                                   summary_table$feedback])
)
rownames(energy_matrix) <- c("No feedback", "Feedback")
colnames(energy_matrix) <- c("16 MHz", "2.097 MHz")
png(file.path(pde_output_dir, "energy_comparison.png"),
    width = 1400, height = 1000, res = 200)
par(mar = c(5, 5, 7, 2))
bar_positions <- barplot(
  energy_matrix,
  beside = TRUE,
  col = c("gray55", "darkorange"),
  names.arg = colnames(energy_matrix),
  xlab = "Mode",
  ylab = "Total energy (J)",
  main = "Energy by mode and feedback scenario",
  ylim = c(0, max(energy_matrix) * 1.18)
)
text(
  bar_positions,
  as.vector(energy_matrix),
  labels = sprintf("%.4f", as.vector(energy_matrix)),
  pos = 3, cex = 0.8
)
legend(
  "top", inset = c(0, -0.12), xpd = TRUE, horiz = TRUE,
  legend = rownames(energy_matrix),
  fill = c("gray55", "darkorange"), bty = "n"
)
dev.off()

uniform_validation <- list()
analytic_uniform_end <- function(start_K, power_W, duration_s) {
  tau_s <- R_th_K_per_W * c_v_J_per_m3K * volume_m3
  T_ambient_K + R_th_K_per_W * power_W +
    (start_K - T_ambient_K - R_th_K_per_W * power_W) *
    exp(-duration_s / tau_s)
}

simulate_uniform_transient <- function(start_temperature_K, power_W,
                                       duration_s, dt_s) {
  total_heat_capacity_J_per_K <- c_v_J_per_m3K * volume_m3
  temperature_K <- start_temperature_K
  elapsed_s <- 0
  max_error_K <- 0
  while (elapsed_s < duration_s) {
    actual_dt_s <- min(dt_s, duration_s - elapsed_s)
    temperature_K <- T_ambient_K +
      (total_heat_capacity_J_per_K / actual_dt_s *
         (temperature_K - T_ambient_K) + power_W) /
      (total_heat_capacity_J_per_K / actual_dt_s + 1 / R_th_K_per_W)
    elapsed_s <- elapsed_s + actual_dt_s
    exact_temperature_K <-
      analytic_uniform_end(start_temperature_K, power_W, elapsed_s)
    max_error_K <- max(max_error_K, abs(temperature_K - exact_temperature_K))
    if (!is.finite(temperature_K) || !is.finite(max_error_K)) {
      stop("Non-finite value in uniform transient validation.")
    }
  }
  max_error_K
}

transient_validation_rows <- list()
transient_time_steps_s <- c(0.2, 0.1, 0.05, 0.025)
for (mode in mode_order) {
  row <- if (mode == "high") high_row else low_row
  startup_errors_K <- vapply(
    transient_time_steps_s,
    function(dt_s) simulate_uniform_transient(
      T_ambient_K, paper_data$power_W[row], 5, dt_s
    ),
    numeric(1)
  )
  startup_monotonic <- all(diff(startup_errors_K) < 0)
  transient_validation_rows[[length(transient_validation_rows) + 1L]] <-
    data.frame(
      phase = "startup_first_5_seconds",
      mode = mode,
      dt_s = transient_time_steps_s,
      max_absolute_temperature_error_K = startup_errors_K,
      error_decreases_with_refinement = startup_monotonic
    )
}
high_active_duration_s <- task_count * paper_data$task_time_s[high_row]
high_analytic_switch_K <- analytic_uniform_end(
  T_ambient_K, paper_data$power_W[high_row], high_active_duration_s
)
cooldown_errors_K <- vapply(
  transient_time_steps_s,
  function(dt_s) simulate_uniform_transient(
    high_analytic_switch_K, idle_reference_power_W, 5, dt_s
  ),
  numeric(1)
)
cooldown_monotonic <- all(diff(cooldown_errors_K) < 0)
transient_validation_rows[[length(transient_validation_rows) + 1L]] <-
  data.frame(
    phase = "high_idle_first_5_seconds",
    mode = "high",
    dt_s = transient_time_steps_s,
    max_absolute_temperature_error_K = cooldown_errors_K,
    error_decreases_with_refinement = cooldown_monotonic
  )
transient_validation <- do.call(rbind, transient_validation_rows)
row.names(transient_validation) <- NULL
if (!all(transient_validation$error_decreases_with_refinement)) {
  stop("Transient validation error did not decrease with time-step refinement.")
}
write.csv(
  transient_validation,
  file.path(pde_output_dir, "transient_validation.csv"),
  row.names = FALSE
)

for (mode in mode_order) {
  uniform_result <- simulate_mode(
    mode, nx = 21L, dt_nominal_s = 1, feedback = FALSE,
    spatial_profile = "uniform"
  )
  row <- if (mode == "high") high_row else low_row
  active_duration_s <- task_count * paper_data$task_time_s[row]
  analytic_active_K <- analytic_uniform_end(
    T_ambient_K, paper_data$power_W[row], active_duration_s
  )
  analytic_end_K <- analytic_uniform_end(
    analytic_active_K, idle_reference_power_W,
    observation_time_s - active_duration_s
  )
  uniform_spread_K <- uniform_result$maximum_spatial_spread_K
  analytic_error_K <- mean(uniform_result$temperature_final_K) - analytic_end_K
  if (uniform_spread_K > 1e-9) {
    stop(sprintf("Uniform-source spatial invariance failed for %s mode.", mode))
  }
  uniform_validation[[mode]] <- data.frame(
    study = "uniform_analytic_validation",
    mode = mode,
    nx = 21L,
    dt_s = 1,
    peak_temperature_C = uniform_result$peak_temperature_C,
    total_energy_J = uniform_result$total_energy_J,
    peak_temperature_difference = NA_real_,
    energy_difference_J = NA_real_,
    uniform_spatial_spread_K = uniform_spread_K,
    analytic_final_temperature_error_K = analytic_error_K,
    energy_balance_residual_J = uniform_result$energy_balance_residual_J
  )
}

convergence_runs <- list()
for (mode in mode_order) {
  for (dt_value in c(2, 1, 0.5)) {
    run <- simulate_mode(mode, nx = 21L, dt_nominal_s = dt_value, feedback = TRUE)
    convergence_runs[[paste(mode, "time", dt_value, sep = "_")]] <- run
  }
  for (nx_value in c(21L, 41L, 81L)) {
    run <- simulate_mode(mode, nx = nx_value, dt_nominal_s = 0.5, feedback = TRUE)
    convergence_runs[[paste(mode, "space", nx_value, sep = "_")]] <- run
  }
}

convergence_rows <- list()
for (mode in mode_order) {
  time_reference <- convergence_runs[[paste(mode, "time", 0.5, sep = "_")]]
  for (dt_value in c(2, 1, 0.5)) {
    run <- convergence_runs[[paste(mode, "time", dt_value, sep = "_")]]
    convergence_rows[[length(convergence_rows) + 1L]] <- data.frame(
      study = "time_step",
      mode = mode,
      nx = 21L,
      dt_s = dt_value,
      peak_temperature_C = run$peak_temperature_C,
      total_energy_J = run$total_energy_J,
      peak_temperature_difference =
        run$peak_temperature_C - time_reference$peak_temperature_C,
      energy_difference_J = run$total_energy_J - time_reference$total_energy_J,
      uniform_spatial_spread_K = NA_real_,
      analytic_final_temperature_error_K = NA_real_,
      energy_balance_residual_J = run$energy_balance_residual_J
    )
  }
  space_reference <- convergence_runs[[paste(mode, "space", 81L, sep = "_")]]
  for (nx_value in c(21L, 41L, 81L)) {
    run <- convergence_runs[[paste(mode, "space", nx_value, sep = "_")]]
    convergence_rows[[length(convergence_rows) + 1L]] <- data.frame(
      study = "spatial_grid",
      mode = mode,
      nx = nx_value,
      dt_s = 0.5,
      peak_temperature_C = run$peak_temperature_C,
      total_energy_J = run$total_energy_J,
      peak_temperature_difference =
        run$peak_temperature_C - space_reference$peak_temperature_C,
      energy_difference_J = run$total_energy_J - space_reference$total_energy_J,
      uniform_spatial_spread_K = NA_real_,
      analytic_final_temperature_error_K = NA_real_,
      energy_balance_residual_J = run$energy_balance_residual_J
    )
  }
}
convergence_table <- do.call(rbind, c(convergence_rows, uniform_validation))
row.names(convergence_table) <- NULL
write.csv(
  convergence_table,
  file.path(pde_output_dir, "convergence.csv"),
  row.names = FALSE
)

baseline_energy_J <- c(
  high = paper_data$power_W[high_row] *
    task_count * paper_data$task_time_s[high_row] +
    idle_reference_power_W *
    (observation_time_s - task_count * paper_data$task_time_s[high_row]),
  low = paper_data$power_W[low_row] *
    task_count * paper_data$task_time_s[low_row] +
    idle_reference_power_W *
    (observation_time_s - task_count * paper_data$task_time_s[low_row])
)
no_feedback_energy <- summary_table$total_energy_J[!summary_table$feedback]
no_feedback_modes <- summary_table$mode[!summary_table$feedback]
baseline_difference_J <- no_feedback_energy - baseline_energy_J[no_feedback_modes]
if (max(abs(baseline_difference_J)) > 1e-8) {
  stop("No-feedback energy does not match the existing baseline model.")
}

feedback_high_energy <- summary_table$total_energy_J[
  summary_table$mode == "high" & summary_table$feedback
]
feedback_low_energy <- summary_table$total_energy_J[
  summary_table$mode == "low" & summary_table$feedback
]
feedback_ranking <- if (feedback_high_energy < feedback_low_energy) {
  "high"
} else if (feedback_high_energy > feedback_low_energy) {
  "low"
} else {
  "equal"
}
no_feedback_ranking <- if (baseline_energy_J["high"] < baseline_energy_J["low"]) {
  "high"
} else {
  "low"
}

time_convergence <- convergence_table[
  convergence_table$study == "time_step" & convergence_table$dt_s %in% c(1, 0.5),
]
space_convergence <- convergence_table[
  convergence_table$study == "spatial_grid" &
    convergence_table$nx %in% c(41, 81),
]

note_lines <- c(
  "# 一维芯片热传导与使用阶段排放模型说明",
  "",
  "## 物理模型与边界",
  "",
  "模型把芯片表示为长度方向上的一维均质介质。温度变化由热传导、体积发热和向环境的等效体积散热共同决定：",
  "",
  "`c_v * ∂T/∂t = k * ∂²T/∂x² + q_gen - h_v * (T - T_a)`",
  "",
  "初始条件为 `T(x, 0) = T_a`；两端采用绝热 Neumann 边界 `∂T/∂x = 0`。体积散热项是集中热阻 `R_th` 的分布式等效表示，满足 `h_v = 1 / (R_th * volume)`。",
  "",
  "## 数据来源与假设",
  "",
  "`paper_data.csv` 提供 high/low 模式的供电电压、电流和单任务时间；参考功率由电压与电流计算。任务数为 `floor(3600 / low_task_time_s)`，两种模式完成相同任务数，工作完成后 idle 至 3600 s。idle 参考功率 0.10 mW、排放因子 0.40 kg CO2e/kWh 均为情景假设。",
  "",
  "几何、热容、导热率、热阻、环境温度、反馈系数均为未实测校准的假设参数。active 热源固定在 [L/3, 2L/3]，按各单元与该物理区域的重叠长度分配并归一化；idle 热源均匀分布。空间发热分布是假设，不是芯片版图实测。",
  "",
  "温度反馈为 `q_gen = q_ref + a_mode * (exp(beta * (T - T_ref)) - 1)`，其中每个时间步的 `a_mode = alpha * reference_power_W / volume`，`reference_power_W` 是当前 operating state（active 或 idle）的参考功率体积积分。反馈只表示相对参考温度的增量，不重复添加完整漏电功耗。`alpha` 与 `beta` 未经数据拟合。",
  "",
  "## 能量与排放边界",
  "",
  "每步将总发热功率积分为用电能量，并以固定排放因子计算间接使用阶段排放。仅包括处理器用电对应的间接排放；忽略电源转换损耗，不包括制造、运输、更换和报废排放。",
  "",
  "## 数值方法",
  "",
  "采用守恒有限体积离散，内部相邻单元通过面热导交换热量，两端无热通量。后向 Euler 处理时间；无反馈时解线性三对角系统，有反馈时用 Picard 迭代，温度容差 1e-8 K、最多 100 次。时间步在 active/idle 切换时及终止时刻精确缩短。程序逐步检查有限体积导热项积分为零，并以生成能量减去散热和储热变化计算能量平衡残差。",
  "",
  "## 验证结果",
  "",
  sprintf("无反馈总能量与基础模型的最大差异为 %.6g J（要求不超过 1e-8 J）。", max(abs(baseline_difference_J))),
  sprintf("均匀热源额外验证的最大空间温差为 %.6g K；与分段整体热模型解析终温相比，最大绝对误差为 %.6g K。", max(vapply(uniform_validation, function(x) x$uniform_spatial_spread_K, numeric(1))), max(abs(vapply(uniform_validation, function(x) x$analytic_final_temperature_error_K, numeric(1))))),
  sprintf("短时瞬态验证前 5 秒采用 dt = 0.2、0.1、0.05、0.025 s；最大解析温度误差为 %.6g K，误差随步长减小而下降：%s。high 切换后的前 5 秒从解析切换温度开始，最大误差为 %.6g K，随步长减小而下降：%s。详见 transient_validation.csv。", max(transient_validation$max_absolute_temperature_error_K[transient_validation$phase == "startup_first_5_seconds"]), if (all(vapply(split(transient_validation$error_decreases_with_refinement[transient_validation$phase == "startup_first_5_seconds"], transient_validation$mode[transient_validation$phase == "startup_first_5_seconds"]), all, logical(1)))) "是" else "否", max(transient_validation$max_absolute_temperature_error_K[transient_validation$phase == "high_idle_first_5_seconds"]), if (cooldown_monotonic) "是" else "否"),
  sprintf("主情景反馈模拟最大绝对能量平衡残差为 %.6g J。", max(abs(summary_table$energy_balance_residual_J))),
  sprintf("时间步从 1 s 减至 0.5 s：%s；峰值温度变化 %.6g K，总能量变化 %.6g J。", if (abs(time_convergence$peak_temperature_difference[1]) < 1e-3 && abs(time_convergence$energy_difference_J[1]) < 1e-6) "数值变化较小" else "仍有可见数值变化，需结合数据谨慎判断", time_convergence$peak_temperature_difference[1], time_convergence$energy_difference_J[1]),
  sprintf("空间网格从 41 增至 81 单元：%s；峰值温度变化 %.6g K，总能量变化 %.6g J。", if (abs(space_convergence$peak_temperature_difference[1]) < 1e-3 && abs(space_convergence$energy_difference_J[1]) < 1e-6) "数值变化较小" else "仍有可见数值变化，需结合数据谨慎判断", space_convergence$peak_temperature_difference[1], space_convergence$energy_difference_J[1]),
  sprintf("未反馈时能量较低的策略为 %s；反馈后为 %s。", no_feedback_ranking, feedback_ranking),
  "固定电网排放因子为正且对两种模式相同，因此碳排放排序与能量排序一致。",
  "时间和空间网格的细化差值列于 convergence.csv。空间源项按固定区域与网格单元的重叠长度分配；峰值温度的网格变化仍可能非单调。这里报告数值敏感性，不把运行成功等同于严格的收敛证明。",
  "",
  "这些结果是探索性假设情景，不是经过实测校准的 STM32L0 数字孪生。温升和策略差异应按模型输出如实解读，不对温度变化作人为放大。",
  ""
)
writeLines(note_lines, file.path(pde_output_dir, "model_notes.md"), useBytes = TRUE)

cat("PDE summary:\n")
print(summary_table)
cat(sprintf("No-feedback baseline energy maximum difference: %.9g J\n",
            max(abs(baseline_difference_J))))
cat(sprintf("Maximum absolute energy-balance residual: %.9g J\n",
            max(abs(summary_table$energy_balance_residual_J))))
cat(sprintf("Uniform validation maximum spatial spread: %.9g K\n",
            max(vapply(uniform_validation, function(x) x$uniform_spatial_spread_K, numeric(1)))))
cat(sprintf("Uniform validation maximum analytic final-temperature error: %.9g K\n",
            max(abs(vapply(uniform_validation, function(x) x$analytic_final_temperature_error_K, numeric(1))))))
cat(sprintf("Startup transient maximum analytic temperature error: %.9g K\n",
            max(transient_validation$max_absolute_temperature_error_K[
              transient_validation$phase == "startup_first_5_seconds"
            ])))
cat(sprintf("High-mode cooldown maximum analytic temperature error: %.9g K\n",
            max(transient_validation$max_absolute_temperature_error_K[
              transient_validation$phase == "high_idle_first_5_seconds"
            ])))
cat(sprintf("Transient error decreases with refinement: startup %s, cooldown %s\n",
            if (all(vapply(split(
              transient_validation$error_decreases_with_refinement[
                transient_validation$phase == "startup_first_5_seconds"
              ],
              transient_validation$mode[
                transient_validation$phase == "startup_first_5_seconds"
              ]
            ), all, logical(1)))) "yes" else "no",
            if (cooldown_monotonic) "yes" else "no"))
cat(sprintf("Temperature feedback strategy ranking: %s (without feedback: %s)\n",
            feedback_ranking, no_feedback_ranking))
cat(sprintf("Time refinement 1 s to 0.5 s: peak delta %.9g K, energy delta %.9g J\n",
            time_convergence$peak_temperature_difference[1],
            time_convergence$energy_difference_J[1]))
cat(sprintf("Spatial refinement 41 to 81 cells: peak delta %.9g K, energy delta %.9g J\n",
            space_convergence$peak_temperature_difference[1],
            space_convergence$energy_difference_J[1]))
baseline_saving_percent <-
  (summary_table$total_energy_J[summary_table$mode == "low" &
                                  !summary_table$feedback] -
     summary_table$total_energy_J[summary_table$mode == "high" &
                                    !summary_table$feedback]) /
  summary_table$total_energy_J[summary_table$mode == "low" &
                                 !summary_table$feedback] * 100
feedback_saving_percent <-
  (feedback_low_energy - feedback_high_energy) / feedback_low_energy * 100
cat(sprintf("High-mode energy saving: no feedback %.6f%%, feedback %.6f%%\n",
            baseline_saving_percent, feedback_saving_percent))
cat(sprintf("High-mode carbon reduction with fixed intensity: no feedback %.6f%%, feedback %.6f%%\n",
            baseline_saving_percent, feedback_saving_percent))
cat(sprintf("New PDE outputs written to %s\n", pde_output_dir))
