# ===============================================================================
# ADSORPTION ISOTHERM CALCULATOR & DATA ANALYZER (R Version)
# Models: Freundlich & Langmuir Adsorption Isotherms
# ===============================================================================
#
# Mathematical Formulations:
# 1. Freundlich Isotherm:
#      Non-linear:  q_e = K_f * C_e^(1/n)
#      Linear form: log10(q_e) = (1/n) * log10(C_e) + log10(K_f)
#      Parameters:  Slope = 1/n  -->  n = 1/Slope
#                   Intercept = log10(K_f)  -->  K_f = 10^Intercept
#
# 2. Langmuir Isotherm:
#      Non-linear:  q_e = (q_max * K_L * C_e) / (1 + K_L * C_e)
#      Linear form: (C_e / q_e) = (1 / (q_max * K_L)) + (1 / q_max) * C_e
#      Parameters:  Slope = 1/q_max  -->  q_max = 1/Slope
#                   Intercept = 1/(q_max * K_L)  -->  K_L = Slope / Intercept
#      Separation Factor: R_L = 1 / (1 + K_L * C_0)
# ===============================================================================

# --- 1. FREUNDLICH ISOTHERM FIT ---
freundlich_fit <- function(C_e, q_e) {
  # Input validation
  valid_idx <- which(C_e > 0 & q_e > 0)
  if (length(valid_idx) < 2) {
    stop("At least 2 data points with C_e > 0 and q_e > 0 are required.")
  }
  
  C_val <- C_e[valid_idx]
  q_val <- q_e[valid_idx]
  
  # Log10 transformation
  log_C <- log10(C_val)
  log_q <- log10(q_val)
  
  # Linear regression model
  fit_lm <- lm(log_q ~ log_C)
  intercept <- coef(fit_lm)[1]
  slope <- coef(fit_lm)[2]
  
  if (slope == 0) {
    stop("Fitted slope is zero; cannot determine Freundlich exponent n.")
  }
  
  one_over_n <- unname(slope)
  n <- 1 / one_over_n
  K_f <- 10^unname(intercept)
  r2_linear <- summary(fit_lm)$r.squared
  
  # Non-linear predictions and real-space R^2
  q_pred <- K_f * (C_val^(1 / n))
  ss_res <- sum((q_val - q_pred)^2)
  ss_tot <- sum((q_val - mean(q_val))^2)
  r2_nonlinear <- ifelse(ss_tot != 0, 1 - (ss_res / ss_tot), 0)
  rmse <- sqrt(mean((q_val - q_pred)^2))
  
  return(list(
    model = "Freundlich",
    K_f = K_f,
    n = n,
    one_over_n = one_over_n,
    slope = slope,
    intercept = intercept,
    r2_linear = r2_linear,
    r2_nonlinear = r2_nonlinear,
    rmse = rmse,
    C_e = C_val,
    q_e = q_val,
    q_pred = q_pred
  ))
}

# --- 2. LANGMUIR ISOTHERM FIT ---
langmuir_fit <- function(C_e, q_e) {
  valid_idx <- which(C_e > 0 & q_e > 0)
  if (length(valid_idx) < 2) {
    stop("At least 2 data points with C_e > 0 and q_e > 0 are required.")
  }
  
  C_val <- C_e[valid_idx]
  q_val <- q_e[valid_idx]
  
  # Linearization: (C_e / q_e) vs C_e
  y_lin <- C_val / q_val
  fit_lm <- lm(y_lin ~ C_val)
  intercept <- coef(fit_lm)[1]
  slope <- coef(fit_lm)[2]
  
  if (slope <= 0) {
    stop("Slope is non-positive; monolayer saturation capacity cannot be calculated.")
  }
  
  q_max <- 1 / unname(slope)
  K_L <- unname(slope) / unname(intercept)
  r2_linear <- summary(fit_lm)$r.squared
  
  # Non-linear predictions & metrics
  q_pred <- (q_max * K_L * C_val) / (1 + K_L * C_val)
  ss_res <- sum((q_val - q_pred)^2)
  ss_tot <- sum((q_val - mean(q_val))^2)
  r2_nonlinear <- ifelse(ss_tot != 0, 1 - (ss_res / ss_tot), 0)
  rmse <- sqrt(mean((q_val - q_pred)^2))
  
  return(list(
    model = "Langmuir",
    q_max = q_max,
    K_L = K_L,
    slope = slope,
    intercept = intercept,
    r2_linear = r2_linear,
    r2_nonlinear = r2_nonlinear,
    rmse = rmse,
    C_e = C_val,
    q_e = q_val,
    q_pred = q_pred
  ))
}

# --- 3. PREDICTION FUNCTIONS ---
predict_freundlich <- function(C_new, K_f, n) {
  return(K_f * (pmax(C_new, 0)^(1 / n)))
}

predict_langmuir <- function(C_new, q_max, K_L) {
  return((q_max * K_L * C_new) / (1 + K_L * C_new))
}

# --- 4. PLOTTING FUNCTION (BASE R - NO EXTERNAL PACKAGES NEEDED) ---
plot_isotherms <- function(C_e, q_e, fr_res, lg_res = NULL, save_filename = NULL) {
  if (!is.null(save_filename)) {
    png(save_filename, width = 1100, height = 500, res = 120)
    on.exit(dev.off())
  }
  
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  
  # Plot 1: q_e vs C_e (Real Space Fit)
  C_seq <- seq(min(C_e) * 0.8, max(C_e) * 1.15, length.out = 200)
  q_fr_seq <- predict_freundlich(C_seq, fr_res$K_f, fr_res$n)
  
  y_max <- max(q_e, q_fr_seq) * 1.15
  plot(C_e, q_e, pch = 19, col = "red", cex = 1.3,
       xlab = expression(bold(paste("Equilibrium Concentration ", C[e], " (mg/L)"))),
       ylab = expression(bold(paste("Adsorption Capacity ", q[e], " (mg/g)"))),
       main = "Adsorption Isotherms Fitting",
       ylim = c(0, y_max), las = 1)
  grid(nx = NULL, ny = NULL, col = "gray85", lty = 2)
  points(C_e, q_e, pch = 19, col = "red", cex = 1.3)
  
  lines(C_seq, q_fr_seq, col = "blue", lwd = 2.5)
  
  legend_labels <- c(
    "Experimental Data",
    sprintf("Freundlich (R² = %.4f)", fr_res$r2_nonlinear)
  )
  legend_colors <- c("red", "blue")
  legend_ltys <- c(NA, 1)
  legend_pchs <- c(19, NA)
  
  if (!is.null(lg_res)) {
    q_lg_seq <- predict_langmuir(C_seq, lg_res$q_max, lg_res$K_L)
    lines(C_seq, q_lg_seq, col = "forestgreen", lwd = 2.5, lty = 2)
    legend_labels <- c(legend_labels, sprintf("Langmuir (R² = %.4f)", lg_res$r2_nonlinear))
    legend_colors <- c(legend_colors, "forestgreen")
    legend_ltys <- c(legend_ltys, 2)
    legend_pchs <- c(legend_pchs, NA)
  }
  
  legend("bottomright", legend = legend_labels, col = legend_colors,
         lty = legend_ltys, pch = legend_pchs, lwd = 2, bty = "o", bg = "white", cex = 0.85)
  
  # Plot 2: Linearized Freundlich Plot (log10(q_e) vs log10(C_e))
  log_C <- log10(fr_res$C_e)
  log_q <- log10(fr_res$q_e)
  
  plot(log_C, log_q, pch = 17, col = "purple", cex = 1.3,
       xlab = expression(bold(paste(log[10], "(", C[e], ")"))),
       ylab = expression(bold(paste(log[10], "(", q[e], ")"))),
       main = expression(bold(paste("Linearized Freundlich: ", log[10](q[e]), " vs ", log[10](C[e])))),
       las = 1)
  grid(nx = NULL, ny = NULL, col = "gray85", lty = 2)
  points(log_C, log_q, pch = 17, col = "purple", cex = 1.3)
  
  abline(a = fr_res$intercept, b = fr_res$slope, col = "darkmagenta", lwd = 2)
  
  legend("topleft",
         legend = c(
           "Transformed Data",
           sprintf("Linear Fit: y = %.3fx + %.3f\n(Linear R² = %.4f)", 
                   fr_res$slope, fr_res$intercept, fr_res$r2_linear)
         ),
         col = c("purple", "darkmagenta"),
         pch = c(17, NA), lty = c(NA, 1), lwd = 2, bty = "o", bg = "white", cex = 0.85)
  
  par(mfrow = c(1, 1))
}


# ===============================================================================
# 5. EXECUTION & DEMONSTRATION WITH EXPERIMENTAL DATA
# ===============================================================================

# Experimental adsorption data
# C_e : equilibrium concentration (mg/L)
# q_e : amount adsorbed per unit mass of adsorbent (mg/g)
C_e_data <- c(10, 20, 30, 40, 50, 60, 70, 80)
q_e_data <- c(2.5, 3.8, 4.9, 5.7, 6.4, 7.0, 7.5, 8.0)

cat(rep("=", 65), "\n", sep = "")
cat("         ADSORPTION ISOTHERM ANALYSIS IN R (FREUNDLICH & LANGMUIR)\n")
cat(rep("=", 65), "\n", sep = "")

# 1. Fit Freundlich Model
fr_res <- freundlich_fit(C_e_data, q_e_data)
cat("\n--- FREUNDLICH ISOTHERM RESULTS ---\n")
cat(sprintf("  Freundlich Constant (K_f)    : %.6f (mg/g)(L/mg)^(1/n)\n", fr_res$K_f))
cat(sprintf("  Freundlich Exponent (n)      : %.6f\n", fr_res$n))
cat(sprintf("  Heterogeneity Factor (1/n)   : %.6f\n", fr_res$one_over_n))
cat(sprintf("  Linear R² (in log10 space)   : %.6f\n", fr_res$r2_linear))
cat(sprintf("  Non-linear R² (real space)   : %.6f\n", fr_res$r2_nonlinear))
cat(sprintf("  Root Mean Squared Error(RMSE): %.6f mg/g\n", fr_res$rmse))

cat("\n  Adsorption Evaluation (Freundlich):\n")
if (fr_res$n > 1) {
  cat(sprintf("  -> n = %.2f > 1 : Favorable adsorption process.\n", fr_res$n))
} else if (fr_res$n == 1) {
  cat("  -> n = 1.0 : Linear partition adsorption.\n")
} else {
  cat(sprintf("  -> n = %.2f < 1 : Unfavorable adsorption process.\n", fr_res$n))
}

# 2. Fit Langmuir Model
lg_res <- langmuir_fit(C_e_data, q_e_data)
cat("\n--- LANGMUIR ISOTHERM RESULTS ---\n")
cat(sprintf("  Monolayer Capacity (q_max)   : %.6f mg/g\n", lg_res$q_max))
cat(sprintf("  Langmuir Constant (K_L)      : %.6f L/mg\n", lg_res$K_L))
cat(sprintf("  Linear R²                    : %.6f\n", lg_res$r2_linear))
cat(sprintf("  Non-linear R² (real space)   : %.6f\n", lg_res$r2_nonlinear))
cat(sprintf("  Root Mean Squared Error(RMSE): %.6f mg/g\n", lg_res$rmse))

# Separation factor at test initial concentration C0 = 80 mg/L
C0_test <- 80
R_L <- 1 / (1 + lg_res$K_L * C0_test)
cat(sprintf("  Separation Factor (R_L @ C0=%.1f mg/L): %.4f (0 < R_L < 1 indicates favorable)\n", C0_test, R_L))

# 3. Model Comparison
cat("\n--- MODEL COMPARISON ---\n")
best_model <- ifelse(fr_res$r2_nonlinear > lg_res$r2_nonlinear, "Freundlich", "Langmuir")
cat(sprintf("  Best fitting model: %s (Freundlich R² = %.4f vs Langmuir R² = %.4f)\n", 
            best_model, fr_res$r2_nonlinear, lg_res$r2_nonlinear))

# 4. Predictions for New Concentrations
cat("\n--- PREDICTIONS FOR NEW CONCENTRATIONS ---\n")
new_Ce <- c(15, 25, 45, 90)
pred_fr <- predict_freundlich(new_Ce, fr_res$K_f, fr_res$n)
pred_lg <- predict_langmuir(new_Ce, lg_res$q_max, lg_res$K_L)

pred_df <- data.frame(
  `C_e (mg/L)` = new_Ce,
  `Freundlich q_e` = round(pred_fr, 4),
  `Langmuir q_e` = round(pred_lg, 4),
  check.names = FALSE
)
print(pred_df, row.names = FALSE)

cat("\n", rep("=", 65), "\n", sep = "")

# 5. Generate plots (saves to isotherm_plot.png & displays in R device)
plot_isotherms(C_e_data, q_e_data, fr_res, lg_res, save_filename = "isotherm_plot.png")
cat("\n[Success] Isotherm fitting plots saved to 'isotherm_plot.png'.\n")
