library(dplyr)

# =========================
# READ + BASIC CLEANING
# =========================
data <- read.csv("~/Downloads/Agg_Percentage_Analysis (Dec 2025).csv")

# Ensure numeric columns (except x)
data <- data %>% mutate(across(-x, as.numeric))

# =========================
# GLOBAL k (fixed)
# =========================
k_global <- 2

cat("######################################\n")
cat("GLOBAL PARAMETER\n")
cat("######################################\n")
cat("k (fixed globally) =", k_global, "\n")
cat("This assumes all groups asymptote to ~1% of Q0\n")
cat("k is treated as a property of the measurement instrument\n\n")

# =========================
# DEMAND FUNCTION (Koff), fixed k
# =========================
Koff <- function(x, alpha, Qo, k_fixed) {
  Qo * 10^(k_fixed * (exp(-alpha * Qo * x) - 1))
}

# =========================
# P50 CALCULATION (fixed k)
# =========================
calculate_p50 <- function(alpha, Qo, k_fixed, x_values, target_q = 50) {
  if (is.na(k_fixed) || is.na(alpha) || is.na(Qo) || alpha <= 0 || Qo <= 0 || k_fixed <= 0) return(NA)
  if (Qo < target_q) return(NA)
  
  objective <- function(price) Koff(price, alpha = alpha, Qo = Qo, k_fixed = k_fixed) - target_q
  
  lower_bound <- min(x_values[x_values > 0], na.rm = TRUE)
  upper_bound <- max(x_values, na.rm = TRUE)
  
  if (!is.finite(lower_bound) || !is.finite(upper_bound) || lower_bound >= upper_bound) return(NA)
  
  f_lower <- objective(lower_bound)
  f_upper <- objective(upper_bound)
  if (!is.finite(f_lower) || !is.finite(f_upper) || f_lower * f_upper > 0) return(NA)
  
  result <- tryCatch(
    uniroot(objective, interval = c(lower_bound, upper_bound), tol = 0.01),
    error = function(e) NULL
  )
  
  if (is.null(result)) NA else result$root
}

# =========================
# SINGLE-GROUP FIT + R2 + P50
# =========================
fit_and_calculate_rsq <- function(condition, data, k_fixed, start_alpha = 1e-7, start_Qo = 90) {
  valid_data <- data %>%
    select(x, y = all_of(condition)) %>%
    filter(x > 0, !is.na(y))
  
  if (nrow(valid_data) < 3) {
    message("Not enough data for ", condition)
    return(NULL)
  }
  
  fit <- tryCatch({
    nls(
      y ~ Koff(x, alpha, Qo, k_fixed),
      data = valid_data,
      start = list(alpha = start_alpha, Qo = start_Qo),
      algorithm = "port",
      lower = c(alpha = 0, Qo = 0),
      upper = c(alpha = 0.1, Qo = 100),
      control = nls.control(maxiter = 50000, warnOnly = TRUE)
    )
  }, error = function(e) {
    message("Error fitting model for ", condition, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(fit)) return(NULL)
  
  resid <- residuals(fit)
  tss <- sum((valid_data$y - mean(valid_data$y, na.rm = TRUE))^2, na.rm = TRUE)
  rss <- sum(resid^2)
  r_squared <- ifelse(is.finite(tss) && tss > 0, 1 - (rss / tss), NA)
  
  params <- coef(fit)
  p50 <- calculate_p50(params["alpha"], params["Qo"], k_fixed, x_values = valid_data$x, target_q = 50)
  
  list(
    fit = fit,
    r_squared = r_squared,
    p50 = p50,
    alpha = unname(params["alpha"]),
    Qo = unname(params["Qo"]),
    k = k_fixed
  )
}

# =========================
# NESTED MODEL F-TESTS (fit full once)
# =========================
perform_parameter_F_tests <- function(data, col_1, col_0, k_fixed, start_alpha = 1e-7, start_Qo = 90) {
  
  # Require paired rows for fair comparison
  valid_data <- data %>%
    filter(x > 0, !is.na(.data[[col_1]]), !is.na(.data[[col_0]]))
  
  if (nrow(valid_data) < 3) {
    message("Not enough paired data for F-tests.")
    return(list(alpha_test = NULL, Q0_test = NULL, model_full = NULL))
  }
  
  long_data <- data.frame(
    x = rep(valid_data$x, 2),
    y = c(valid_data[[col_1]], valid_data[[col_0]]),
    condition = factor(rep(c("group_1", "group_0"), each = nrow(valid_data)))
  )
  
  # FULL model: different alpha AND different Qo
  model_full <- tryCatch({
    nls(
      y ~ ifelse(condition == "group_1",
                 Koff(x, alpha_1, Qo_1, k_fixed),
                 Koff(x, alpha_0, Qo_0, k_fixed)),
      data = long_data,
      start = list(alpha_1 = start_alpha, alpha_0 = start_alpha, Qo_1 = start_Qo, Qo_0 = start_Qo),
      algorithm = "port",
      lower = c(alpha_1 = 0, alpha_0 = 0, Qo_1 = 0, Qo_0 = 0),
      upper = c(alpha_1 = 0.1, alpha_0 = 0.1, Qo_1 = 100, Qo_0 = 100),
      control = nls.control(maxiter = 50000, warnOnly = TRUE)
    )
  }, error = function(e) NULL)
  
  if (is.null(model_full)) {
    cat("\nFull model failed to converge; cannot run nested F-tests.\n")
    return(list(alpha_test = NULL, Q0_test = NULL, model_full = NULL))
  }
  
  # Restricted for alpha: same alpha, different Qo
  model_same_alpha <- tryCatch({
    nls(
      y ~ ifelse(condition == "group_1",
                 Koff(x, alpha, Qo_1, k_fixed),
                 Koff(x, alpha, Qo_0, k_fixed)),
      data = long_data,
      start = list(alpha = start_alpha, Qo_1 = start_Qo, Qo_0 = start_Qo),
      algorithm = "port",
      lower = c(alpha = 0, Qo_1 = 0, Qo_0 = 0),
      upper = c(alpha = 0.1, Qo_1 = 100, Qo_0 = 100),
      control = nls.control(maxiter = 50000, warnOnly = TRUE)
    )
  }, error = function(e) NULL)
  
  # Restricted for Q0: same Qo, different alpha
  model_same_Q0 <- tryCatch({
    nls(
      y ~ ifelse(condition == "group_1",
                 Koff(x, alpha_1, Qo, k_fixed),
                 Koff(x, alpha_0, Qo, k_fixed)),
      data = long_data,
      start = list(Qo = start_Qo, alpha_1 = start_alpha, alpha_0 = start_alpha),
      algorithm = "port",
      lower = c(Qo = 0, alpha_1 = 0, alpha_0 = 0),
      upper = c(Qo = 100, alpha_1 = 0.1, alpha_0 = 0.1),
      control = nls.control(maxiter = 50000, warnOnly = TRUE)
    )
  }, error = function(e) NULL)
  
  # Helper to compute nested-model F test
  nested_F <- function(model_restricted, model_full) {
    RSS_r <- sum(residuals(model_restricted)^2)
    RSS_f <- sum(residuals(model_full)^2)
    
    df_r <- df.residual(model_restricted)
    df_f <- df.residual(model_full)
    
    if (!is.finite(RSS_r) || !is.finite(RSS_f) || df_r <= df_f) return(NULL)
    
    num_df <- df_r - df_f
    den_df <- df_f
    F_stat <- ((RSS_r - RSS_f) / num_df) / (RSS_f / den_df)
    if (!is.finite(F_stat) || F_stat < 0) return(NULL)
    
    p_value <- pf(F_stat, num_df, den_df, lower.tail = FALSE)
    list(F_stat = F_stat, p_value = p_value, df1 = num_df, df2 = den_df)
  }
  
  alpha_test <- if (!is.null(model_same_alpha)) nested_F(model_same_alpha, model_full) else NULL
  Q0_test    <- if (!is.null(model_same_Q0))    nested_F(model_same_Q0,    model_full) else NULL
  
  cat("\nAlpha (Elasticity) Comparison:\n")
  if (!is.null(alpha_test)) {
    cat("F statistic =", round(alpha_test$F_stat, 3), "\n")
    cat("df =", alpha_test$df1, ",", alpha_test$df2, "\n")
    cat("p-value =", sprintf("%.3f", alpha_test$p_value), "\n\n")
    cat("Restricted (same alpha) coefficients:\n")
    print(coef(model_same_alpha))
    cat("\nFull (different alpha) coefficients:\n")
    print(coef(model_full))
  } else {
    cat("Could not compute alpha comparison (restricted model failed or invalid).\n")
  }
  
  cat("\nQ0 (Intensity) Comparison:\n")
  if (!is.null(Q0_test)) {
    cat("F statistic =", round(Q0_test$F_stat, 3), "\n")
    cat("df =", Q0_test$df1, ",", Q0_test$df2, "\n")
    cat("p-value =", sprintf("%.3f", Q0_test$p_value), "\n\n")
    cat("Restricted (same Q0) coefficients:\n")
    print(coef(model_same_Q0))
    cat("\nFull (different Q0) coefficients:\n")
    print(coef(model_full))
  } else {
    cat("Could not compute Q0 comparison (restricted model failed or invalid).\n")
  }
  
  cat("\nNote: P50 is derived from alpha and Q0 (k fixed), not directly F-tested here.\n")
  
  list(
    alpha_test = alpha_test,
    Q0_test = Q0_test,
    model_full = model_full,
    model_same_alpha = model_same_alpha,
    model_same_Q0 = model_same_Q0
  )
}

# =========================
# DEFINE DEMOGRAPHIC PAIRS
# =========================
demographic_pairs <- list(
  list(var = "Economy", label1 = "Open economy", label0 = "Closed economy"),
  list(var = "GENDER", label1 = "Female", label0 = "Male"),
  list(var = "AGE_binary", label1 = "Above median age", label0 = "Below median age"),
  list(var = "RELIGION_binary", label1 = "High religiosity", label0 = "Low religiosity"),
  list(var = "INC", label1 = "Higher income", label0 = "Lower income"),
  list(var = "HINS", label1 = "Insured", label0 = "Uninsured"),
  list(var = "ALC", label1 = "History of alcohol use", label0 = "No history of alcohol use"),
  list(var = "GB", label1 = "Gambling history", label0 = "No gambling history")
)

# =========================
# MAIN ANALYSIS LOOP
# =========================
all_results <- list()

# FIT SUMMARY TABLE
summary_table <- data.frame(
  Variable = character(),
  Group = character(),
  Alpha = numeric(),
  Q0 = numeric(),
  k = numeric(),
  P50 = numeric(),
  R_squared = numeric(),
  stringsAsFactors = FALSE
)

# F-TEST TABLE (SEPARATE)
ftest_table <- data.frame(
  Variable = character(),
  Test = character(),      # "alpha" or "Q0"
  F_stat = numeric(),
  df1 = numeric(),
  df2 = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(demographic_pairs)) {
  demo_pair <- demographic_pairs[[i]]
  
  cat("\n\n######################################\n")
  cat("Analyzing:", demo_pair$var, "\n")
  cat("######################################\n")
  
  # Get column names
  col_1 <- paste0(demo_pair$var, "_1")
  col_0 <- paste0(demo_pair$var, "_0")
  
  # Check if columns exist
  if (!col_1 %in% names(data) || !col_0 %in% names(data)) {
    cat("Columns", col_1, "or", col_0, "not found in data. Skipping.\n")
    next
  }
  
  # Fit individual models
  fit_1 <- fit_and_calculate_rsq(col_1, data, k_global)
  fit_0 <- fit_and_calculate_rsq(col_0, data, k_global)
  
  # INLINE CONSOLE OUTPUT
  if (!is.null(fit_1)) {
    cat("\nGroup 1:", demo_pair$label1, "\n")
    cat("R² =", round(fit_1$r_squared, 4), 
        " alpha =", format(fit_1$alpha, scientific = TRUE),
        " Q0 =", round(fit_1$Qo, 2),
        " P50 =", ifelse(is.na(fit_1$p50), "NA", 
                         paste0("₦", formatC(fit_1$p50, format = "f", big.mark = ",", digits = 2))),
        "\n")
  } else {
    cat("\nGroup 1:", demo_pair$label1, "fit failed.\n")
  }
  
  if (!is.null(fit_0)) {
    cat("Group 0:", demo_pair$label0, "\n")
    cat("R² =", round(fit_0$r_squared, 4), 
        " alpha =", format(fit_0$alpha, scientific = TRUE),
        " Q0 =", round(fit_0$Qo, 2),
        " P50 =", ifelse(is.na(fit_0$p50), "NA", 
                         paste0("₦", formatC(fit_0$p50, format = "f", big.mark = ",", digits = 2))),
        "\n")
  } else {
    cat("Group 0:", demo_pair$label0, "fit failed.\n")
  }
  
  # Add to summary table
  if (!is.null(fit_1)) {
    summary_table <- rbind(summary_table, 
                           data.frame(
                             Variable = demo_pair$var,
                             Group = demo_pair$label1,
                             Alpha = fit_1$alpha,
                             Q0 = fit_1$Qo,
                             k = fit_1$k,
                             P50 = ifelse(is.na(fit_1$p50), NA, fit_1$p50),
                             R_squared = fit_1$r_squared
                           ))
  }
  
  if (!is.null(fit_0)) {
    summary_table <- rbind(summary_table,
                           data.frame(
                             Variable = demo_pair$var,
                             Group = demo_pair$label0,
                             Alpha = fit_0$alpha,
                             Q0 = fit_0$Qo,
                             k = fit_0$k,
                             P50 = ifelse(is.na(fit_0$p50), NA, fit_0$p50),
                             R_squared = fit_0$r_squared
                           ))
  }
  
  # Perform F-tests
  cat("\n=== Parameter Comparison Tests ===\n")
  cat("Note: k =", k_global, "(fixed globally for all groups)\n")
  parameter_tests <- perform_parameter_F_tests(data, col_1, col_0, k_global)
  
  # Add to F-test table
  if (!is.null(parameter_tests$alpha_test)) {
    ftest_table <- rbind(ftest_table, data.frame(
      Variable = demo_pair$var,
      Test = "alpha",
      F_stat = parameter_tests$alpha_test$F_stat,
      df1 = parameter_tests$alpha_test$df1,
      df2 = parameter_tests$alpha_test$df2,
      p_value = parameter_tests$alpha_test$p_value
    ))
  } else {
    ftest_table <- rbind(ftest_table, data.frame(
      Variable = demo_pair$var,
      Test = "alpha",
      F_stat = NA,
      df1 = NA,
      df2 = NA,
      p_value = NA
    ))
  }
  
  if (!is.null(parameter_tests$Q0_test)) {
    ftest_table <- rbind(ftest_table, data.frame(
      Variable = demo_pair$var,
      Test = "Q0",
      F_stat = parameter_tests$Q0_test$F_stat,
      df1 = parameter_tests$Q0_test$df1,
      df2 = parameter_tests$Q0_test$df2,
      p_value = parameter_tests$Q0_test$p_value
    ))
  } else {
    ftest_table <- rbind(ftest_table, data.frame(
      Variable = demo_pair$var,
      Test = "Q0",
      F_stat = NA,
      df1 = NA,
      df2 = NA,
      p_value = NA
    ))
  }
  
  # Store results
  all_results[[demo_pair$var]] <- list(
    group1 = fit_1,
    group0 = fit_0,
    tests = parameter_tests
  )
}

# =========================
# PRINT SUMMARY TABLES
# =========================
cat("\n\n######################################\n")
cat("FIT SUMMARY TABLE (per group)\n")
cat("######################################\n")
cat("Note: k =", k_global, "(fixed globally for all groups)\n")
cat("All P50 values are directly comparable on the same scale\n\n")
print(summary_table, row.names = FALSE)

cat("\n\n######################################\n")
cat("NESTED F-TEST TABLE (per comparison)\n")
cat("######################################\n")
cat("Note: NA means the relevant restricted/full model failed to converge.\n\n")
print(ftest_table, row.names = FALSE)

