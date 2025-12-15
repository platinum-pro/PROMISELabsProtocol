library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(gridExtra)

# Read the data
data <- read.csv("~/Downloads/Agg_Percentage_Analysis (Dec 2025).csv")

# Convert all columns to numeric (except x which should already be)
data <- data %>% mutate(across(-x, as.numeric))

# Define the demographic pairs to analyze
demographic_pairs <- list(
  list(var = "Economy", label1 = "Access to multiple sellers", label0 = "Access to single seller"),
  list(var = "GENDER", label1 = "Female", label0 = "Male"),
  list(var = "AGE_binary", label1 = "Above median age", label0 = "Below median age"),
  list(var = "RELIGION_binary", label1 = "High religiosity", label0 = "Low religiosity"),
  list(var = "INC", label1 = "Higher income", label0 = "Lower income"),
  list(var = "HINS", label1 = "Insured", label0 = "Uninsured"),
  list(var = "ALC", label1 = "History of drinking alcohol", label0 = "No history of drinking alcohol"),
  list(var = "GB", label1 = "History of gambling", label0 = "No history of gambling")
)

# GLOBAL k: Fixed at 2 for all analyses
k_global <- 2

cat("######################################\n")
cat("GLOBAL PARAMETER\n")
cat("######################################\n")
cat("k (fixed globally) =", k_global, "\n")
cat("This sets the implied floor near ~1% of Q0 under the model form\n")
cat("k is treated as a property of the measurement instrument\n\n")

# Define the Koff function with FIXED k
Koff <- function(x, alpha, Qo, k_fixed) {
  Qo * 10^(k_fixed * (exp(-alpha * Qo * x) - 1))
}

# Function to calculate P50 using numerical root-finding with FIXED k
calculate_p50 <- function(alpha, Qo, k_fixed, data_subset, target_q = 50) {
  if(is.na(k_fixed) || alpha <= 0 || Qo <= 0 || k_fixed <= 0) {
    return(NA)
  }
  
  objective_function <- function(price) {
    q_pred <- Qo * 10^(k_fixed * (exp(-alpha * Qo * price) - 1))
    return(q_pred - target_q)
  }
  
  q_at_zero <- Qo * 10^(k_fixed * (exp(0) - 1))
  
  if(q_at_zero < target_q) {
    warning("Q0 is below target Q - P50 cannot be calculated for this group")
    return(NA)
  }
  
  lower_bound <- min(data_subset$x[data_subset$x > 0], na.rm = TRUE)
  upper_bound <- max(data_subset$x, na.rm = TRUE)
  
  f_lower <- objective_function(lower_bound)
  f_upper <- objective_function(upper_bound)
  
  if(f_lower * f_upper > 0) {
    warning("Target Q is not crossed within the observed price range for this group")
    return(NA)
  }
  
  result <- tryCatch({
    uniroot(objective_function, interval = c(lower_bound, upper_bound), tol = 0.01)
  }, error = function(e) {
    warning("Error in root-finding: ", e$message)
    return(NULL)
  })
  
  if(!is.null(result)) {
    return(result$root)
  } else {
    return(NA)
  }
}

# Function to fit model for a demographic subgroup
fit_demographic_group <- function(data_subset, group_label, k_fixed) {
  valid_data <- data_subset[!is.na(data_subset$value) & data_subset$x > 0, ]
  
  if(nrow(valid_data) < 3) {
    message("Not enough data for ", group_label)
    return(NULL)
  }
  
  cat("\n=== ", group_label, " ===\n", sep = "")
  
  fit <- tryCatch({
    nls(formula = value ~ Koff(x, alpha, Qo, k_fixed),
        data = valid_data,
        start = list(alpha = 0.0000001, Qo = 100),
        algorithm = "port",
        lower = c(alpha = 0, Qo = 0),
        upper = c(alpha = 0.1, Qo = 100),
        control = nls.control(maxiter = 50000))
  }, error = function(e) {
    message("Error fitting model for ", group_label, ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(fit)) {
    residuals <- residuals(fit)
    tss <- sum((valid_data$value - mean(valid_data$value, na.rm = TRUE))^2, na.rm = TRUE)
    rss <- sum(residuals^2)
    r_squared <- 1 - (rss / tss)
    
    params <- coef(fit)
    p50 <- calculate_p50(params["alpha"], params["Qo"], k_fixed, valid_data, target_q = 50)
    
    cat("R-squared =", round(r_squared, 4), "\n")
    cat("alpha =", format(params["alpha"], scientific = TRUE), "\n")
    cat("Q0 =", round(params["Qo"], 2), "\n")
    
    if(!is.na(p50)) {
      cat("P50 = ₦", formatC(p50, format = "f", big.mark = ",", digits = 2), "\n", sep = "")
    } else {
      cat("P50 could not be calculated\n")
    }
    
    return(list(fit = fit, r_squared = r_squared, p50 = p50, label = group_label, 
                alpha = params["alpha"], Qo = params["Qo"], k = k_fixed))
  }
  
  return(NULL)
}

# Function to create comparison plot WITH OVERALL CURVE
create_comparison_plot <- function(data, demographic_pair, fit_results_1, fit_results_0, k_fixed, 
                                   overall_data, show_title = TRUE) {
  col_1 <- paste0(demographic_pair$var, "_1")
  col_0 <- paste0(demographic_pair$var, "_0")
  
  data_1 <- data.frame(x = data$x, value = data[[col_1]]) %>% filter(!is.na(value) & x > 0)
  data_0 <- data.frame(x = data$x, value = data[[col_0]]) %>% filter(!is.na(value) & x > 0)
  
  x_range <- 10^seq(log10(min(data$x[data$x > 0])), log10(max(data$x)), length.out = 100)
  plot_data <- data.frame(x = x_range)
  
  # Add overall curve prediction
  if(!is.null(overall_data$fit)) {
    params_overall <- coef(overall_data$fit)
    plot_data$overall_pred <- Koff(x_range, alpha = params_overall["alpha"], 
                                   Qo = params_overall["Qo"], k_fixed = k_fixed)
  }
  
  if (!is.null(fit_results_1$fit)) {
    params_1 <- coef(fit_results_1$fit)
    plot_data$group1_pred <- Koff(x_range, alpha = params_1["alpha"], Qo = params_1["Qo"], k_fixed = k_fixed)
  }
  
  if (!is.null(fit_results_0$fit)) {
    params_0 <- coef(fit_results_0$fit)
    plot_data$group0_pred <- Koff(x_range, alpha = params_0["alpha"], Qo = params_0["Qo"], k_fixed = k_fixed)
  }
  
  p <- ggplot()
  
  # Add OVERALL curve FIRST (so it's in the background) - BLACK DASHED
  if(!is.null(overall_data$fit)) {
    p <- p + geom_line(data = plot_data, aes(x = x, y = overall_pred), 
                       color = "black", linewidth = 0.6, alpha = 1, linetype = "dashed")
  }
  
  if(nrow(data_1) > 0) {
    p <- p + geom_point(data = data_1, aes(x = x, y = value, color = demographic_pair$label1), 
                        shape = 16, alpha = 0.6, size = 1.5)
  }
  
  if(nrow(data_0) > 0) {
    p <- p + geom_point(data = data_0, aes(x = x, y = value, color = demographic_pair$label0), 
                        shape = 17, alpha = 0.6, size = 1.5)
  }
  
  if(!is.null(fit_results_1$fit)) {
    p <- p + geom_line(data = plot_data, aes(x = x, y = group1_pred, color = demographic_pair$label1), linewidth = 0.8)
  }
  
  if(!is.null(fit_results_0$fit)) {
    p <- p + geom_line(data = plot_data, aes(x = x, y = group0_pred, color = demographic_pair$label0), linewidth = 0.8)
  }
  
  p <- p + geom_hline(yintercept = 50, linetype = "dotted", color = "grey50", alpha = 0.7)
  
  colors <- scales::hue_pal()(2)
  
  if (!is.null(fit_results_1$p50) && !is.na(fit_results_1$p50)) {
    p <- p + geom_vline(xintercept = fit_results_1$p50, linetype = "dashed", color = colors[1], alpha = 0.4)
  }
  
  if (!is.null(fit_results_0$p50) && !is.na(fit_results_0$p50)) {
    p <- p + geom_vline(xintercept = fit_results_0$p50, linetype = "dashed", color = colors[2], alpha = 0.4)
  }
  
  plot_title <- if(show_title) {
    var_titles <- list("Economy" = "Economy Type", "GENDER" = "Gender", "AGE_binary" = "Age", 
                       "RELIGION_binary" = "Religiosity", "INC" = "Income",
                       "ALC" = "Alcohol Use", "GB" = "Gambling")
    var_titles[[demographic_pair$var]]
  } else ""
  
  p <- p +
    scale_x_log10(breaks = c(100, 1000, 10000, 100000), labels = scales::comma) +
    scale_color_manual(values = colors, name = "") +
    labs(title = plot_title, x = "Price (₦)", y = "% Purchasing") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
          axis.title = element_text(size = 9),
          axis.text = element_text(size = 8),
          legend.position = "bottom",
          legend.text = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"),
          plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")) +
    annotation_logticks(sides = "b", size = 0.3)
  
  return(p)
}

# CREATE THE "ALL" PLOT FIRST (need it for overlay)
cat("\n\n######################################\n")
cat("Creating 'All' plot\n")
cat("######################################\n")

valid_data_all <- data[!is.na(data$All) & data$x > 0, ]

fit_all <- tryCatch({
  nls(formula = All ~ Koff(x, alpha, Qo, k_global),
      data = valid_data_all,
      start = list(alpha = 0.0000001, Qo = 100),
      algorithm = "port",
      lower = c(alpha = 0, Qo = 0),
      upper = c(alpha = 0.1, Qo = 100),
      control = nls.control(maxiter = 50000))
}, error = function(e) NULL)

overall_fit_data <- NULL

if (!is.null(fit_all)) {
  params_all <- coef(fit_all)
  p50_all <- calculate_p50(params_all["alpha"], params_all["Qo"], k_global, valid_data_all, target_q = 50)
  
  # Store overall fit data for use in subgroup plots
  overall_fit_data <- list(fit = fit_all, alpha = params_all["alpha"], Qo = params_all["Qo"])
}

# Main analysis loop for demographic pairs
all_results <- list()

for(i in seq_along(demographic_pairs)) {
  demo_pair <- demographic_pairs[[i]]
  
  cat("\n\n######################################\n")
  cat("Analyzing:", demo_pair$var, "\n")
  cat("######################################\n")
  
  col_1 <- paste0(demo_pair$var, "_1")
  col_0 <- paste0(demo_pair$var, "_0")
  
  if(!col_1 %in% names(data) || !col_0 %in% names(data)) {
    cat("Columns", col_1, "or", col_0, "not found in data. Skipping.\n")
    next
  }
  
  data_1 <- data.frame(x = data$x, value = data[[col_1]])
  data_0 <- data.frame(x = data$x, value = data[[col_0]])
  
  fit_1 <- fit_demographic_group(data_1, demo_pair$label1, k_global)
  fit_0 <- fit_demographic_group(data_0, demo_pair$label0, k_global)
  
  if(!is.null(fit_1) || !is.null(fit_0)) {
    # Pass overall_fit_data to create_comparison_plot
    p <- create_comparison_plot(data, demo_pair, fit_1, fit_0, k_global, overall_fit_data, show_title = TRUE)
    all_results[[demo_pair$var]] <- list(group1 = fit_1, group0 = fit_0, plot = p)
  }
}

# Create grid with selected demographic plots (excluding HINS)
selected_vars <- c("Economy", "GENDER", "AGE_binary", "RELIGION_binary", "INC", "ALC", "GB")
grid_plots <- list()

# Add demographic plots
for(var in selected_vars) {
  if(var %in% names(all_results)) {
    grid_plots[[length(grid_plots) + 1]] <- all_results[[var]]$plot
  }
}

# Arrange in grid - 7 plots: 4 rows x 2 columns (last row has 1 plot)
if(length(grid_plots) > 0) {
  grid_combined <- grid.arrange(grobs = grid_plots, ncol = 2, nrow = 4)
  
  # Save the grid
  #ggsave("~/Downloads/demand_curves_grid_with_overall.png", grid_combined, 
         #width = 10, height = 16, dpi = 300)
  
  cat("\n\nGrid plot with", length(grid_plots), "panels saved and displayed\n")
}

# Create summary table of P50 values
p50_summary <- data.frame(
  Variable = character(), Group = character(), P50 = numeric(),
  Alpha = numeric(), Q0 = numeric(), k = numeric(), R_squared = numeric(),
  stringsAsFactors = FALSE
)

for(var_name in names(all_results)) {
  result <- all_results[[var_name]]
  
  if(!is.null(result$group1)) {
    p50_summary <- rbind(p50_summary, 
                         data.frame(Variable = var_name, Group = result$group1$label,
                                    P50 = ifelse(is.na(result$group1$p50), NA, result$group1$p50),
                                    Alpha = result$group1$alpha, Q0 = result$group1$Qo,
                                    k = result$group1$k, R_squared = result$group1$r_squared))
  }
  
  if(!is.null(result$group0)) {
    p50_summary <- rbind(p50_summary,
                         data.frame(Variable = var_name, Group = result$group0$label,
                                    P50 = ifelse(is.na(result$group0$p50), NA, result$group0$p50),
                                    Alpha = result$group0$alpha, Q0 = result$group0$Qo,
                                    k = result$group0$k, R_squared = result$group0$r_squared))
  }
}

cat("\n\n######################################\n")
cat("P50 SUMMARY TABLE\n")
cat("######################################\n")
cat("Note: k = ", k_global, " (fixed globally for all groups)\n", sep = "")
cat("All P50 values are directly comparable on the same scale\n\n")
print(p50_summary, row.names = FALSE)