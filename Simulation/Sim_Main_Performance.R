library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

set.seed(123)
sim_data <- readRDS("Data/Sim_Data_Main.rds")

scenarios <- c("a", "b", "c")
N <- c(1000, 2000, 5000)
alpha_gamma_values <- c(alpha_u = 0.5, gamma_u = -0.5)

estimate_bitsls <- function(data) {
  fit_w <- lm(W ~ Z + V + I(Z*V), data = data)
  fit_z <- lm(Z ~ W + V + I(W*V), data = data)
  
  model_x_z <- lm(X ~ Z + fit_w$fitted.values + V, data = data)
  model_y_z <- lm(Y ~ Z + fit_w$fitted.values + V, data = data)
  model_x_w <- lm(X ~ fit_z$fitted.values + W + V, data = data)
  model_y_w <- lm(Y ~ fit_z$fitted.values + W + V, data = data)
  
  beta_xy <- coef(model_y_z)["Z"] / coef(model_x_z)["Z"]
  beta_yx <- coef(model_x_w)["W"] / coef(model_y_w)["W"]
  
  c(beta_xy, beta_yx)
}

estimate_iv <- function(data) {
  model_x <- lm(X ~ Z + W + V, data = data)
  model_y <- lm(Y ~ Z + W + V, data = data)
  
  beta_xy <- coef(model_y)["Z"] / coef(model_x)["Z"]
  beta_yx <- coef(model_x)["W"] / coef(model_y)["W"]
  
  c(beta_xy, beta_yx)
}

estimate_ols <- function(data) {
  model_x <- lm(X ~ Y + Z + W + V, data = data)
  model_y <- lm(Y ~ X + Z + W + V, data = data)
  
  c(coef(model_y)["X"], coef(model_x)["Y"])
}

estimators <- list(
  OLS = estimate_ols,
  IV = estimate_iv,
  Bi_TSLS = estimate_bitsls
)

run_simulation <- function(iterations_data) {
  combo_key <- paste0("alpha_u_", alpha_gamma_values["alpha_u"], "_gamma_u_", alpha_gamma_values["gamma_u"])
  
  results <- array(NA, dim = c(length(scenarios), length(N), length(estimators), 2, length(iterations_data)))
  
  for (i in seq_along(iterations_data)) {
    for (s in seq_along(scenarios)) {
      for (n in seq_along(N)) {
        n_key <- paste0("N_", N[n])
        current_data <- iterations_data[[i]][[scenarios[s]]][[combo_key]][[n_key]]
        
        for (e in seq_along(estimators)) {
          estimates <- estimators[[e]](current_data)
          results[s, n, e, , i] <- estimates
        }
      }
    }
  }
  
  results
}

results_array <- run_simulation(sim_data)

process_results <- function(results_array) {
  df_list <- list()
  
  for (s in seq_along(scenarios)) {
    for (n in seq_along(N)) {
      for (e in seq_along(estimators)) {
        beta_xy <- results_array[s, n, e, 1, ]
        beta_yx <- results_array[s, n, e, 2, ]
        
        df_list[[length(df_list) + 1]] <- data.frame(
          scenario = scenarios[s],
          N = N[n],
          estimator = names(estimators)[e],
          beta_xy = beta_xy,
          beta_yx = beta_yx
        )
      }
    }
  }
  
  do.call(rbind, df_list)
}

df_results <- process_results(results_array)

create_plot <- function(df, direction, true_value, manual_ylim, sample_size = 5000) {
  beta_var <- ifelse(direction == "XY", "beta_xy", "beta_yx")
  
  main_title <- if(direction == "XY") {
    bquote(hat(beta)[X %->% Y] ~ "across different sample sizes and methods")
  } else {
    bquote(hat(beta)[Y %->% X] ~ "across different sample sizes and methods")
  }
  
  subtitle_text <- bquote(alpha[u] == 0.5 ~ ";" ~ gamma[u] == -0.5)
  
  df$scenario <- factor(df$scenario, levels = c("a", "b", "c"), labels = c("Scenario (a)", "Scenario (b)", "Scenario (c)"))
  df$estimator <- factor(df$estimator, levels = c("OLS", "IV", "Bi_TSLS"))
  df$N <- factor(df$N)
  
  ylim_range <- manual_ylim
  
  ggplot(df, aes(x = N, y = !!sym(beta_var), fill = estimator)) +
    stat_boxplot(geom = 'errorbar', width = 0.4, position = position_dodge(width = 0.7)) +
    stat_boxplot(geom = "boxplot", position = position_dodge(width = 0.7), width = 0.6, outlier.shape = NA) +
    geom_hline(yintercept = true_value, linetype = "dashed", color = "black", linewidth = 0.5) +
    facet_wrap(~scenario, nrow = 1) +
    scale_x_discrete(name = "Sample sizes", expand = expansion(add = c(0.5, 0.5))) +
    scale_y_continuous(name = "Effect estimates", labels = scales::label_number(accuracy = 0.01)) +
    scale_fill_manual(values = c("#a6cee3", "#b2df8a", "#fb9a99"), labels = c("OLS", "IV", "Bi-TSLS")) +
    theme_minimal(base_size = 20) +
    theme(
      panel.background = element_rect(fill = "#F5F5F5"),
      panel.grid.major = element_line(colour = "lightgray"),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 25),
      axis.text = element_text(size = 12),
      legend.position = "right",
      legend.title = element_blank(),
      strip.text = element_text(size = 14),
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 15, hjust = 0.5),
      plot.margin = margin(10, 20, 10, 20)
    ) +
    labs(title = main_title, subtitle = subtitle_text) +
    coord_cartesian(ylim = ylim_range, clip = "off")
}

plot_xy <- create_plot(df_results, "XY", 0.5, manual_ylim = c(-0.1, 1))
plot_yx <- create_plot(df_results, "YX", -0.5, manual_ylim = c(-0.8, 0.1))

combined_plot <- grid.arrange(plot_xy, plot_yx, nrow = 2)
ggsave("Figure/Fig_Main_Results.eps", combined_plot, width = 18, height = 12, units = "in", dpi = 600, device = "eps")