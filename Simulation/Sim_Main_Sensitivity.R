library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(stringr)

set.seed(123)
data <- readRDS("Data/Sim_Data_Sensitivity.rds")

scenarios <- c("a", "b", "c")
N <- c(5000)
alpha_gamma_values <- c(alpha_u = 0.5, gamma_u = -0.5)

R_combinations <- list()
for(i in 0:10) {  
  R_w_val <- -0.5 + i * 0.1  
  R_z_val <- 0.5 - i * 0.1   
  R_combinations[[i+1]] <- c(R_w = R_w_val, R_z = R_z_val)
}

estimate_bitsls <- function(data, R_w, R_z){
  fit_w <- lm(W ~ Z + V + I(Z * V), data = data)
  fit_z <- lm(Z ~ W + V + I(W * V), data = data)
  
  model_x_z <- lm(X ~ Z + fit_w$fitted.values + V, data = data)
  model_y_z <- lm(Y ~ Z + fit_w$fitted.values + V, data = data)
  model_x_w <- lm(X ~ fit_z$fitted.values + W + V, data = data)
  model_y_w <- lm(Y ~ fit_z$fitted.values + W + V, data = data)
  
  k1 <- coef(model_x_w)["W"] / coef(model_y_w)["W"]
  k2 <- coef(model_y_z)["Z"] / coef(model_x_z)["Z"]
  
  beta_xy <- ((k2 * (1 + k1 * R_z - R_w * R_z)) - R_z) / (1 - k1 * k2 * R_w * R_z)
  beta_yx <- ((k1 * (1 + k2 * R_w - R_w * R_z)) - R_w) / (1 - k1 * k2 * R_w * R_z)
  
  c(beta_xy, beta_yx)
}

estimate_iv <- function(data, R_w, R_z){
  model_x <- lm(X ~ Z + W + V, data = data)
  model_y <- lm(Y ~ Z + W + V, data = data)
  
  beta_xy <- coef(model_y)["Z"] / coef(model_x)["Z"]
  beta_yx <- coef(model_x)["W"] / coef(model_y)["W"]
  
  c(beta_xy, beta_yx)
}

estimate_ols <- function(data, R_w, R_z){
  model_x <- lm(X ~ Y + Z + W + V, data = data)
  model_y <- lm(Y ~ X + Z + W + V, data = data)
  
  beta_xy <- coef(model_y)["X"]
  beta_yx <- coef(model_x)["Y"]
  
  c(beta_xy, beta_yx)
}

estimators <- list(
  OLS = estimate_ols,
  IV = estimate_iv,
  Bi_TSLS = estimate_bitsls
)

run_sensitivity_simulation <- function(all_iterations) {
  df_list <- list()
  combo_key <- paste0("alpha_u_", alpha_gamma_values["alpha_u"], "_gamma_u_", alpha_gamma_values["gamma_u"])
  
  for(i in seq_along(all_iterations)) {
    for(scenario in scenarios) {
      for(r_idx in seq_along(R_combinations)) {
        r_values <- R_combinations[[r_idx]]
        r_key <- paste0("R_w_", r_values["R_w"], "_R_z_", r_values["R_z"])
        
        for(n_size in N) {
          n_key <- paste0("N_", n_size)
          current_data <- all_iterations[[i]][[scenario]][[combo_key]][[r_key]][[n_key]]
          
          for(e in seq_along(estimators)) {
            estimates <- estimators[[e]](current_data, r_values["R_w"], r_values["R_z"])
            
            df_list[[length(df_list) + 1]] <- data.frame(
              scenario = scenario,
              R_w = unname(r_values["R_w"]),
              R_z = unname(r_values["R_z"]),
              estimator = names(estimators)[e],
              beta_xy = estimates[1],
              beta_yx = estimates[2]
            )
          }
        }
      }
    }
  }
  
  do.call(rbind, df_list)
}

df_results <- run_sensitivity_simulation(data)

create_sensitivity_plot <- function(df, direction, manual_ylim) {
  beta_col <- ifelse(direction == "X to Y", "beta_xy", "beta_yx")
  true_value <- ifelse(direction == "X to Y", 0.5, -0.5)
  
  main_title <- if (direction == "X to Y") {
    bquote(hat(beta)[X %->% Y] ~ "across different levels of proxy structural conditions violation")
  } else {
    bquote(hat(beta)[Y %->% X] ~ "across different levels of proxy structural conditions violation")
  }
  
  subtitle_text <- "N = 5000"
  
  df$scenario <- factor(df$scenario, 
                        levels = c("a", "b", "c"),
                        labels = c("Scenario (a)", "Scenario (b)", "Scenario (c)"))
  df$estimator <- factor(df$estimator, 
                         levels = c("OLS", "IV", "Bi_TSLS"))
  
  ggplot(df, aes(x = factor(R_w), y = !!sym(beta_col), colour = estimator)) +
    stat_boxplot(geom = "errorbar", width = 0.5, position = position_dodge(width = 0.75)) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.5, outlier.shape = NA) +
    geom_hline(yintercept = true_value, linetype = "dashed", colour = "black", linewidth = 0.5) +
    facet_grid(. ~ scenario) +
    scale_colour_manual(values = c("#2E86C1", "#28B463", "#E74C3C"), labels = c("OLS", "IV", "Bi-TSLS")) +
    scale_x_discrete(name = bquote("Sensitivity parameter values"),
                     labels = as.character(seq(-0.5, 0.5, 0.1))) +
    scale_y_continuous(name = "Effect estimates", labels = scales::label_number(accuracy = 0.01)) +
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
    coord_cartesian(ylim = manual_ylim, clip = "off")
}

plot_xy <- create_sensitivity_plot(df_results, "X to Y", manual_ylim = c(-0.2, 1.2))
plot_yx <- create_sensitivity_plot(df_results, "Y to X", manual_ylim = c(-1.5, 0.1))

combined_plot <- grid.arrange(plot_xy, plot_yx, nrow = 2)
ggsave("Figure/Fig_Main_Sensitivity.eps", combined_plot, width = 18, height = 12, units = "in", dpi = 600, device = "eps")