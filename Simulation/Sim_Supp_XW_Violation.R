library(ggplot2)
library(gridExtra)

set.seed(123)
data <- readRDS("Data/Sim_Data_Supp_XW_Violation.rds")

alpha_u_fixed <- 0.5
gamma_u_fixed <- -0.5
N <- 5000

beta_xw_values <- seq(-0.5, 0.5, by = 0.1)

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

run_beta_xw_simulation <- function(iterations_data) {
  df_list <- list()
  
  for (i in seq_along(iterations_data)) {
    combo_key <- paste0("alpha_u_", alpha_u_fixed, "_gamma_u_", gamma_u_fixed)
    
    for (beta_xw_idx in seq_along(beta_xw_values)) {
      beta_xw_val <- beta_xw_values[beta_xw_idx]
      beta_key <- paste0("beta_XW_", beta_xw_val)
      n_key <- paste0("N_", N)
      
      current_data <- iterations_data[[i]][[combo_key]][[beta_key]][[n_key]]
      
      for (e in seq_along(estimators)) {
        estimates <- estimators[[e]](current_data)
        
        df_list[[length(df_list) + 1]] <- data.frame(
          iteration = i,
          beta_xw = beta_xw_val,
          estimator = names(estimators)[e],
          beta_xy = estimates[1],
          beta_yx = estimates[2]
        )
      }
    }
  }
  
  do.call(rbind, df_list)
}

df_results <- run_beta_xw_simulation(data)

plot_beta_xw <- function(df, direction, manual_ylim) {
  beta_col <- ifelse(direction == "X to Y", "beta_xy", "beta_yx")
  true_value <- ifelse(direction == "X to Y", 0.5, -0.5)
  
  main_title <- if(direction == "X to Y") {
    bquote(hat(beta)[X %->% Y] ~ "across different" ~ beta[X %->% W] ~ "values")
  } else {
    bquote(hat(beta)[Y %->% X] ~ "across different" ~ beta[X %->% W] ~ "values")
  }
  
  subtitle_text <- bquote(alpha[u] == .(alpha_u_fixed) ~ ";" ~ 
                            gamma[u] == .(gamma_u_fixed) ~ ";" ~
                            N == .(N))
  
  df$estimator <- factor(df$estimator, 
                         levels = c("OLS", "IV", "Bi_TSLS"))
  
  ggplot(df, aes(x = factor(beta_xw), y = !!sym(beta_col), color = estimator)) +
    stat_boxplot(geom = 'errorbar', width = 0.4, position = position_dodge(width = 0.7)) +
    stat_boxplot(geom = "boxplot", position = position_dodge(width = 0.7), width = 0.6, outlier.shape = NA) +
    geom_hline(yintercept = true_value, linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_x_discrete(
      name = bquote("Direct Effect from X to W"),
      labels = as.character(beta_xw_values),
      expand = expansion(add = c(0.5, 0.5))
    ) +
    scale_y_continuous(
      name = "Effect estimates",
      labels = scales::label_number(accuracy = 0.01)
    ) +
    scale_color_manual(values = c("#2E86C1", "#28B463", "#E74C3C"), labels = c("OLS", "IV", "Bi-TSLS")) +
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

plot_xy <- plot_beta_xw(df_results, "X to Y", manual_ylim = c(-0.2, 1))
plot_yx <- plot_beta_xw(df_results, "Y to X", manual_ylim = c(-0.7, 0.1))

combined_plot <- grid.arrange(plot_xy, plot_yx, nrow = 2)
ggsave("Figure/Fig_Supp_XW_Variation.eps", combined_plot, width = 15, height = 12, units = "in", dpi = 600, device = "eps")