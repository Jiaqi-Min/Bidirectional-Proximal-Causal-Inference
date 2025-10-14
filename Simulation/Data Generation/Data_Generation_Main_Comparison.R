library(foreach)
library(doParallel)

set.seed(123)

generate_data <- function(N, scenario, alpha_u = 0.5, gamma_u = -0.5, n_time = 5000) {
  chunk_size <- ceiling(N/getDoParWorkers())
  
  alpha_0 <- 1;      gamma_0 <- -1
  beta_yx <- -0.5;   beta_xy <- 0.5
  alpha_v <- 1;      gamma_v <- -1
  alpha_z <- 1;      gamma_w <- 1
  delta_u <- 1;      eta_u <- -1
  delta_v <- -0.5;   eta_v  <- 0.5
  
  v_vec <- rnorm(N + n_time)
  u_vec <- exp(v_vec) + 2 * rbinom(N + n_time, size = 1, prob = 0.5) - 1
  
  if (scenario == "a") {
    epsilon_z <- rnorm(N + n_time)
    epsilon_w <- rnorm(N + n_time)
  } else if (scenario == "b") {
    epsilon_z <- runif(N + n_time, -1, 1)
    epsilon_w <- runif(N + n_time, -1, 1)
  } else if (scenario == "c") {
    epsilon_z <- 2 * rbinom(N + n_time, size = 1, prob = 0.5) - 1
    epsilon_w <- 2 * rbinom(N + n_time, size = 1, prob = 0.5) - 1
  }
  
  z_vec <- 1 + delta_u * u_vec + delta_v * v_vec + epsilon_z
  w_vec <- 1 + eta_u * u_vec + eta_v * v_vec + epsilon_w
  
  error_y <- 2 * rbinom(N + n_time, size = 1, prob = 0.5) - 1
  error_x <- 2 * rbinom(N + n_time, size = 1, prob = 0.5) - 1
  
  simulation_results <- foreach(i = seq(1, N, by = chunk_size), 
                                .combine = rbind, 
                                .packages = c("stats")) %dopar% {
                                  current_chunk <- min(chunk_size, N - i + 1)
                                  chunk_indices <- i:(i + current_chunk - 1)
                                  
                                  results <- matrix(NA, nrow = current_chunk, ncol = 6)
                                  colnames(results) <- c("Y", "X", "Z", "W", "V", "U")
                                  
                                  for (j in 1:current_chunk) {
                                    idx <- chunk_indices[j]
                                    y_series <- numeric(n_time)
                                    x_series <- numeric(n_time)
                                    
                                    y_series[1] <- rnorm(1)
                                    x_series[1] <- rnorm(1)
                                    
                                    for (t in 2:n_time) {
                                      x_series[t] <- alpha_0 + beta_yx * y_series[t-1] + alpha_v * v_vec[idx] + alpha_z * z_vec[idx] + alpha_u * u_vec[idx] + error_x[idx]
                                      y_series[t] <- gamma_0 + beta_xy * x_series[t-1] + gamma_v * v_vec[idx] + gamma_w * w_vec[idx] + gamma_u * u_vec[idx] + error_y[idx]
                                    }
                                    
                                    results[j,] <- c(y_series[n_time], x_series[n_time], z_vec[idx], w_vec[idx], v_vec[idx], u_vec[idx])
                                  }
                                  
                                  results
                                }
  
  as.data.frame(simulation_results)
}

scenarios <- c("a", "b", "c")
N <- c(1000, 2000, 5000)
n_iterations <- 200
alpha_u_val <- 0.5
gamma_u_val <- -0.5

all_iterations <- vector("list", n_iterations)

n_cores <- parallel::detectCores() - 1
cluster <- makeCluster(n_cores)
registerDoParallel(cluster)
clusterExport(cluster, c("generate_data"))

for(i in 1:n_iterations) {
  current_data <- list()
  
  for (scenario in scenarios) {
    current_data[[scenario]] <- list()
    combo_key <- paste0("alpha_u_", alpha_u_val, "_gamma_u_", gamma_u_val)
    current_data[[scenario]][[combo_key]] <- list()
    
    for (n_val in N) {
      results <- generate_data(
        N = n_val,
        scenario = scenario,
        alpha_u = alpha_u_val,
        gamma_u = gamma_u_val
      )
      
      current_data[[scenario]][[combo_key]][[paste0("N_", n_val)]] <- results
    }
  }
  
  all_iterations[[i]] <- current_data
}

stopCluster(cluster)
saveRDS(all_iterations, file = "Data/Sim_Data_Main_Comparison.rds")