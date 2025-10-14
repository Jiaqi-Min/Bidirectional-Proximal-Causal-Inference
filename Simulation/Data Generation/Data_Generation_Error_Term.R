library(foreach)
library(doParallel)
library(MASS)

set.seed(123)

generate_data <- function(N, scenario, alpha_u, gamma_u, time = 5000) {
  chunk_size <- ceiling(N/getDoParWorkers())
  
  alpha_0 <- 1;      gamma_0 <- -1
  beta_yx <- -0.5;   beta_xy <- 0.5
  alpha_v <- 1;      gamma_v <- -1
  alpha_z <- 1;      gamma_w <- 2
  delta_u <- 1;      eta_u <- -1
  delta_v <- -0.5;   eta_v  <- 0.5
  
  V <- rnorm(N + time)
  U <- exp(V) + 2 * rbinom(N + time, size = 1, prob = 0.5) - 1
  
  if (scenario == "a") {
    Sigma <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
    epsilon_zw <- mvrnorm(n = N + time, mu = c(0, 0), Sigma = Sigma)
    epsilon_z <- epsilon_zw[, 1]
    epsilon_w <- epsilon_zw[, 2]
  } else if (scenario == "b") {
    angles <- runif(N + time, min = 0, max = 2*pi)
    epsilon_z <- cos(angles)
    epsilon_w <- sin(angles)
  } else if (scenario == "c") {
    points <- matrix(c(0, 1, 0, -1, -1, 0, 1, 0), nrow = 4, byrow = TRUE)
    selected_indices <- sample(1:4, size = N + time, replace = TRUE, prob = rep(0.25, 4))
    epsilon_z <- points[selected_indices, 1]
    epsilon_w <- points[selected_indices, 2]
  }
  
  Z <- 1 + delta_u * U + delta_v * V + epsilon_z
  W <- 1 + eta_u * U + eta_v * V + epsilon_w
  
  u_Y <- 2 * rbinom(N + time, size = 1, prob = 0.5) - 1
  u_X <- 2 * rbinom(N + time, size = 1, prob = 0.5) - 1
  
  simulation_results <- foreach(i = seq(1, N, by = chunk_size), 
                                .combine = rbind, 
                                .packages = c("stats")) %dopar% {
                                  current_chunk <- min(chunk_size, N - i + 1)
                                  chunk_indices <- i:(i + current_chunk - 1)
                                  
                                  results <- matrix(NA, nrow = current_chunk, ncol = 6)
                                  colnames(results) <- c("Y", "X", "Z", "W", "V", "U")
                                  
                                  for (j in 1:current_chunk) {
                                    idx <- chunk_indices[j]
                                    Y <- numeric(time)
                                    X <- numeric(time)
                                    
                                    Y[1] <- rnorm(1)
                                    X[1] <- rnorm(1)
                                    
                                    for (t in 2:time) {
                                      X[t] <- alpha_0 + beta_yx * Y[t-1] + alpha_v * V[idx] + alpha_z * Z[idx] + alpha_u * U[idx] + u_X[idx]
                                      Y[t] <- gamma_0 + beta_xy * X[t-1] + gamma_v * V[idx] + gamma_w * W[idx] + gamma_u * U[idx] + u_Y[idx]
                                    }
                                    
                                    results[j,] <- c(Y[time], X[time], Z[idx], W[idx], V[idx], U[idx])
                                  }
                                  
                                  results
                                }
  
  as.data.frame(simulation_results)
}

scenarios <- c("a", "b", "c")
alpha_gamma_combinatioN <- list(c(alpha_u = 0.5, gamma_u = -0.5))
N <- c(1000, 2000, 5000)
iter <- 200

all_iteratioN <- vector("list", iter)

cores <- parallel::detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterExport(cl, c("generate_data"))
clusterEvalQ(cl, library(MASS))

for(i in 1:iter) {
  current_data <- list()
  
  for (scenario in scenarios) {
    current_data[[scenario]] <- list()
    
    for (combo in alpha_gamma_combinatioN) {
      alpha_u_val <- combo["alpha_u"]
      gamma_u_val <- combo["gamma_u"]
      combo_key <- paste0("alpha_u_", alpha_u_val, "_gamma_u_", gamma_u_val)
      current_data[[scenario]][[combo_key]] <- list()
      
      for (N_val in N) {
        results <- generate_data(
          N = N_val,
          scenario = scenario,
          alpha_u = alpha_u_val,
          gamma_u = gamma_u_val
        )
        
        current_data[[scenario]][[combo_key]][[paste0("N_", N_val)]] <- results
      }
    }
  }
  
  all_iteratioN[[i]] <- current_data
}

stopCluster(cl)
saveRDS(all_iteratioN, file = "Data/Sim_Data_Error_Term.rds")