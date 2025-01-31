library(here)
library(beepr)

source(here("src", "simulation.R"))

LOAD <- TRUE
SAVE <- FALSE

start_time <- proc.time()

if(!LOAD){
  
  set.seed(123)
  N <- 2000
  P <- 100
  A <- 3  # number of active features
  N_SURROGATES <- c(5,5,5)
  SIGMA <- 2
  RHO <- c(rep(0.6, 10), rep(1.2, 5))# 0.59 # control surrogate correlation 0.6
  SNR_STRENGTH_VECTOR <- c(12, 12, 12, 12, 12) # 12
  K <- 10
  KKT_COUNT <- 20 # 30 # at kkt 5 max should be at0.556 <- something is wrong
  LAMBDA_TYPE <- "lambda.min"   # <- sparser model helps, so possibly clustering actually helps
  
  
  output_sim <- simulation_param_grid(param_name="p",               # The name of the parameter to vary
                                 # param_values=c(500,10000),
                                  param_values=c(50, 100, 200, 500, 1000),             # The values of the parameter to iterate over
                                  n_trial=5,                  # Number of trials per simulation
                                  n=N, p=P, a=A, sigma=SIGMA, rho=RHO, snr_strength_vector=SNR_STRENGTH_VECTOR, lambda_type=LAMBDA_TYPE, k=K, kkt_count=KKT_COUNT, pval = FALSE,
                                  n_surrogates=N_SURROGATES)
} else {
  
  load(here("data/simulation_results", "simulation_output_P_5.RData"))
  output_sim <- sim_output$runs 

}

output <- output_sim$final_results

# output <- output %>%
#   select(-c("Sensitivity_Baseline_75", "Precision_Baseline_75"))

# Generate the plots
sensitivity_plot <- plot_metrics_by_param(output, "p", "Sensitivity", log_scale = TRUE) + theme(legend.position = "none")  # Remove legend
precision_plot <- plot_metrics_by_param(output, "p", "Precision", log_scale = TRUE) + theme(legend.position = "none")

# Combine the plots side by side
combined_plot <- sensitivity_plot + precision_plot

# Display the combined plot
print(combined_plot)
beep(2)

ggsave(here("data/plots/combined_plot_p_5.png"), plot = combined_plot, width = 8, height = 8, dpi=600)

if (SAVE){
  sim_output <- list()
  sim_output$runs <- output_sim
  save(sim_output, file = here("data/simulation_results", "simulation_output_P_5.RData"))
  
}


end_time <- proc.time()

execution_time <- end_time - start_time
print(execution_time)

