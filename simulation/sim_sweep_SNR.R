library(here)
library(beepr)

source(here("src", "simulation.R"))
source(here("src", "plotting.R"))


LOAD <- FALSE
SAVE <- TRUE

start_time <- proc.time()

if(!LOAD){
  
  set.seed(123)
  N <- 1000
  P <- 100
  A <- 3  # number of active features
  N_SURROGATES <- c(5,5,5)
  SIGMA <- 2
  RHO <- c(rep(0.6, 10), rep(0.9, 5))# 0.59 # control surrogate correlation 0.6
  SNR_STRENGTH_VECTOR <- c(5, 10, 15, 20, 30) # 12
  K <- 10
  KKT_COUNT <- 20 # 30 # at kkt 5 max should be at0.556 <- something is wrong
  LAMBDA_TYPE <- "lambda.min"   # <- sparser model helps, so possibly clustering actually helps
  
  
  output_sim <- simulation_param_grid(param_name="n",               # The name of the parameter to vary
                                 # param_values=c(500,10000),
                                  param_values=c(1000, 1001, 1002, 1003, 1004),             # The values of the parameter to iterate over
                                  n_trial=5,                  # Number of trials per simulation
                                  n=N, p=P, a=A, sigma=SIGMA, rho=RHO, snr_strength_vector=SNR_STRENGTH_VECTOR, lambda_type=LAMBDA_TYPE, k=K, kkt_count=KKT_COUNT, pval = FALSE,
                                  n_surrogates=N_SURROGATES)
} else {
  
  load(here("data/simulation_results", "simulation_output_SNR_sweep_5.RData"))
  output_sim <- sim_output$runs 

}

output <- output_sim$final_results

output <- output %>%
  select(-c("Precision_Baseline_75"))



# Generate the plots
sensitivity_plot <- plot_metrics_by_param_snr(output, "n", "Sensitivity", log_scale = TRUE, custom_x_labels=c("1000"="0.1","1001"="0.3", "1002"="0.6","1003"= "1.2", "1004"="2.4"))  + theme(legend.position = "none")  # Remove legend
precision_plot <- plot_metrics_by_param_snr(output, "n", "Precision", log_scale = TRUE, custom_x_labels=c("1000"="0.1","1001"="0.3", "1002"="0.6","1003"= "1.2", "1004"="2.4"))  + theme(legend.position = "none")

# Combine the plots side by side
combined_plot <- sensitivity_plot + precision_plot

# Display the combined plot
print(combined_plot)
beep(2)

ggsave(here("data/plots/combined_plot_snr_5.png"), plot = combined_plot, width = 8, height = 8, dpi=600)


if (SAVE){
  sim_output <- list()
  sim_output$runs <- output_sim
  save(sim_output, file = here("data/simulation_results", "simulation_output_SNR_sweep_5.RData"))
  
}


end_time <- proc.time()

execution_time <- end_time - start_time
print(execution_time)


# estimated SNR 
snr1 <- output_sim$sim_data[[1]][[1]]$snr
snr2 <- output_sim$sim_data[[2]][[1]]$snr
snr3 <- output_sim$sim_data[[3]][[1]]$snr
snr4 <- output_sim$sim_data[[4]][[1]]$snr
snr5 <- output_sim$sim_data[[5]][[1]]$snr
snrs <- round(c(snr1, snr2, snr3, snr4, snr5),2) 
snrs
# 0.1, 0.3, 0.7, 1.2, 2.6


