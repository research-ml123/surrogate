source(here("src", "simulation_delta.R"))

SAVE <- TRUE
LOAD <- FALSE

if (!LOAD){
  n <- 5000
  p <- 100
  a <- 3  # number of active features
  n_surrogates <- c(5,5,5)
  sigma <- 2
  rho <- 1.8 # control surrogate correlation 0.6 to 1.8
  snr_strength <- 50 # 12
  k <- 10
  kkt_count <- 20 # at kkt 5 max should be at0.556 <- something is wrong
  lambda_type <- "lambda.min" 
  # 
  
  n_trial <- 5
  rhos <- seq(from=0, to=1.5, by=0.1)
  snr_strength_vector <- rep(50, length(rhos))
  out <- simulation_param_grid(param_name="rho",               # The name of the parameter to vary
                                param_values=rhos,             # The values of the parameter to iterate over
                                n_trial,                  # Number of trials per simulation
                                n, p, a, sigma, rho, snr_strength_vector, lambda_type, k, kkt_count, 
                                n_surrogates)
  if (SAVE){
    save(out, file = here("data/simulation_results", "simulation_output_rho_sweep_5.RData"))
  }
  
} else {
  load(here("data/simulation_results", "simulation_output_rho_sweep_5.RData"))
}



prop_plot <- plot_metrics_by_param_rho_sweep(out, "rho", "Prop", log_scale = TRUE)  + theme(legend.position = "none")
# perc_plot <- plot_metrics_by_param(out, "rho", "Perc", log_scale = TRUE)
# 
# combined_plot <- prop_plot + perc_plot
# 
# # Display the combined plot
# print(combined_plot)

prop_plot 
ggsave(here("data/plots/proportion_plot_5.png"), plot = prop_plot, width = 8, height = 8, dpi=600)
beep(2)

