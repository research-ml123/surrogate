library(here)
library(glmnet)
library(parallel)
library(knockoff)

source(here("src", "generate_data.R"))
source(here("src", 'surrogate_lasso.R'))
source(here("src", "baseline.R"))
source(here("src", "old/helper.R"))
source(here("src", "plotting.R"))


simulation_param_grid <- function(
    param_name,               # The name of the parameter to vary
    param_values,             # The values of the parameter to iterate over
    n_trial,                  # Number of trials per simulation
    n, p, a, sigma, rho, snr_strength_vector, lambda_type, k, kkt_count, pval = FALSE, 
    n_surrogates
) {
  # Initialize a list to store results for each parameter value
  all_results <- list()
  sim_data <- list()
  
  # Loop over the parameter values
  for (i in seq_along(param_values)) {
    param_value <- param_values[i]
    cat("Running simulations for", param_name, "=", param_value, "(", i, "of", length(param_values), ")\n")
    
    # Update the parameter dynamically
    if (param_name == "n") {
      n <- param_value
    } else if (param_name == "p") {
      p <- param_value
    } else if (param_name == "a") {
      a <- param_value
    } else if (param_name == "sigma") {
      sigma <- param_value
    } else if (param_name == "rho") {
      rho <- param_value
    } else if (param_name == "snr_strength") {
      snr_strength <- param_value
    } else if (param_name == "lambda_type") {
      lambda_type <- param_value
    } else if (param_name == "k") {
      k <- param_value
    } else if (param_name == "kkt_count") {
      kkt_count <- param_value
    } else if (param_name == "n_surrogates") {
      n_surrogates <- as.data.frame(param_value)  # Ensure it's a data frame
    } else {
      stop("Unknown parameter name: ", param_name)
    }
    
    # Run simulations for the current parameter value
    results <- simulation_n_trial(n_trial, n, p, a, n_surrogates, sigma, rho, snr_strength_vector[i], lambda_type, k, kkt_count, pval)
    
    # Add an additional column to store the parameter value
    results$results[[param_name]] <- param_value
    sim_data[[i]] <- results$sim_data
    
    # Store the results
    all_results[[i]] <- results$results
  }
  
  # Combine all results into a single data frame
  final_results <- do.call(rbind, all_results)
  final_output <- list()
  final_output$final_results <- final_results
  final_output$sim_data <- sim_data
  
  return(final_output)
}

 
simulation_n_trial <- function(n_trial, n, p, a, n_surrogates, sigma, rho, snr_strength, lambda_type = "lambda.min", k, kkt_count, pval = FALSE) {
  # Number of available cores
  num_cores <- detectCores() - 1  # Reserve one core for other processes
  cat("Using", num_cores, "cores for parallel processing.\n")

  # Parallel processing for trials
  results_list <- mclapply(1:n_trial, function(trial) {
    cat("Running trial", trial, "of", n_trial, "\n")

    # Call the single trial function
    output <- simulation_single_trial(
      n = n, p = p, a = a, n_surrogates = n_surrogates,
      sigma = sigma, rho = rho, snr_strength = snr_strength,
      lambda_type = lambda_type, k = k, kkt_count = kkt_count,
      pval = pval
    )

    # Assume `compute_metrics` in `simulation_single_trial` returns a list with sensitivity and precision
    list(
      Trial = trial,
      Sensitivity_CV = output$sensitivity_cv,
      Precision_CV = output$ppv_cv,
      Sensitivity_Baseline_25 = output$sensitivity_baseline_25,
      Precision_Baseline_25 = output$ppv_baseline_25,
      Sensitivity_Baseline_50 = output$sensitivity_baseline_50,
      Precision_Baseline_50 = output$ppv_baseline_50,
      Sensitivity_Baseline_75 = output$sensitivity_baseline_75,
      Precision_Baseline_75 = output$ppv_baseline_75,
      Sim_Data = output$sim_data
    )
  }, mc.cores = num_cores)  # Specify the number of cores to use

  # Combine results into a single data frame
  results_df <- do.call(rbind, lapply(results_list, function(res) {
    data.frame(
      Trial = res$Trial,
      Sensitivity_CV = res$Sensitivity_CV,
      Precision_CV = res$Precision_CV,
      Sensitivity_Baseline_25 = res$Sensitivity_Baseline_25,
      Precision_Baseline_25 = res$Precision_Baseline_25,
      Sensitivity_Baseline_50 = res$Sensitivity_Baseline_50,
      Precision_Baseline_50 = res$Precision_Baseline_50,
      Sensitivity_Baseline_75 = res$Sensitivity_Baseline_75,
      Precision_Baseline_75 = res$Precision_Baseline_75,
      stringsAsFactors = FALSE
    )
  }))

  # Combine simulation data
  sim_data <- lapply(results_list, function(res) res$Sim_Data)

  # Final output
  list(
    results = results_df,
    sim_data = sim_data
  )
}


simulation_single_trial <- function(n, p, a, n_surrogates, sigma, rho, snr_strength, lambda_type="lambda.min", k, kkt_count, pval=FALSE){
  
  
  beta <- rep(0, p)  # use beta to control if surrogate is in Y or not
  beta[1:a] <- sigma * snr_strength / sqrt(n) 
  total_n_surrogates <- sum(n_surrogates)
  beta_surrogate <- rep(0, total_n_surrogates) # assign beta to surrogates
  
  
  sim_data <- generate_data(n=n,
                            p=p,
                            a=a,
                            n_surrogates=n_surrogates,
                            rho=rho,
                            sigma=sigma,
                            beta=beta,
                            beta_surrogate=beta_surrogate)
  
  # Unpack simulated data
  X <- sim_data$X
  y <- sim_data$y
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  surrogates <- sim_data$surrogates
  BETA <- sim_data$params$beta
  print(paste0("SNR: ", round(sim_data$snr,2)))
  
  # check correlation for X1
  for (i in 1:n_surrogates[1]){
    print(paste0("cor surrog: ", round(cor(X[,1], X[,surrogates$`1`[i]]), 2) ))
  }
  
  
  #########################
  #   Include knockoffs   #
  #########################
  
  
  if (n > p){
    knockoffs <- create.fixed(X, randomize=F)
    X <- cbind(X, knockoffs$Xk)

    # knockoffs <- create.gaussian(X, mu=colMeans(X), Sigma=cov(X), method="equi")
    # X <- cbind(X, knockoffs)
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  } else {

    # knockoffs <- create.gaussian(X, mu=colMeans(X), Sigma=cov(X), method="equi")
    # X <- cbind(X, knockoffs)
    # colnames(X) <- paste0("X", seq_len(ncol(X)))

  }

  #cor(X[,1],knockoffs$Xk[,1])
  # cor(X[,1],knockoffs[,1])

  
  ######################### 
  # Surrogate search - CV #
  ######################### 
  
  # Generate fold id for cv
  foldid <- generate_folds(data=X, k=k)
  
  # Find lambda0 through cross-validation
  cvfit0 <- cv.glmnet(x=X,
                      y=y,
                      family="gaussian",
                      alpha=1,
                      foldid=foldid,
                      standardize=FALSE)
  # plot(cvfit0)
  
  # Fix lambda0 for all the following steps
  if (lambda_type == "lambda.1se"){
    lambda0 <- cvfit0$lambda.1se
  }
  
  if (lambda_type == "lambda.min"){
    lambda0 <- cvfit0$lambda.min
  }
  
  lambdas <- cvfit0$lambda
  
  # Get nonzero beta indices
  betas <- coef(cvfit0, s=lambda0)[-1]
  idx_nonzero <- which(abs(betas) > 0)
  
  # Define feature sets
  active_set <- colnames(X)[idx_nonzero]
  candidate_set <- setdiff(colnames(X), active_set)
  
  ###########################
  #   Run surrogate search  #
  ###########################
  
  # hclsut apply if correlation is past 90? 
  
  surrogate_result <- surrogate_search(target_set=active_set,   # active_set if you want all
                                       active_set=active_set,
                                       candidate_set=candidate_set,
                                       x=X,
                                       y=y,
                                       lambdas=lambdas,
                                       lambda0=lambda0,
                                       alpha=1,
                                       k=k,
                                       foldid= NULL,#foldid,#NULL, <<<
                                       family="gaussian",
                                       verbose=FALSE,
                                       kkt=TRUE,
                                       kkt_eps=NULL,
                                       kkt_count=kkt_count,
                                       hclust=FALSE, 
                                       hclust_rho=NULL)
  
  surrogate_list <- get_surrogates(surrogate_result)
  
  ######################################
  #   Run baseline correlation search  #
  ######################################
  baseline_result_75 <- baseline(active_set=active_set, x=X, threshold=0.75)
  baseline_result_50 <- baseline(active_set=active_set, x=X, threshold=0.5)
  baseline_result_25 <- baseline(active_set=active_set, x=X, threshold=0.25)
  
  #######################
  # SELECTIVE INFERENCE #
  #######################
  
  if (pval){
    source(here("src", "selective_inference.R"))
    
    # get active set feature signs
    signs <- as.vector(sign(betas[idx_nonzero]))
    
    # estimate Sigma
    # sigma <- sigma # error_variance_estimate(x=X, y=y, beta=lasso_coefs)
    
    # Correct path to the Python executable in the Conda environment
    temp_dir_path <- "/Users/minwoosun/Desktop/temp" #<- need to delete for new run
    python_path <- "/Users/minwoosun/mambaforge/envs/pyg/bin/python"
    
    if (!dir.exists(temp_dir_path)) {
      dir.create(temp_dir_path, showWarnings = FALSE, recursive = TRUE)
    }
    
    surrogates_pval <- get_pval_selected_surrogates(x=X,
                                                    y=y,
                                                    active=active_set,
                                                    signs=signs,
                                                    lambda=lambda0,
                                                    sigma=sigma,
                                                    temp_dir_path=temp_dir_path,
                                                    python_path=python_path,
                                                    correction_method="bonferroni",
                                                    alpha=0.05)
    surrogates_pval
  }
  
  #######################
  #      Evaluation     #
  #######################
  
  # <- need to remove activeset overlap
  ground_truth <- convert_list(surrogates)
  ground_truth <- lapply( ground_truth, function(x) x[!x %in% active_set])
  ground_truth <- flatten_list_to_vector(ground_truth)

  
  
  
  
  # print(paste0("ground truth count: ", length(ground_truth)))
  
  pred_cv <- flatten_list_to_vector(surrogate_list)
  pred_base_25 <- flatten_list_to_vector(baseline_result_25)
  pred_base_50 <- flatten_list_to_vector(baseline_result_50)
  pred_base_75 <- flatten_list_to_vector(baseline_result_75)
  
  metrics_baseline_25 <- compute_metrics(ground_truth, pred_base_25, candidate_set = candidate_set, surrogates = surrogates)
  metrics_baseline_50 <- compute_metrics(ground_truth, pred_base_50, candidate_set = candidate_set, surrogates = surrogates)
  metrics_baseline_75 <- compute_metrics(ground_truth, pred_base_75, candidate_set = candidate_set, surrogates = surrogates)
  metrics_cv <- compute_metrics(ground_truth, pred_cv, candidate_set = candidate_set, surrogates = surrogates)

  
  return(list(
    sensitivity_cv = metrics_cv$sensitivity,
    ppv_cv = metrics_cv$precision,
    
    sensitivity_baseline_25 = metrics_baseline_25$sensitivity,
    ppv_baseline_25 = metrics_baseline_25$precision,
    sensitivity_baseline_50 = metrics_baseline_50$sensitivity,
    ppv_baseline_50 = metrics_baseline_50$precision,
    sensitivity_baseline_75 = metrics_baseline_75$sensitivity,
    ppv_baseline_75 = metrics_baseline_75$precision,
    
    # specificity_cv = metrics_cv$specificity,
    ground_truth = ground_truth,
    selected_features = pred_cv,
    sim_data = sim_data,
    active_set = active_set
  ))
  
}




