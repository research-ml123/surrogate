library(here)
library(glmnet)
library(parallel)

source(here("src", "generate_data.R"))
source(here("src", 'surrogate_lasso.R'))
source(here("src", "baseline.R"))
source(here("src", "old/helper.R"))
source(here("src", "plotting.R"))


simulation_param_grid <- function(
    param_name,               # The name of the parameter to vary
    param_values,             # The values of the parameter to iterate over
    n_trial,                  # Number of trials per simulation
    n, p, a, sigma, rho, snr_strength_vector, lambda_type, k, kkt_count, 
    n_surrogates
) {
  # Initialize a list to store results for each parameter value
  all_results <- list()

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
    results <- simulation_n_trial(n_trial, n, p, a, n_surrogates, sigma, rho, snr_strength_vector[i], lambda_type, k, kkt_count)
    
    # Add an additional column to store the parameter value
    results$results[[param_name]] <- param_value
    
    # Store the results
    all_results[[i]] <- results$results
  }
  
  # Combine all results into a single data frame
  final_results <- do.call(rbind, all_results)

  return(final_results)
}


simulation_n_trial <- function(n_trial, n, p, a, n_surrogates, sigma, rho, snr_strength, lambda_type = "lambda.min", k, kkt_count) {

  # Number of available cores
  num_cores <- detectCores() - 1  # Reserve one core for other processes
  cat("Using", num_cores, "cores for parallel processing.\n")

  # Parallel execution using mclapply
  results_list <- mclapply(1:n_trial, function(trial) {
    cat("Running trial", trial, "of", n_trial, "\n")

    # Call the single trial function
    output <- simulation_single_trial(n = n, p = p, a = a, n_surrogates = n_surrogates,
                                      sigma = sigma, rho = rho, snr_strength = snr_strength,
                                      lambda_type = lambda_type, k = k, kkt_count = kkt_count)

    # Store results in a list
    list(
      Trial = trial,
      Prop_ss = output$prop_ss,
      Prop_25 = output$prop_25,
      Prop_50 = output$prop_50,
      Prop_75 = output$prop_75,
      Perc_ss = output$perc_ss,
      Perc_25 = output$perc_25,
      Perc_50 = output$perc_50,
      Perc_75 = output$perc_75
    )
  }, mc.cores = num_cores)  # Use multiple cores

  # Combine results into a single data frame
  results_df <- do.call(rbind, lapply(results_list, function(res) {
    data.frame(
      Trial = res$Trial,
      Prop_ss = res$Prop_ss,
      Prop_25 = res$Prop_25,
      Prop_50 = res$Prop_50,
      Prop_75 = res$Prop_75,
      Perc_ss = res$Perc_ss,
      Perc_25 = res$Perc_25,
      Perc_50 = res$Perc_50,
      Perc_75 = res$Perc_75,
      stringsAsFactors = FALSE
    )
  }))

  # Return results
  list(results = results_df)
}




simulation_single_trial <- function(n, p, a, n_surrogates, sigma, rho, snr_strength, lambda_type="lambda.min", k, kkt_count){
  
  
  beta <- rep(0, p)  # use beta to control if surrogate is in Y or not
  beta[1:a] <- sigma * snr_strength / sqrt(n) 
  total_n_surrogates <- sum(n_surrogates)
  beta_surrogate <- rep(0, total_n_surrogates) # assign beta to surrogates
  rho <- rep(rho, total_n_surrogates)
  
  sim_data <- generate_data(n=n,
                            p=p,
                            a=a,
                            n_surrogates=n_surrogates,
                            rho=rho,
                            sigma=sigma,
                            beta=beta,
                            beta_surrogate=beta_surrogate)
  
  # <- TRAIN - TEST SPLIT
  
  # Unpack simulated data
  X <- sim_data$X
  y <- sim_data$y
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  surrogates <- sim_data$surrogates
#  BETA <- sim_data$params$beta
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
  

  ##########################
  #    train-test split    #
  ##########################
  train_proportion <- 0.8
  n_samples <- nrow(X)
  train_indices <- sample(seq_len(n_samples), size = floor(train_proportion * n_samples))
  X_train <- X[train_indices, , drop = FALSE]
  y_train <- y[train_indices]
  X_test <- X[-train_indices, , drop = FALSE]
  y_test <- y[-train_indices]
  
  
  
  
  ######################### 
  # Surrogate search - CV #
  ######################### 
  
  # Generate fold id for cv
  foldid <- generate_folds(data=X_train, k=k)
  
  # Find lambda0 through cross-validation
  cvfit0 <- cv.glmnet(x=X_train,
                      y=y_train,
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
  active_set <- colnames(X_train)[idx_nonzero]
  candidate_set <- setdiff(colnames(X_train), active_set)
  
  ###########################
  #   Run surrogate search  #
  ###########################
  
  surrogate_result <- surrogate_search(target_set=active_set,   # active_set if you want all
                                       active_set=active_set,
                                       candidate_set=candidate_set,
                                       x=X_train,
                                       y=y_train,
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
  baseline_result_75 <- baseline(active_set=active_set, x=X_train, threshold=0.75)
  baseline_result_50 <- baseline(active_set=active_set, x=X_train, threshold=0.5)
  baseline_result_25 <- baseline(active_set=active_set, x=X_train, threshold=0.25)
  
  
  ##################
  #    Test time   #
  ##################
  
  if (length(surrogate_list) > 0) {
    df_ss <- evaluate_on_test(x_train=X_train, y_train=y_train, x_test=X_test, y_test=y_test, family="gaussian", standardize=TRUE, lambda=lambda0, active_set=active_set, surrogate_list=surrogate_list )
    prop_ss <- sum(df_ss$good_surrogate) / dim(df_ss)[1]
    perc_ss <- mean( (df_ss$loss_null- df_ss$loss_surrogate) / df_ss$loss_null)    
  } else {
    prop_ss <- 0
    perc_ss <- 0
  }
  
  if (length(baseline_result_25 ) > 0) {
    df_25 <- evaluate_on_test(x_train=X_train, y_train=y_train, x_test=X_test, y_test=y_test, family="gaussian", standardize=TRUE, lambda=lambda0, active_set=active_set, surrogate_list=baseline_result_25 )
    prop_25 <- sum(df_25$good_surrogate) / dim(df_25)[1]
    perc_25 <- mean( (df_25$loss_null- df_25$loss_surrogate) / df_25$loss_null)
  } else {
    prop_25 <- 0
    perc_25 <- 0
  }
  
  if (length(baseline_result_50 ) > 0) {
    df_50 <- evaluate_on_test(x_train=X_train, y_train=y_train, x_test=X_test, y_test=y_test, family="gaussian", standardize=TRUE, lambda=lambda0, active_set=active_set, surrogate_list=baseline_result_50 )
    prop_50 <- sum(df_50$good_surrogate) / dim(df_50)[1]
    perc_50 <- mean( (df_50$loss_null- df_50$loss_surrogate) / df_50$loss_null)
  } else {
    prop_50 <- 0
    perc_50 <- 0
  }
  
  if (length(baseline_result_75) > 0) {
    df_75 <- evaluate_on_test(x_train=X_train, y_train=y_train, x_test=X_test, y_test=y_test, family="gaussian", standardize=TRUE, lambda=lambda0, active_set=active_set, surrogate_list=baseline_result_75)
    prop_75 <- sum(df_75$good_surrogate) / dim(df_75)[1]
    perc_75 <- mean( (df_75$loss_null- df_75$loss_surrogate) / df_75$loss_null)
  } else {
    prop_75 <- 0
    perc_75 <- 0
  }

  return(list(
    prop_ss = prop_ss,
    prop_25 = prop_25,
    prop_50 = prop_50,
    prop_75 = prop_75,
    
    perc_ss = perc_ss,
    perc_25 = perc_25,
    perc_50 = perc_50,
    perc_75 = perc_75
  ))
}


