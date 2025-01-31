
compute_pvalues <- function(x, y, active, signs, lambda, sigma, temp_dir_path, python_path){
  
  # Save data to CSV files in the temporary directory
  write.csv(x, file.path(temp_dir_path, "X.csv"), row.names = FALSE)
  write.csv(data.frame(y = y), file.path(temp_dir_path, "y.csv"), row.names = FALSE)
  write.csv(data.frame(active_set = active), file.path(temp_dir_path, "active_set.csv"), row.names = FALSE)
  write.csv(data.frame(signs = signs), file.path(temp_dir_path, "signs.csv"), row.names = FALSE)
  writeLines(as.character(lambda), file.path(temp_dir_path, "lambda0.txt"))
  writeLines(as.character(sigma), file.path(temp_dir_path, "sigma.txt"))
  
  # Run python script with data from temp dir
  python_script_path <- here("src", "pvalue_og.py")
  command <- paste(shQuote(python_path), shQuote(python_script_path), shQuote(temp_dir_path))
  system_output <- system(command, intern = TRUE)
  
  # Load result back into R
  pvals_path <- paste0(temp_dir_path, "/pvals.csv")
  pvals_df <- read.csv(pvals_path, stringsAsFactors = FALSE)
  # names(pvals_df) <- paste0("X", 1:ncol(pvals_df))
  
  return(pvals_df)
}


get_pval_selected_surrogates <- function(x, y, active, signs, lambda, sigma, temp_dir_path, python_path, correction_method="BH", alpha=0.05){
  
  # convert active_set into python zero index
  active_index <- which(colnames(x) %in% active_set) -1 
  
  pvals <- compute_pvalues(x=x,
                           y=y,
                           active=active_index,
                           signs=signs,
                           lambda=lambda,
                           sigma=sigma,
                           temp_dir_path=temp_dir_path,
                           python_path=python_path)
  
  corrected_pval <- adjust_pval_by_column(pvals, method=correction_method) 
  
  total_selected <- sum(corrected_pval < alpha, na.rm=TRUE)
  
  print(total_selected)
  # print(paste0("Total selected: ", total_selected))
  indices <- which(corrected_pval < alpha, arr.ind = TRUE)
  indices <- cbind(indices, corrected_pval[indices])
  result <- data.frame(indices)
  names(result) <- c("surrogate", "target", "pval")
  
  # convert indices back to feature name
  result$surrogate <- colnames(x)[result$surrogate]
  result$target <- active[result$target]
  
  return(result)
}


create_surrogate_list <- function(surrogate_df){
  
  # convert dataframe to organized list
  selected_list <- split(surrogate_df$surrogate, surrogate_df$target)
  
  # # add in missing 
  # indices <- names(surrogate_list)
  # selected_list <- ensure_indices_sorted(selected_list, indices)
  
  return(selected_list)
}


## estimating error variance
# from: A Study of Error Variance Estimation in Lasso Regression
error_variance_estimate <- function(x, y, beta){
  n <- nrow(x)
  s <- sum(abs(beta) > 0)
  sigma <- (1 / (n-s)) * t(y - x %*% beta) %*% (y - x %*% beta)
  return(sigma)
}




adjust_pval <- function(pval_df, method){
  
  # Flatten the dataframe to a vector
  data_vector <- as.vector(as.matrix(pval_df))
  
  # apply correction
  adjusted_pvals <- p.adjust(data_vector, method = method)
  
  # Reshape the ranked vector back into the original dataframe's dimensions
  adjusted_pvals <- matrix(adjusted_pvals, nrow = nrow(pval_df), byrow = FALSE)
  adjusted_pvals <- as.data.frame(adjusted_pvals)
  
  return(adjusted_pvals)
}


adjust_pval_by_column <- function(pval_df, method){
  
  adjusted_columns <- list()
  
  for (j in 1:ncol(pval_df)){
    
    data_vector <- as.vector(as.matrix(pval_df[,j]))
    adjusted_pvals <- p.adjust(data_vector, method = method)
    adjusted_columns[[j]] <- adjusted_pvals
  }
  
  adjusted_pvals <- data.frame(adjusted_columns)
  
  return(adjusted_pvals)
}

