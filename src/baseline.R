baseline <- function(active_set, x, threshold = 0.5) {
  # Prepare an empty list to store results
  results <- list()
  
  # Identify columns not in the active_set
  other_vars <- setdiff(colnames(x), active_set)
  
  # Loop over each variable in the active set
  for (v in active_set) {
    # Compute correlation with all non-active columns
    cor_vals <- cor(x[, v], x[, other_vars])
    
    # Choose variables whose absolute correlation is >= threshold
    keep_vars <- other_vars[abs(cor_vals) >= threshold]
    
    # Only store non-empty vectors in 'results'
    if (length(keep_vars) > 0) {
      results[[v]] <- keep_vars
    }
  }
  
  return(results)
}

