#' Generate Fold IDs for k-fold Cross Validation
#'
#' This function assigns fold IDs to each row of a data matrix or data frame 
#' for the purpose of k-fold cross validation. Rows are randomly shuffled 
#' before assigning fold IDs. In cases where the total number of rows isn't 
#' perfectly divisible by `k`, extra rows are randomly distributed among the folds.
#'
#' @param data A matrix or data frame for which fold IDs need to be generated.
#' @param k An integer specifying the number of folds for cross validation.
#'
#' @return A numeric vector of fold IDs with the same length as the number of rows in `data`.
#'
#' @examples
#' \dontrun{
#' data <- matrix(1:100, ncol=2)
#' k <- 5
#' fold_ids <- generate_folds(data, k)
#' table(fold_ids)
#' }
#' 
generate_folds <- function(data, k) {
  # Randomly shuffle the row indices
  shuffled_indices <- sample(1:nrow(data))
  
  # Calculate the size of each fold
  fold_sizes <- rep(floor(nrow(data) / k), k)
  
  # Distribute the remainder to random folds
  remainder <- nrow(data) %% k
  if (remainder > 0) {
    fold_sizes[sample(1:k, remainder, replace = FALSE)] <- fold_sizes[sample(1:k, remainder, replace = FALSE)] + 1
  }
  
  # Assign fold IDs
  fold_ids <- rep(NA, nrow(data))
  start_idx <- 1
  for (i in 1:k) {
    end_idx <- start_idx + fold_sizes[i] - 1
    fold_ids[shuffled_indices[start_idx:end_idx]] <- i
    start_idx <- end_idx + 1
  }
  
  return(fold_ids)
}



convert_list <- function(input_list) {
  output_list <- list()
  
  for (i in seq_along(input_list)) {
    # Prefix all elements with "X" regardless of the index
    output_list[[paste0("X", i)]] <- paste0("X", input_list[[i]])
  }
  
  return(output_list)
}


flatten_list_to_vector <- function(input_list) {
  # Initialize an empty vector
  result <- c()
  
  # Iterate through each element in the list
  for (name in names(input_list)) {
    # Combine the list name with each value in the vector
    result <- c(result, paste0(name, input_list[[name]]))
  }
  
  return(result)
}



compute_metrics <- function(ground_truth, pred_cv, candidate_set = NULL, surrogates = NULL) {
  # Compute confusion metrics
  ground_truth_set <- unique(ground_truth)
  pred_cv_set <- unique(pred_cv)
  
  TP <- length(intersect(ground_truth_set, pred_cv_set))
  FP <- length(setdiff(pred_cv_set, ground_truth_set))
  FN <- length(setdiff(ground_truth_set, pred_cv_set))
  
  # Compute TN if a candidate set is provided
  if (!is.null(candidate_set)) {
    total_set <- unique(candidate_set)
    negatives <- setdiff(total_set, union(ground_truth_set, pred_cv_set))
    TN <- length(negatives)
  } else {
    TN <- NA
  }
  
  # True Positive Rate (TPR) and False Positive Rate (FPR)
  TPR <- ifelse((TP + FN) > 0, TP / (TP + FN), NA) # Recall
  FPR <- ifelse((FP + TN) > 0, FP / (FP + TN), NA) # Fall-out
  
  # Positive Predictive Value (PPV) - Precision
  PPV <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
  
  # Return all metrics as a list
  list(
    # TP = TP,
    # FP = FP,
    # FN = FN,
    # TN = TN,
    sensitivity = TPR,    # Recall
 #   specificity = 1 - FPR, # True Negative Rate
    precision = PPV       # Positive Predictive Value
  )
}

