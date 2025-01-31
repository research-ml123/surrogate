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


#' Select Specific Features from a Data Frame or Matrix
#'
#' This function takes a data frame or matrix `X` and a vector of feature names `features`. 
#' It returns a subset of `X` containing only the specified features.
#'
#' @param X A data frame or matrix from which features are to be selected.
#' @param features A character vector specifying the names of the features to be selected from `X`.
#'
#' @return A subset of `X` containing only the columns with names matching the `features` vector.
#' If no features match, an empty data frame or matrix is returned.
select_columns <- function(X, features){
  return(X[, colnames(X) %in% features])
}


#' Translate Feature Names to Standard Format
#'
#' This function converts a vector of arbitrary feature names into a vector of standardized names
#' formatted as "X1", "X2", ..., "Xp", where p is the length of the input vector.
#'
#' @param features A character vector of feature names.
#'
#' @return A character vector with names formatted as "X" followed by the index of each feature.
#'
#' @examples
#' feature_names <- c("feature_1", "feature_2", "feature_3")
#' translate_to_X(feature_names)
#'
#' @export
translate_to_X <- function(features) {
  seq_along(features) %>%
    paste0("X", .)
}


#' Translate Standard Format Back to Feature Names
#'
#' This function converts a vector of standardized names ("X1", "X2", ..., "Xp") back to their
#' original arbitrary feature names based on a predefined mapping.
#'
#' @param x_names A character vector of names in the format "X1", "X2", ..., "Xp".
#' @param original_features A character vector of original feature names that correspond to "X" names.
#'
#' @return A character vector of the original feature names corresponding to the input "X" names.
#'
#' @examples
#' feature_names <- c("feature_1", "feature_2", "feature_3")
#' x_names <- translate_to_X(feature_names)
#' translate_to_features(x_names, feature_names)
#'
#' @export
translate_to_features <- function(x_names, original_features) {
  indices <- as.numeric(sub("X", "", x_names))
  original_features[indices]
}


# Function to check and add missing indices with empty vectors
ensure_indices_sorted <- function(lst, idx) {
  # Ensure all indices are present, adding empty if missing
  for (i in idx) {
    if (!i %in% names(lst)) {
      lst[[i]] <- character(0)  # Adding an empty vector
    }
  }
  # Reorder list according to idx
  lst <- lst[idx]
  return(lst)
}

########################## 
#   Evaluation Metrics   #
########################## 

log_loss <- function(y, yhat) {
  epsilon <- 1e-15 # avoid log(0) errors
  yhat <- pmin(pmax(yhat, epsilon), 1 - epsilon)
  
  ll <- -mean(y * log(yhat) + (1 - y) * log(1 - yhat))
  return(ll)
}


mean_squared_error <- function(y, yhat){
  mse <- mean((y - yhat)^2)
  return(mse)
}


# Define a function to calculate metrics for a single contingency table
calculate_metrics <- function(results) {
  TP <- results$TP
  TN <- results$TN
  FP <- results$FP
  FN <- results$FN
  
  accuracy <- (TP + TN) / (TP + TN + FP + FN) 
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  
  return(list(accuracy = accuracy, sensitivity = sensitivity, specificity = specificity))
}


##########################
#   P value correction   #
########################## 
bonferroni <- function(pval_df, n=NULL){
  
  if (is.null(n)){
    # default, correct by total number of experiments
    n = nrow(pval_df) * ncol(pval_df)
  }
  
  adjusted_pvals = pval_df * n
  return(adjusted_pvals)
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


