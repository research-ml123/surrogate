library(here)
library(glmnet)
library(proxy)
library(pROC)


source(here('src', 'helper.R'))


#' Create Feature Sets for Surrogate Modeling
#'
#' This function generates a list of feature sets for building surrogate models,
#' where one or more features from a candidate set are used to approximate a
#' target feature from an active set.
#'
#' @param active_set A vector of feature names in the active set.
#' @param candidate_set A vector of feature names in the candidate set.
#' @param target_feature The name of the target feature to be approximated.
#' @param drop_set A vector of feature names to be excluded from the active set.
#'
#' @return A list of feature sets, with keys:
#'   \describe{
#'     \item{__base__}{The original active set, including the target feature.}
#'     \item{__null__}{The active set without the target feature or features in drop_set.}
#'     \item{<feature_name>}{The active set without the target feature or features
#'       in drop_set, but with the corresponding feature from the candidate set.}
#'   }
#'
#' @examples
#' active_set <- c("A", "B", "C", "D")
#' candidate_set <- c("X", "Y", "Z")
#' target_feature <- "C"
#' drop_set <- "D"
#' surrogate_sets <- create_surrogate_set(active_set, candidate_set, target_feature, drop_set)
#' str(surrogate_sets)
#'
#' @export
create_surrogate_set <- function(active_set, candidate_set, target_feature, drop_set){
  
  if (target_feature %in% candidate_set){
    stop("Input target_feature cannot be in candidate_set")
  }
  
  if (!(target_feature %in% active_set)){
    stop("Input target_feature is missing in active_set")
  }
  
  # always contain target_feature in dropset
  drop_set <- c(drop_set, target_feature)
  
  surrogate_sets <- list()
  
  # include model with target_feature (base)
  surrogate_sets[['__base__']] <- active_set
  
  # selected features as surrogate (null)
  # if this wins, better to not have surrogate
  surrogate_sets[['__null__']] <- active_set[!(active_set %in% drop_set)]
  
  # remove drop_set features from active_set
  active_set <- active_set[!(active_set %in% drop_set)]
  
  # excluded features as surrogate
  for (i in 1:length(candidate_set)){
    surrogate <- candidate_set[i]
    features <- c(surrogate, active_set)
    surrogate_sets[[surrogate]] <- features 
  }
  
  return(surrogate_sets)
}


#' Perform Hierarchical Clustering on Variables
#'
#' The `clusterList` function performs hierarchical clustering on a set of variables
#' and returns a list containing the feature clusters. The function computes the
#' distances between variables using correlation or Jaccard distance, depending on
#' whether the data is binary or not. The height of the dendrogram cut is determined
#' by the specified `rho` (correlation) value.
#'
#' @param x A numeric data frame or matrix containing the variables to be clustered.
#'   Each row represents an observation, and each column represents a variable.
#' @param rho A numeric value between 0 and 1, specifying the desired correlation
#'   threshold for cutting the dendrogram. Variables with a correlation greater than
#'   or equal to `rho` will be grouped into the same cluster.
#' @param plotting A logical value indicating whether to plot the dendrogram with
#'   the cut height marked by a red line. Default is `FALSE`.
#' @param binary A logical value indicating whether the data is binary (TRUE) or
#'   not (FALSE). If `TRUE`, the Jaccard distance is used for calculating the
#'   distances between variables; otherwise, the correlation-based distance is used.
#'
#' @return A list, where each element represents a cluster of features (variables).
#'   The names of the list elements correspond to the feature names, and the values
#'   are the lists of features in the respective clusters.
#'
#' @examples
#' # Load the iris dataset
#' data(iris)
#'
#' # Convert the iris dataset to a numeric matrix
#' iris_numeric <- as.matrix(iris[, 1:4])
#'
#' # Cluster the variables with a correlation threshold of 0.7
#' clusters <- clusterList(iris_numeric, rho = 0.7, plotting = TRUE)
#'
#' # Print the cluster assignments
#' print(clusters)
#'
#' @export
clusterList <- function(x, rho, plotting=FALSE, binary=FALSE){
  
  # Use correlations between variables "as distance"
  if (binary){
    transposed_data <- t(x)
    dd <- proxy::dist(transposed_data, method = "Jaccard")
    cut_height <- 1 - rho
  } else {
    # dd <- as.dist((1 - cor(x))/2)
    dd <- as.dist((1-abs(cor(x))))
    cut_height <- 1 - rho
  }
  
  hclust_obj <- hclust(dd)
  
  # Cut the dendrogram at the specified height
  clusters <- cutree(hclust_obj, h = cut_height)
  
  # Create an empty list to store the result
  cluster_list <- list()
  
  # Loop through each unique cluster
  for (cluster_number in unique(clusters)) {
    
    # Find the features that belong to this cluster
    features_in_cluster <- names(clusters[clusters == cluster_number])
    
    # For each feature in the cluster, assign the cluster's feature list to the result list
    for (feature in features_in_cluster) {
      cluster_list[[feature]] <- features_in_cluster
    }
  } 
  
  # plot cluster
  if(plotting){
    plot(hclust_obj)
    abline(h=cut_height, col='red')
  }
  
  return(cluster_list)
}


#' Extract Surrogate Feature Names
#'
#' This function extracts the surrogate feature names from a vector of feature names.
#' The vector should contain the strings "__base__" and "__null__" in a specific order,
#' and the surrogate feature names should be prefixed with "X_".
#'
#' @param vec A vector of feature names, including "__base__", "__null__", and
#'   surrogate feature names prefixed with "X_".
#'
#' @return A vector of surrogate feature names, or NULL if the required conditions
#'   are not met.
#'
#' @examples
#' feature_names <- c("X_feature1", "__base__", "X_feature2", "__null__", "X_feature3")
#' extract_surrogates(feature_names)
#'
#' @export
extract_surrogates <- function(vec) {
  base_index <- which(vec == "__base__")
  null_index <- which(vec == "__null__")
  
  # Check if "__base__" comes before "__null__"
  if (base_index < null_index) {
    # Extract "X_" elements before "__base__" and between "__base__" and "__null__"
    xs_before_base <- vec[1:(base_index-1)]
    xs_between_base_null <- vec[(base_index+1):(null_index-1)]
    out <- c(xs_before_base, xs_between_base_null)
    out <- out[!(out == "__base__")]
    if ((length(out)==1) & any(out %in% "__null__")){out <- NULL}
    return(out)
  } else {
    return(NULL) # Condition not met
  }
}


#' Extract Surrogate Features from Search Results
#'
#' This function processes a list of search results containing potential surrogate
#' features. Each element in the list is expected to be a list that includes a
#' sub-list or vector under the "surrogate" key, from which surrogate features are extracted.
#' The function utilizes `extract_surrogates` internally to parse each vector and extract
#' feature names that are prefixed with "X_" and situated between "__base__" and "__null__" markers.
#'
#' @param search_result Direct output from the function surrogate_search().
#'
#' @return A list where each element is named after a sequential "X" prefix followed by an
#'   integer (starting from 1), and contains a vector of extracted surrogate feature names
#'   for the corresponding search result. If no valid surrogates are found, the corresponding
#'   list element will be NULL.
#'
#' @export
get_surrogates <- function(search_result){
    
  result <- search_result[["search_result"]]
  target_names <- names(result)
  surrogate_list <- list()
  
  for (i in 1:length(result)){
    surrogs <- extract_surrogates(result[[i]][["surrogate"]])
    target_name <- target_names[i]
    surrogate_list[[target_name]] <- surrogs
  }
  
  # remove __base__
  # surrogate_list <- lapply(surrogate_list, function(x) x[x != "__base__"])
  
  return(surrogate_list)
}


#' Create a Contingency Table for Feature Selection
#'
#' This function constructs a contingency table based on true features, selected features,
#' active features, and candidate features. It helps in evaluating the effectiveness
#' of the feature selection process by calculating true positives, false positives,
#' true negatives, and false negatives.
#'
#' @param true A vector of true features.
#' @param selected A vector of selected features, which are filtered to exclude '__null__' and '__base__'.
#' @param active A vector of currently active features, which are excluded from the candidate set.
#' @param candidate A vector of potential candidate features.
#' @param verbose A boolean flag that, when TRUE, prints detailed information about the feature sets.
#'
#' @return A list containing counts of true positives (TP), true negatives (TN),
#'         false positives (FP), and false negatives (FN).
#'
#' @details
#' The function first filters out non-candidate and active features from the selected and candidate lists, respectively.
#' It then identifies true positives, false positives, true negatives, and false negatives based on the
#' intersections and differences between the true, selected, and not selected sets. The function also
#' checks if the total counts from the contingency table equal the length of the candidate set, throwing
#' an error if they do not match.
#'
#' @examples
#' true0 <- c("X4", "X5", "X6", "X50")
#' active0 <- c("X1", "X2,", "X3", "X50")
#' selected0 <- c("X4", "X5", "X100")
#' candidate0 <- paste0("X", 4:100)
#'
#' table0 <- contingency_table(true0, selected0, active0, candidate0, verbose = TRUE)
#' print(table0)
#'
#' @export
contingency_table <- function(true, selected, active, candidate, verbose=FALSE){
  
  # remove noncandidates
  selected <- selected[!(selected %in% c('__null__', '__base__'))]
  
  # remove active features in candidate set
  candidate <- candidate[!c(candidate %in% active)]
  
  # remove true surrogates in active set
  true <- true[!(true %in% active)]
  
  # define set of candidates that were not selected
  not_selected <- candidate[!(candidate %in% selected)]
  
  TP <- selected[selected %in% true]
  FP <- selected[!(selected %in% true)]
  TN <- not_selected[!(not_selected %in% true)]
  FN <- true[!(true %in% selected)]
  
  TP_count <- length(TP)
  FP_count <- length(FP)
  TN_count <- length(TN)
  FN_count <- length(FN)
  
  contingency <- list(TP=TP_count,
                      TN=TN_count,
                      FP=FP_count,
                      FN=FN_count)
  
  # total sum of the contingency table should equal to the length of the candidate set
  if (sum(TP_count, TN_count, FP_count, FN_count) != length(candidate)){
    stop("Total sum of the contingency table != length of candidate set.")
  }
  
  if (verbose){
    print(paste0("True set: ", length(true)))
    print(paste0("Active set: ", length(active)))
    print(paste0("Candidate set: ", length(candidate)))
  }
  
  return(contingency)
}
