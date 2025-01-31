
#' Generate Synthetic Data with Surrogate Features
#'
#' This function creates a synthetic dataset with `n` observations and `p` features.
#' Among these features, `a` are active (i.e., contribute to the response), and each
#' active feature has a specified number of surrogates correlated with it.
#'
#' @param n Integer. Number of observations.
#' @param p Integer. Total number of features (columns in X).
#' @param a Integer. Number of active features.
#' @param n_surrogates Integer vector of length \code{a}, specifying how many surrogates each active feature has.
#' @param rho Numeric. Correlation coefficient used when creating surrogates.
#' @param sigma Numeric. Standard deviation of the noise term added to the response.
#' @param beta Numeric vector of length \code{p}, coefficients for generating the response \code{y}.
#'
#' @details
#' \enumerate{
#'   \item \strong{Surrogates}: For each active feature, the function samples as many surrogate features
#'         as specified in \code{n_surrogates}, making them correlated with the original active feature.
#'   \item \strong{Response \code{y}}: Generated as \eqn{X \beta + \sigma \epsilon}, where \eqn{\epsilon}
#'         is standard normal noise. The resulting \code{y} is centered by subtracting its mean.
#' }
#' 
#' #'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\code{X}}{An \code{n x p} matrix of standardized features.}
#'   \item{\code{y}}{A numeric vector of length \code{n}, the generated response.}
#'   \item{\code{surrogates}}{A list indicating which features are surrogates of each active feature.}
#'   \item{\code{snr}}{The signal-to-noise ratio of the generated data.}
#'   \item{\code{params}}{A list of input parameters used to generate the data.}
#' }
#'
#' @export
generate_data <- function(n, p, a, n_surrogates, rho, sigma, beta, beta_surrogate){
  
  ##########
  # checks #
  ##########
  if (length(n_surrogates) != a){
    stop("n_surrogates must have 'a' elements.")
  }
  
  if (length(beta_surrogate) != sum(n_surrogates)){
    stop("The length of beta_surrogate has to match the total number of surrogates.")
  }
  
  ##############
  # Generate X #
  ##############
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  
  #####################
  # Create surrogates #
  #####################
  
  # # NOT ACCOUNTING FOR OVERLAP
  # surrogates <- list()
  # 
  # for (i in 1:a){  # for each active feature
  #   for (j in 1:n_surrogates[i]){  # assign n_surrog[i] surrogates
  #     s <- sample((a + 1):p, 1)
  #     X[, s] <- X[, i] * rho + rnorm(n)
  #     surrogates <- append(surrogates, list(c(i, s)))
  #   }
  # }
  # 
  surrogates <- list()
  available_indices <- setdiff((a + 1):p, 1:a)  # Indices available for surrogates
  
  rho_tally <- 1
  
  for (i in 1:a) {  # for each active feature
    if(n_surrogates[i] > 0){
      for (j in 1:n_surrogates[i]) {  # assign n_surrog[i] surrogates
        
        if (length(available_indices) == 0) {
          stop("Not enough available indices to assign surrogates. Check p and n_surrogates.")
        }
        
        s <- sample(available_indices, 1)  # Sample from available indices
        X[, s] <- X[, i] * rho[rho_tally] + rnorm(n)
        surrogates <- append(surrogates, list(c(i, s)))
        available_indices <- setdiff(available_indices, s)  # Remove assigned index
        rho_tally <- rho_tally + 1
      }
    } else {
      next
    }

  }
  
  
  
  # Reorganize surrogates
  df <- data.frame(do.call(rbind, surrogates))  # each row is one pair
  colnames(df) <- c("first", "second")
  surrogates <- split(df$second, df$first)
  
  # update beta for surrogates (default is 0)
  all_surrogate_idx <- unlist(surrogates)
  beta[all_surrogate_idx] <- beta_surrogate
  
  ########
  # Data #
  ########
  # Standardize X
  X <- scale(X)
  
  # Generate y
  y <- X %*% beta + sigma * rnorm(n)
  y <- y - mean(y)
  
  # Calculate signal-to-noise ratio (SNR)
  snr <- sum((X %*% beta)^2) / (n * sigma^2)
  
  ##########
  # Output #
  ##########
  params = list()
  params$n = n
  params$p = p
  params$a = a
  params$n_surrogates = n_surrogates
  params$beta = beta
  
  output = list()
  output$X = X
  output$y = y
  output$surrogates = surrogates
  output$snr = snr
  output$params = params
  
  return(output)
}







