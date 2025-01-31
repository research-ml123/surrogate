source(here("src", 'surrogate_core.R'))
source(here("src", 'helper.R'))


kkt_screening <- function(target_feature, active_set, drop_set, x, y, lambdas, lambda0, kkt_eps=1e-3, kkt_count=NULL, family, standardize=FALSE, verbose=FALSE){
  
  if(verbose){cat("\n", paste0("-----> KKT screening for ", target_feature, " candidate surrogates..."))}
  
  # # get target feature index
  # target_idx <- which(target_feature == colnames(x))
  # 
  # # exclude active cluster set
  # cluster_set_k <- cluster_list[[target_idx]]
  # drop_set_idx <- which(active_set %in% cluster_set_k) 
  # # default would be the target feature b/c included in cluster set k
  
  drop_set_idx <- which(active_set %in% drop_set) 
  
  # refit model without the active cluster features
  fit_exclude_target <- glmnet(x=x[,-drop_set_idx], 
                               y=y, 
                               alpha=1, 
                               lambda=lambdas, 
                               family=family, 
                               standardize=standardize)
  
  betas <- coef(fit_exclude_target, s=lambda0)
  ones <- rep(1, nrow(x))
  x <- cbind(ones, x[,-drop_set_idx])
  
  # compute screening score (kkt)
  if (family == 'gaussian'){
    grad_without_target <- t(x) %*% (y - x %*% betas)
    grad_without_target <- abs(grad_without_target) / nrow(x)
  }
  
  if (family == 'binomial'){
    xb <- x %*% betas
    probs <- 1 / (1 + exp(-xb))
    grad_without_target <- t(x) %*% (y - probs)
    grad_without_target <- abs(grad_without_target) / nrow(x)
  }
  
  # create score data frame
  score <- as.vector(grad_without_target)
  feature_name <- rownames(grad_without_target)
  df_score <- data.frame(cbind(feature_name, score))
  df_score[,2] <- as.numeric(df_score[,2]) 
  # df_score <- df_score[order(df_score[,2], decreasing = TRUE),]
  
  # remove features in the active set
  df_score_active <- df_score
  df_score <- df_score[!(df_score[,1] %in% active_set),]
  
  # grab only w/n kkt threshold
  # kkt_eps <- sd(df_score$score)
  if (is.null(kkt_count)){
    idx_screen <- which(abs(df_score$score - lambda0) < kkt_eps)
    df_score_surrogate <- df_score[idx_screen,]
  } else {
    # grab top k candidate surrogates based on the sorted score
    idx_screen <- order(df_score$score, decreasing=TRUE)
    df_score_sorted <- df_score[idx_screen,] 
    df_score_surrogate <- df_score_sorted[1:kkt_count,]
  }
  
  # # plot screening
  # if(verbose){
  #   plot(sort(df_score_active$score),
  #        xlab="feature",
  #        ylab="gradient",
  #        main="Plot of gradient screening (target excluded)")
  #   abline(h=lambda0, col='red')
  #   abline(h=lambda0-kkt_eps, col='blue')
  # }
  
  # select features with score > threshold
  #  df_score_surrogate <- df_score[df_score$score > threshold,]
  kkt_set <- df_score_surrogate$feature_name
  
  out <- list(kkt_set=kkt_set, df_score=df_score)
  
  return(out)
}


# make fold_ids an input, if missing, then generate
surrogate_cross_validation_one_target <- function(target_feature, active_set, candidate_set, drop_set, x, y, lambdas, lambda0, alpha=1, family='gaussian', k=10, foldid=NULL, standardized=FALSE, verbose=FALSE){
  
  if(verbose){cat("\n", paste0("-----> Cross-Validation for ", target_feature, " candidate surrogates...", "\n"))}
  
  # get list of vectors of feature names for candidate surrogate models
  surrogate_set <- create_surrogate_set(active_set=active_set,
                                        candidate_set=candidate_set, 
                                        target_feature=target_feature,
                                        drop_set=drop_set)
  
  #--- Cross-validation for each surrogate model ---#
  
  # set up cross-validation folds if missing
  if (is.null(foldid)){
    foldid <- generate_folds(data=x, k=k)
  }
  
  cv_errors <- list()
  cv_sds <- c()
  surrogate_coefs <- list()
  
  for (i in 1:length(surrogate_set)){
    if(verbose){cat(paste0(i,'.'))}
    surrogate_name_i <- names(surrogate_set)[i]
    surrogate_set_i <- surrogate_set[[i]]
    cverr_i <- c()
    surrogate_coef <- c()
    
    for (j in 1:k){
      
      # train-test split
      idx_train <- (foldid != j)
      idx_test <- (foldid == j)
      x_ <- select_columns(x, surrogate_set_i)
      x_train <- x_[idx_train, ]
      x_test <- x_[idx_test, ]
      y_train <- y[idx_train]
      y_test <- y[idx_test]
      
      # train model on k-1 folds
      model <- glmnet(x=x_train, 
                      y=y_train, 
                      alpha=alpha, 
                      lambda=lambdas, 
                      family=family,
                      standardized=standardized)
      
      # test model on hold out fold
      yhat <- predict(model, newx=x_test, s=lambda0, type='response')
      
      # grab coefficient for the surrogate
      coefs <- coef(model, s=lambda0)
      surrog_idx <- which(rownames(coefs) == surrogate_name_i)
      surrogate_coef<- c(surrogate_coef, coefs[surrog_idx])
      
      # compute MSE
      if (family == "binomial"){
        # auroc
        roc_result <- suppressMessages(roc(response = as.vector(y_test), predictor = as.vector(yhat)))
        cverr_i[j] <- auc(roc_result)
      }
      
      if (family == "gaussian"){
        cverr_i[j] <- mean_squared_error(y_test, yhat)
      }
      
    }
    
    # average error across folds to get CV error
    cv_errors[[i]] <- c(surrogate_name_i, mean(cverr_i))
    cv_sds[i] <- sd(cverr_i)
    surrogate_coefs[[i]] <- list(surrogate_name_i, surrogate_coef)
  }
  #---------------------------------------------------#
  
  if(verbose){cat("\n")}
  
  # convert to data frame and sort by CV loss
  cv_errors_df <- as.data.frame(do.call(rbind, cv_errors))
  colnames(cv_errors_df) <- c("surrogate", "CV_loss")
  cv_errors_df[['CV_loss']] <- as.numeric(cv_errors_df[['CV_loss']])
  
  if (family == 'gaussian'){ 
    cv_errors_df <- cv_errors_df[order(cv_errors_df[['CV_loss']]),]
  }
  
  if (family == 'binomial'){ # for auroc: larger better 
    cv_errors_df <- cv_errors_df[order(cv_errors_df[['CV_loss']], decreasing=TRUE),]
  }
  
  #  Compute percent change from the base model
  # i.e. how much error goes up after replacing with surrogate variable
  null_loss <- cv_errors_df['CV_loss'][cv_errors_df['surrogate'] == '__base__']
  percent_change <- (cv_errors_df['CV_loss'] - null_loss) / null_loss
  cv_errors_df['percent_change_from_base'] <- percent_change
  
  # concat cv sd
  cv_errors_df['cv_sd'] <- cv_sds
  
  # output
  out <- list(cv_errors_df=cv_errors_df, surrogate_coefs=surrogate_coefs)
  
  return(out)
}


surrogate_search <- function(target_set, active_set, candidate_set, x, y, lambdas, lambda0, alpha=1, k=10, foldid=NULL, family, standardize=FALSE, kkt_eps=0.01, kkt_count=NULL, verbose=TRUE, hclust=TRUE, hclust_rho=0.7, kkt=TRUE, drop_one=FALSE){
  
  if(verbose){cat('Active set: ', active_set, '\n')}
  if(verbose){cat("\n")}
  
  # run hcluster and define cluster sets
  if (hclust){
    cluster_list <- clusterList(x=x, rho=hclust_rho, plotting=FALSE)
  }
  
  search_output <- list()
  search_results <- list()
  
  for (target_feature in target_set){
    
    if(verbose){cat('> Surrogate search for:', paste0(target_feature,'...'))}
    
    # 1) define drop set
    # default is the target feature b/c target is always included in cluster set k
    if (hclust){
      cluster_set_k <- cluster_list[[target_feature]]
      drop_set_idx <- which(active_set %in% cluster_set_k)
      drop_set <- active_set[drop_set_idx]
    } else {
      drop_set <- target_feature
    }
    
    # 2) run kkt screen
    if (kkt){
      screening_result <- kkt_screening(target_feature=target_feature, 
                                  x=x, 
                                  y=y, 
                                  lambdas=lambdas,
                                  lambda0=lambda0, 
                                  active_set=active_set, 
                                  drop_set=drop_set,
                                  kkt_eps=kkt_eps, 
                                  kkt_count=kkt_count,
                                  family=family, 
                                  standardize=standardize,
                                  verbose=verbose)
      
      candidate_set <- screening_result[['kkt_set']]
    }
    
    # 3) run CV
    result <- surrogate_cross_validation_one_target(x=x,
                                                    y=y,
                                                    lambda=lambda0,
                                                    lambdas=lambdas,
                                                    alpha=alpha,
                                                    family=family,
                                                    standardize=standardize,
                                                    active_set=active_set,
                                                    candidate_set=candidate_set,
                                                    drop_set=drop_set,
                                                    target_feature=target_feature,
                                                    foldid=foldid,
                                                    k=k,
                                                    verbose=verbose)
    
    search_results[[target_feature]] <- result[["cv_errors_df"]]
    if(verbose){cat("\n")}
  }
  
  ## output 
  search_output[["search_result"]] <- search_results
  search_output[["surrogate_coefficients"]] <- result[["surrogate_coefs"]]
  if (kkt){search_output[['screen_detail']] <- screening_result}
  
  return(search_output)
}






# active_set, target_feature, surrogate_feature should be index w.r.t x? 
# x has to be matrix with colnames

fit_model_base <- function(x, y, lambda, family, standardize=FALSE, active_set){
  
  idx_active <- which(colnames(x) %in% active_set)
  
  model <- glmnet(x=x[, idx_active],
                  y=y,
                  alpha=1,
                  lambda=lambda,
                  family=family,
                  standardize=standardize)
  
  out <- list(base_model=model, idx_active=idx_active)
  
  return(out)
}


fit_model_null <- function(x, y, lambda, family, standardize=FALSE, active_set, target_feature){
  
  # remove target feature
  active_set <- active_set[!(active_set %in% target_feature)]
  idx_active <- which(colnames(x) %in% active_set)
  
  model <- glmnet(x=x[, idx_active],
                  y=y,
                  alpha=1,
                  lambda=lambda,
                  family=family,
                  standardize=standardize)
  
  out <- list(null_model=model, idx_active=idx_active)
  
  return(out)
}



fit_model_surrogate <- function(x, y, lambda, family, standardize=FALSE, active_set, target_feature, surrogate_feature){

  # remove target feature
  active_set <- active_set[!(active_set %in% target_feature)]
  
  # include surrogate feature
  active_set <- c(active_set, surrogate_feature)   # <-- problem to append at the end?
  
  idx_active <- which(colnames(x) %in% active_set)
  idx_active <- sort(idx_active)
  
  model <- glmnet(x=x[, idx_active],
                  y=y,
                  alpha=1,
                  lambda=lambda,
                  family=family,
                  standardize=standardize)
  
  out <- list(surrogate_model=model, idx_active=idx_active)
  
  return(out)  
}


# for a given target and surrogate, fits base, null, surrogate model on train data
# then predicts on test data 
predict_model_all <- function(target_feature, surrogate_feature, active_set, x_train, y_train, x_test, lambda, family){
  
  model_base <- fit_model_base(x=x_train, 
                               y=y_train, 
                               lambda=lambda, 
                               family=family, 
                               active_set=active_set)
  
  model_null <- fit_model_null(x=x_train, 
                               y=y_train, 
                               lambda=lambda, 
                               family=family, 
                               active_set=active_set, 
                               target_feature=target_feature)
  
  model_surrogate <- fit_model_surrogate(x=x_train, 
                                         y=y_train, 
                                         lambda=lambda, 
                                         family=family, 
                                         active_set=active_set, 
                                         target_feature=target_feature,
                                         surrogate_feature=surrogate_feature)
  
  yhat_base <- predict(model_base$base_model, 
                       x_test[,model_base$idx_active], 
                       s=lambda, 
                       type = "response")
  
  yhat_null <- predict(model_null$null_model, 
                       x_test[,model_null$idx_active], 
                       s=lambda, 
                       type = "response")
  
  yhat_surrogate <- predict(model_surrogate$surrogate_model, 
                            x_test[,model_surrogate$idx_active], 
                            s=lambda, 
                            type = "response")
  
  out <- list(model_base=model_base,
              model_null=model_null,
              model_surrogate=model_surrogate,
              yhat_base=yhat_base,
              yhat_null=yhat_null,
              yhat_surrogate=yhat_surrogate)
  
  return(out)
}


# compute_mse_all
compute_mse_all <- function(y_test, predict_model_all_output){
  
  yhat_base <- predict_model_all_output$yhat_base
  yhat_null <- predict_model_all_output$yhat_null
  yhat_surrogate <- predict_model_all_output$yhat_surrogate
  
  loss_base <- mean((y_test - yhat_base)^2)
  loss_null <- mean((y_test - yhat_null)^2)
  loss_surrogate <- mean((y_test - yhat_surrogate)^2)
  
  condition_met <- (loss_base <= loss_null) & (loss_surrogate < loss_null)
  
  out <- list(good_surrogate=condition_met,
              loss_base=loss_base,
              loss_null=loss_null,
              loss_surrogate=loss_surrogate)
  
  return(out)
}


# compute_auroc_all
compute_auc_all <- function(y_test, predict_model_all_output){
  
  yhat_base <- predict_model_all_output$yhat_base
  yhat_null <- predict_model_all_output$yhat_null
  yhat_surrogate <- predict_model_all_output$yhat_surrogate
  
  roc_result_base <- roc(Y_test, yhat_base)
  loss_base <- auc(roc_result_base)
  roc_result_null <- roc(Y_test, yhat_null)
  loss_null <- auc(roc_result_null)
  roc_result_surrogate <- roc(Y_test, yhat_surrogate)
  loss_surrogate <- auc(roc_result_surrogate)
  
  # loss_base <- mean((Y_test - yhat_base)^2)
  # loss_null <- mean((Y_test - yhat_null)^2)
  # loss_surrogate <- mean((Y_test - yhat_surrogate)^2)
  
  condition_met <- (loss_base >= loss_null) & (loss_surrogate > loss_null)
  
  out <- list(good_surrogate=condition_met,
              loss_base=loss_base,
              loss_null=loss_null,
              loss_surrogate=loss_surrogate)
  
  return(out)
}


# compute_logloss_all
compute_auc_all <- function(y_test, predict_model_all_output){
  
  yhat_base <- predict_model_all_output$yhat_base
  yhat_null <- predict_model_all_output$yhat_null
  yhat_surrogate <- predict_model_all_output$yhat_surrogate
  
  roc_result_base <- roc(Y_test, yhat_base)
  loss_base <- auc(roc_result_base)
  roc_result_null <- roc(Y_test, yhat_null)
  loss_null <- auc(roc_result_null)
  roc_result_surrogate <- roc(Y_test, yhat_surrogate)
  loss_surrogate <- auc(roc_result_surrogate)
  
  # loss_base <- mean((Y_test - yhat_base)^2)
  # loss_null <- mean((Y_test - yhat_null)^2)
  # loss_surrogate <- mean((Y_test - yhat_surrogate)^2)
  
  condition_met <- (loss_base >= loss_null) & (loss_surrogate > loss_null)
  
  out <- list(good_surrogate=condition_met,
              loss_base=loss_base,
              loss_null=loss_null,
              loss_surrogate=loss_surrogate)
  
  return(out)
}


log_loss <- function(y, yhat) {
  epsilon <- 1e-15 # avoid log(0) errors
  yhat <- pmin(pmax(yhat, epsilon), 1 - epsilon)
  
  ll <- -mean(y * log(yhat) + (1 - y) * log(1 - yhat))
  return(ll)
}


# compute_logloss_all
compute_logloss_all <- function(y_test, predict_model_all_output){
  
  yhat_base <- predict_model_all_output$yhat_base
  yhat_null <- predict_model_all_output$yhat_null
  yhat_surrogate <- predict_model_all_output$yhat_surrogate
  
  loss_base <- log_loss(Y_test, yhat_base)
  loss_null <- log_loss(Y_test, yhat_null)
  loss_surrogate <- log_loss(Y_test, yhat_surrogate)
  
  condition_met <- (loss_base <= loss_null) & (loss_surrogate < loss_null)
  
  out <- list(good_surrogate=condition_met,
              loss_base=loss_base,
              loss_null=loss_null,
              loss_surrogate=loss_surrogate)
  
  return(out)
}


evaluate_on_test <- function(x_train, y_train, x_test, y_test, family, standardize=FALSE, lambda, active_set, surrogate_list){
  
  out_list <- list()
  list_idx <- 1
  true_count <- 0
  null_surrog_diff <- c()
  surrog_base_diff <- c()
  null_base_diff <- c()
  ratios <- c()
  
  surrogate_list_targets <- names(surrogate_list)
  
  for (i in 1:length(surrogate_list)){
    target_i <- surrogate_list_targets[i]
    surrogate_list_i <- surrogate_list[[i]]
    
    for (j in 1:length(surrogate_list_i)){
      
      surrogate_j <- surrogate_list_i[j] 
      
      preds_all <- predict_model_all(x_train=x_train, 
                                     y_train=y_train,
                                     x_test=x_test,
                                     lambda=lambda, 
                                     family=family, 
                                     active_set=active_set, 
                                     target_feature=target_i, 
                                     surrogate_feature=surrogate_j)
      
      if (family=="gaussian"){
        test_mse <- compute_mse_all(y_test=y_test, preds_all)
      }
      
      if (family=="binomial"){
        # test_mse <- compute_auc_all(y_test=y_test, preds_all)
        test_mse <- compute_logloss_all(y_test=y_test, preds_all)
      }

      out_list[[list_idx]] <- test_mse
      out_list[[list_idx]]["target"] <- target_i
      out_list[[list_idx]]$surrogate <- surrogate_j
      
      list_idx <- list_idx + 1     
      true_count <- true_count + test_mse$good_surrogate
    }
  }
  
  out_df <- do.call(rbind, lapply(out_list, as.data.frame))
  
  return(out_df)
  
}


