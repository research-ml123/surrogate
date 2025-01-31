# Download data from: https://archive.ics.uci.edu/dataset/17/breast+cancer+wisconsin+diagnostic

############################
#   Simulation parameters  #
############################
SEED <- 123
N_TRAIN <- 500
KKT_COUNT <- 30
FAMILY <- 'binomial'
LAMBDA_TYPE <- 'lambda.min'    
SAVE <- FALSE
FILENAME <- "breast_cancer_result.RData"

#################################
#   Load scripts and libraries  #
#################################
library(here)
library(glmnet)
library(pROC)
source(here("src", 'surrogate_lasso.R'))
source(here("src", 'old/helper.R'))

################
#   Load data  #
################
set.seed(SEED)

df <- read.csv(here::here('data/real/breast+cancer+wisconsin+diagnostic/wdbc.data'), 
               header=FALSE)
df <- df[,-1]
X <- as.matrix(df[,-1])
Y <- ifelse(df$V2=="M", 1, 0)

# train-test split
idx_train <- sample(1:nrow(X), size=N_TRAIN, replace = FALSE)
X_train <- X[idx_train, ]
X_test <- X[-idx_train, ]
Y_train <- as.matrix(Y[idx_train])
Y_test <- as.matrix(Y[-idx_train])

# standardize X
X_train <- scale(X_train, center=TRUE, scale=TRUE)
X_test <- scale(X_test, center=TRUE, scale=TRUE)


########################
#   Fit initial model  #
########################
# generate fold id for cv
k=10
foldid <- generate_folds(data=X_train, k=k)

# find lambda0 through cross-validation
cvfit0 <- cv.glmnet(x=X_train, 
                    y=Y_train, 
                    family=FAMILY, 
                    alpha=1, 
                    foldid=foldid,
                    standardize=FALSE)
# plot(cvfit0)

# fix lambda0 for all the following steps
if (LAMBDA_TYPE == "lambda.1se"){
  lambda0 <- cvfit0$lambda.1se
} 

if (LAMBDA_TYPE == "lambda.min"){
  lambda0 <- cvfit0$lambda.min 
} 

lambdas <- cvfit0$lambda

# get nonzero beta indices
betas <- coef(cvfit0, s=lambda0)[-1]
idx_nonzero <- which(abs(betas) > 0)

# define feature sets
active_set <- colnames(X)[idx_nonzero]
candidate_set <- setdiff(colnames(X), active_set)

###########################
#   Run surrogate search  #
###########################
surrogate_result <- surrogate_search(target_set=active_set,   # active_set if you want all
                                     active_set=active_set,
                                     candidate_set=candidate_set,
                                     x=X_train,
                                     y=Y_train,
                                     lambdas=lambdas,
                                     lambda0=lambda0,
                                     alpha=1,
                                     k=10,
                                     foldid=NULL,
                                     family=FAMILY,
                                     verbose=FALSE,
                                     kkt=TRUE,
                                     kkt_eps=0.002,
                                     kkt_count=KKT_COUNT,
                                     hclust=FALSE, 
                                     hclust_rho=0.7)

surrogate_list <- get_surrogates(surrogate_result)

#############################################
#   evaluate surrogate models on test set   #
#############################################
test_eval_df <- evaluate_on_test(x_train=X_train, 
                                 y_train=Y_train,
                                 x_test=X_test, 
                                 y_test=Y_test, 
                                 family=FAMILY, 
                                 standardize=FALSE, 
                                 lambda=lambda0, 
                                 active_set=active_set, 
                                 surrogate_list=surrogate_list)

mean((test_eval_df$loss_null - test_eval_df$loss_surrogate ) / test_eval_df$loss_null)

test_eval_df$percent_change <- (test_eval_df$loss_null - test_eval_df$loss_surrogate) / test_eval_df$loss_null

#########################
#      Results          #
#########################
df_true <- test_eval_df[test_eval_df$good_surrogate,]
print(df_true)

# proportion of true surrogates among all surrgoates identified
prop <- sum(test_eval_df$good_surrogate) / nrow(test_eval_df)
print(paste0("Proportion of true surrogates: ", prop))


# AMONG THE SURROGATES THAT LOOK GOOD
loss_perc_change <- (df_true$loss_null - df_true$loss_surrogate ) / df_true$loss_null

count_better <- sum(loss_perc_change  <= 0.01)
prop_better_1pct <- sum(loss_perc_change  >= 0.01) / nrow(test_eval_df)
prop_better_10pct <- sum(loss_perc_change  >= 0.1) / nrow(dtest_eval_df)
mean_better <- mean(loss_perc_change)
print(paste0("Proportion of surrogates better than null (>1percent): ", prop_better_1pct))
print(paste0("Proportion of surrogates better than null (>10percent): ", prop_better_10pct))
print(paste0("Average delta: ", mean_better ))
print(paste0("Number of surrogates"))


#########################
#      OUTPUT           #
#########################
out <- list()
out$result_df <- test_eval_df
out$result_df_true <- df_true
out$count_over1pct <- count_better
out$prop_better_1pct <- prop_better_1pct
out$prop_better_10pct <- prop_better_10pct
out$avg_delta <- mean_better

# if (SAVE){
#   output_dir <- here::here("data/output/")
#   output_file_path <- paste0(output_dir, FILENAME) 
#   save(out, file=output_file_path)
# }
