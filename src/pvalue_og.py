import sys
import numpy as np
import pandas as pd
from glmnet import GaussNet 
from sklearn.preprocessing import StandardScaler
# import matplotlib.pyplot as plt
from scipy.stats import norm as normal_dbn
from glmnet.inference import constraints
# need to install mpmath
# rng = np.random.default_rng(2)

def active_constraints(X, 
                       active, 
                       signs, 
                       lambda_value,
                       penalty_factor=None,
                       weights=None):
    n, p = X.shape
    if weights is None:
        weights = np.ones(n)
    weights /= weights.sum()

    X = X * np.sqrt(weights[:,None])
    
    if penalty_factor is None:
        penalty_factor = np.ones(X.shape[1])
    X_A = X[:, active]
    Q_A = np.linalg.inv(X_A.T @ X_A)
    L_A = -np.diag(signs) @ Q_A @ X_A.T
    O_A = signs * (Q_A @ (penalty_factor[active] * 
                   signs * lambda_value))
    L_A = L_A * np.sqrt(weights[None, :])
    return L_A, O_A, weights


def inactive_constraints(X, 
                         active, 
                         signs, 
                         lambda_value,
                         penalty_factor=None,
                         weights=None):

    n, p = X.shape

    if weights is None:
        weights = np.ones(n)
    weights /= weights.sum()
    X = X * np.sqrt(weights[:,None])

    if penalty_factor is None:
        penalty_factor = np.ones(p)

    inactive = np.ones(p, bool)
    inactive[active] = 0
    X_I = X[:,inactive]

    X_A = X[:, active]
    Q_A = np.linalg.inv(X_A.T @ X_A)

    irrep = X_I.T @ X_A @ Q_A @ (penalty_factor[active] * 
                                 signs * lambda_value)
    X_ImA = X_I - X_A @ Q_A @ X_A.T @ X_I

    RHS = penalty_factor[inactive] * lambda_value
    O_I = np.hstack([RHS - irrep, 
                     RHS + irrep])
    L_I = np.concatenate([X_ImA.T, -X_ImA.T])
    L_I = L_I * np.sqrt(weights[None, :])
    return L_I, O_I


def add_drop_directions(X, active, var_drop):

    v = var_drop
    
    n, p = X.shape
    active_v = np.zeros(p, bool)
    active_v[active] = 1
    active_v[v] = 0
    X_Av = X[:,active_v]

    inactive = np.ones(p, bool)
    inactive[active] = 0
    
    Q_Av = np.linalg.inv(X_Av.T @ X_Av)
    X_add = X[:,inactive]
    X_add = X_add - X_Av @ Q_Av @ X_Av.T @ X_add
    X_add /= np.sqrt((X_add**2).sum(0))[None,:]
    
    return X_add


def truncnorm_pval(L, 
                   Z,
                   U,
                   S, 
                   meanZ=0,
                   alternative='twosided'):

    N = (normal_dbn.cdf((Z-meanZ)/S) - normal_dbn.cdf((L-meanZ)/S))
    D = (normal_dbn.cdf((U-meanZ)/S) - normal_dbn.cdf((L-meanZ)/S))
    P = N / D

    if alternative == 'greater':                                            
        return 1 - P                                                        
    elif alternative == 'less':                                             
        return P                                                            
    else:                                                                   
        return max(2 * min(P, 1-P), 0)                                      

    return N/D


def pvalues(X, 
            y,
            active,
            con):

    n, p = X.shape
    
    inactive = np.ones(p, bool)
    inactive[active] = 0
    pvalues = pd.DataFrame(np.zeros((p, len(active)))* np.nan,
                           columns=active,
                           index=np.arange(p))
    for i, v in enumerate(active):

        X_add = add_drop_directions(X, active, v)
        bounds = [con.bounds(X_add[:,i], y) 
                  for i in range(X_add.shape[1])]
        
        pvalues.iloc[inactive,i] = np.array([truncnorm_pval(*bound_info) for
                                             bound_info in bounds])
    return pvalues


def main(X, 
         y,
         active, 
         signs, 
         lambda_value,
         sigma,
         penalty_factor=None,
         weights=None):
    
    n, _ = X.shape
    
    # compute active constraints
    L_A, O_A, W = active_constraints(X=X,
                                  active=active,
                                  signs=signs,
                                  lambda_value=lambda_value,
                                  penalty_factor=penalty_factor,
                                  weights=weights)
 #   assert((L_A @ y - O_A).max() < 0)   
    
    # compute inactive constraints
    L_I, O_I = inactive_constraints(X=X,
                                active=active,
                                signs=signs,
                                lambda_value=lambda_value,
                                penalty_factor=penalty_factor,
                                weights=W)
 #   assert((L_I @ y - O_I).max() < 0)

    # combine constraints
    L = np.concatenate([L_A, L_I])
    O = np.hstack([O_A, O_I])
 #   assert((L @ y - O).max() < 0)
    con = constraints(L, O, covariance=sigma**2 * np.identity(n))
    
    pval = pvalues(X=X, y=y, active=active, con=con)

    return(pval)


if __name__ == "__main__":

    # Retrieve the directory path from the command line arguments
    data_dir = sys.argv[1]

    # Load data from CSV files
    X = pd.read_csv(f"{data_dir}/X.csv")
    X = X.values
    y = pd.read_csv(f"{data_dir}/y.csv")
    y = y['y'].values
    active_set = pd.read_csv(f"{data_dir}/active_set.csv")['active_set'].values
    signs = pd.read_csv(f"{data_dir}/signs.csv")['signs'].values

    with open(f"{data_dir}/lambda0.txt", 'r') as file:
        lambda0 = float(file.read().strip())

    with open(f"{data_dir}/sigma.txt", 'r') as file:
        sigma = float(file.read().strip())

    # run selective inference    
    pvals = main(X=X, 
                 y=y,
                 active=active_set, 
                 signs=signs, 
                 lambda_value=lambda0,
                 sigma=sigma,
                 penalty_factor=None,
                 weights=None)

    # Define the path to save the p-values CSV
    output_file_path = f"{data_dir}/pvals.csv"

    # Save the DataFrame to a CSV file
    pvals.to_csv(output_file_path, index=False)

    print(f"Saved p-values to {output_file_path}")
