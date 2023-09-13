'''
- (in other file) Estimate Omega_0 and Omega_1 from ldsc results
- (in other file) Estimate mixing parameters pi_0, pi_1

Estimate new effect sizes 

'''


import itertools

import time
import random
from math import sqrt
import math
import numpy as np
from numpy.linalg import inv as np_inv # numpy matrix invert
from numpy.linalg import pinv as np_pinv
from numpy.linalg import cond as np_cond
#from numpy import matmul as np_mul # numpy matrix mult, or use dot
from numpy import transpose as np_t # numpy matrix transpose
import sys


import scipy
from scipy.stats import norm
from scipy.stats import multivariate_normal
from scipy import linalg

from sklearn import mixture

# compute b_MTAG, b_MT from b_GWAS

def get_wts_t_mtag(t, Cov_g, Sigma_jL):
    '''Computes mtag weights for trait t from gwas summary statistics and genetic + environmental covariance'''
    wts_tL = []
    omega_t = Cov_g[:,t]
    omega_tt = Cov_g[t,t]
    for Sigma_j in Sigma_jL:
        G = Cov_g - np.outer(omega_t, omega_t)/omega_tt + Sigma_j
        if np_cond(G) < 1/sys.float_info.epsilon:
            A = np_inv(G)
        else:
            A = np_pinv(G)
        denom =  np.dot( np.dot(np_t(omega_t/omega_tt),A) , omega_t/omega_tt ) 
        wts_t = np.dot(np_t(omega_t/omega_tt), A) / denom
        wts_tL.append(wts_t)
    return wts_tL

# where the mtag prediction for trait t is np.dot(mtag_wts_t, betahat)

def get_gmm_wts(Omega_0, Omega_1, Sigma_jL):
    '''computes weights for posterior estimator of mean for each component, helper for estimating beta_MT'''

    # beta | betahat, gamma=c ~ N(A_c Vinv betahat, A_c)  - posterior dist of beta for gamma = c
    # let AV_c = (Omega_c^-1 + Sigma_j^-1)^-1 * Sigma_j^-1
    # computing AV_0 and AV_1 is same for all snp, so do it in a separate step
    
    Sigma_j_invL = [np_inv(Sigma_j) for Sigma_j in Sigma_jL]
    
    if np.linalg.cond(Omega_0) < 1/sys.float_info.epsilon:
        Omega_0_inv = np.linalg.inv(Omega_0)
    else:
        Omega_0_inv = np.linalg.pinv(Omega_0)
    
    if np.linalg.cond(Omega_1) < 1/sys.float_info.epsilon:
        Omega_1_inv = np.linalg.inv(Omega_1)
    else:
        Omega_1_inv = np.linalg.pinv(Omega_1)

    AV_0L = [np.dot(np_inv(Omega_0_inv + Sigma_j_inv) , Sigma_j_inv) for Sigma_j_inv in Sigma_j_invL] # wt on sample mean for gamma = 0
    AV_1L = [np.dot(np_inv(Omega_1_inv + Sigma_j_inv) , Sigma_j_inv) for Sigma_j_inv in Sigma_j_invL]
    
    return AV_0L, AV_1L


# computes E[beta_j | betahat_j, gamma_j = c]
# note: some values which are the same across SNPs should be pre-computed, then input
# for example, the numpy mvn objects and inverse var matrix
# TODO check whether E_given_betahat_0 is same as mtag wts time betahat
def get_beta_gmm(betahat, Omega_0, Omega_1, Sigma_j, pi_0, pi_1, mvn_0, mvn_1, AV_0, AV_1):
    '''GMM estimator for beta given betahat and mixture model parameters'''

    # find N(betahat_j | 0, Omega_c + V)
    phi0_betahat = mvn_0.pdf(betahat) # density at betahat under component 0
    phi1_betahat = mvn_1.pdf(betahat) # density at betahat under component 1
    
    E_given_betahat_0 = np.dot(AV_0 , betahat) # mean of beta | betahat, gamma = 0
    E_given_betahat_1 = np.dot(AV_1 , betahat)
    
    num = pi_0 * phi0_betahat * E_given_betahat_0 + pi_1 * phi1_betahat * E_given_betahat_1
    denom = pi_0 * phi0_betahat + pi_1*phi1_betahat
    
    posterior_beta = 1.0*num/denom
    return posterior_beta
    # TODO, if density at both models is very small, may get divide by 0


