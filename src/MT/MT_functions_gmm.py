'''
(other file) Estimate Omega_0 and Omega_1 from ldsc results
Estimate mixing parameters pi_0, pi_1

'''


import itertools

import time
import random
from math import sqrt
import math
import numpy as np
from numpy.linalg import inv as np_inv # numpy matrix invert
from numpy import transpose as np_t # numpy matrix transpose


import scipy
from scipy.stats import norm
from scipy.stats import multivariate_normal
from scipy import linalg

from sklearn import mixture



# fit mixture model on all M snps to estimate pi_0, or set pi_0
# compute b_GWAS either by adding noise to true effects, or from geno/pheno
# betahats should only contain approximately independent snps

def update_gamma0( betahats, Omega_0, Omega_1, Cov_eL, pi_0, pi_1, K):

    '''E step, updates membership for vector of SNPs, returns vector for each component'''

    mvn_0L = [multivariate_normal(mean=np.zeros(K), cov=Omega_0 + Cov_e) for Cov_e in Cov_eL]
    phi0_betahats = np.asarray([mvn_0L[i].pdf(betahats[i]) for i in range(len(mvn_0L))])     # density at betahat under component 0


    mvn_1L = [multivariate_normal(mean=np.zeros(K), cov=Omega_1 + Cov_e) for Cov_e in Cov_eL]
    phi1_betahats = np.asarray([mvn_1L[i].pdf(betahats[i]) for i in range(len(mvn_1L))])     # density at betahat under component 1
    
    denoms = pi_0*phi0_betahats + pi_1*phi1_betahats  # each SNP has diff normalization factor
    denoms[denoms < 1e-20] = 1e-20 

# TODO currently set a small value to avoid divide by zero. need to figure out better way to avoid it/if it means model is bad, or if it's a snp that fits poorly for both models


    # TODO need to exclude poorly fitting snps

    gamma_0s = pi_0 * phi0_betahats/denoms
    gamma_1s = pi_1 * phi1_betahats/denoms

    return gamma_0s, gamma_1s


def update_pi( gamma_0s, gamma_1s, M, ignore_nans = True, verbose = False ):

    '''M step, returns updated pi_0 and pi_1'''

    if ignore_nans: # ignores any nans from a divide by 0
        # M2 = # of non-nan snps
        M2 = np.count_nonzero(~np.isnan(gamma_0s))
        if M2 == 0:
            print("TODO no non-nan snps, returning nan")
            return np.nan, np.nan

        pi_0 = 1.0*np.nansum(gamma_0s)/M2
        pi_1 = 1.0*np.nansum(gamma_1s)/M2
        if (M - M2) > 0 and verbose:
            print("ignoring %d SNPs where gamma_j was nan" % (M - M2))

    else:
        pi_0 = 1.0*np.sum(gamma_0s)/M
        pi_1 = 1.0*np.sum(gamma_1s)/M

    return pi_0, pi_1


# currently assume we have the true covariance matrices

def estimatePi(betahats, Omega_0, Omega_1, Cov_eL, K, M):

    maxIter = 1000 # TODO other maxIter values
    pi_0_est = 0.5
    pi_1_est = 1 - pi_0_est
    gamma_0s_est, gamma_1s_est = update_gamma0( betahats, Omega_0, Omega_1, Cov_eL, pi_0_est, pi_1_est, K  )

    for i in range(maxIter):

        if i % 100 == 0:
            print("ITERATION %d"%i)
            print(pi_0_est, pi_1_est, sum(gamma_0s_est), sum(gamma_1s_est))

        #gamma_0s_prev = gamma_0s_est # check for convergence using gamma_0s
        pi_0_prev = pi_0_est # used to check for convergence
        if np.isnan(pi_0_prev): # TODO figure out why nan
            print("pi_0_est is nan, breaking")
            break

        gamma_0s_est, gamma_1s_est = update_gamma0( betahats, Omega_0, Omega_1, Cov_eL, pi_0_est, pi_1_est, K  )
        pi_0_est, pi_1_est = update_pi( gamma_0s_est, gamma_1s_est, M, verbose=(i%100==0) )
        
        if np.abs(pi_0_prev - pi_0_est) < 1e-4:
         #and np.sum(np.abs(gamma_0s_prev - gamma_0s_est)) < 1e-8:
            print("converged 1e-4 at iter %d" % i)
            break

    return (pi_0_est, pi_1_est, gamma_0s_est, gamma_1s_est)





