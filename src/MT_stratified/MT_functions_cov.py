import itertools

import time
import random
from math import sqrt
import math
import numpy as np
from numpy.linalg import inv as np_inv # numpy matrix invert
#from numpy import matmul as np_mul # numpy matrix mult, or use dot
from numpy import transpose as np_t # numpy matrix transpose


import scipy
from scipy.stats import norm
from scipy.stats import multivariate_normal
from scipy import linalg

import sklearn
from sklearn import mixture


def getPairInfo(gcovTableFile, h2TableFile):
    '''given table of ldsc results for each pair, get line corresponding to each pair in a dict, and header for table column names'''

    pairInfoD = {} # pairInfoD[(t1,t2)] = pairInfoD[(t2,t1)] = lineL from table
    with open(gcovTableFile) as f:
        headerL_gcov = f.readline().rstrip().split()

        for line in f:
            lineL = line.rstrip().split()
            a = lineL[0]; b = lineL[1] # add each pair of traits by name
            pairInfoD[(a,b)] = tuple(line.rstrip().split())
            pairInfoD[(b,a)] = tuple(line.rstrip().split())

    # also get a,a pair (e.g. armfat, armfat)?
    with open(h2TableFile) as f:
        headerL_h2 = f.readline().rstrip().split()
        for line in f:
            lineL = line.rstrip().split()
            a = lineL[0] # add trait name
            pairInfoD[(a,a)] = tuple(line.rstrip().split())

    return pairInfoD, headerL_gcov, headerL_h2


def makeMatFromTable(traitsL, pairInfoD, headerL_gcov, headerL_h2, colName):
    '''return genetic correlation and estimatation error correlation (or whatever) from a tsv of trait_a, trait_b, gcorr, se '''

    K = len(traitsL)
    Cov_g = np.zeros((K,K))
    Sigma_LD = np.zeros((K,K))
    missingCodes = ["NA", "nan", "-9"]


    Mat = np.zeros((K,K))
    colIndex_pair = headerL_gcov.index(colName) # e.g. "gcov"
    #"gcov", "gcovSE", "gcovIntercept", "gcovInterceptSE"  
    # not using SE at moment but good to report

    # for cov(a,a), use h2 instead of gcov as header, and intercept instead of gcov intercept

    if colName == 'gcov':
        colName_single = 'h2'
    if colName == 'gcovIntercept':
        colName_single = 'intercept'
    if colName == 'gcovSE':
        colName_single = 'h2SE'
    if colName == 'gcovInterceptSE':
        colName_single = 'interceptSE'

    colIndex_single = headerL_h2.index(colName_single)

    for i in range(K):
        for j in range(K):
            a = traitsL[i]; b = traitsL[j]
            if a == b: # if pair is same trait
                val_ij = pairInfoD[(a,a)][colIndex_single]
            else:
                val_ij = pairInfoD[(a,b)][colIndex_pair]

            if val_ij in missingCodes: 
                Mat[i,j] = np.nan
            else:
                Mat[i,j] = float(val_ij)
    return Mat



def closestPSD(A):
    '''truncated SVD to get closest PSD matrix to A. Assumes A is symmetric'''
    s, V = np.linalg.eigh(A)
    #print(s) # if you want to print eigenvalues before and after
    s = np.maximum(s, 0.0)
    #print(s)
    B = np.matmul(np.matmul(V, np.diag(s)), np.transpose(V))
    return B






''' 
# old code
h2File = '/u/flashscratch/l/lgai/MultiTrait_UKBB/data/Neale_trait_h2.txt'
Corr_g = get_Corr(gName, traitsL)
Corr_e = get_Corr(pheName, traitsL)
tausArr, vsArr = get_ge_sds(h2File, traitsL, M_effective, N_gwas)

epsilon = 1e-12 # check that epsilon is smaller than entries in Cov_g

Cov_g = Corr_g * np.outer(tausArr, tausArr)
Cov_e = Corr_e * np.outer(vsArr, vsArr)


tausqArr = tausArr**2
vsqArr = vsArr**2

# get correlations from pairwise file, get h2 from text file
# e.g. cols of "trait_i trait_j corr"
# if additional columns after that, they are ignored
def get_Corr(corrFile, traitsL):
    K = len(traitsL)
    CorrMat = np.ones((K,K)) # will fill off-diag with corr
    trait2indexD = {} # traitsD[t] = index in traitsL
    for i in range(K):
        trait2indexD[traitsL[i]] = i 

    with open(corrFile) as f:
        for line in f:
            ti, tj, corr = line.rstrip().split()[:3]
            
            if ti in trait2indexD and tj in trait2indexD:
                i = trait2indexD[ti]
                j = trait2indexD[tj]
                corr = max(min(float(corr),1), -1) 
                CorrMat[i,j] = corr
                CorrMat[j,i] = corr

    return CorrMat

# get sqrt of var(beta_j) and var(beta_j|beta_j) in all traits
def get_ge_sds(h2File, traitsL, M_eff, N_gwas):
    K = len(traitsL)
    heritsD = {}
    with open(h2File) as f:
        for line in f:
            lineL = line.rstrip().split() # either tname, h2, se or just tname, h2
            t = lineL[0]; h2 = lineL[1]
            heritsD[t] = min(float(h2),1)

    tausqArr = np.zeros(K) # var of beta_j
    vsqArr = np.zeros(K) # var of betahat_j given beta_j

    for i in range(K):
        tausqArr[i] = 1.0*heritsD[traitsL[i]]/M_eff 
        vsqArr[i] = max(0, 1 - heritsD[traitsL[i]])/N_gwas[i] 
    
    return np.sqrt(tausqArr), np.sqrt(vsqArr)

'''



