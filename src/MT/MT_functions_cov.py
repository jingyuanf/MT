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

    '''return genetic correlation and estimatation error correlation from a tsv of trait_a, trait_b, gcorr, gcorrse, h2, h2se '''

    K = len(traitsL)
    Cov_g = np.zeros((K,K))
    Sigma_LD = np.zeros((K,K))
    missingCodes = ["NA", "nan", "-9"]


    Mat = np.zeros((K,K))
    colIndex_pair = headerL_gcov.index(colName) 

    # e.g. "gcov"
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





