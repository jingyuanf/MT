#!/usr/bin/env python
# coding: utf-8

# In[15]:


import sys
import getopt
import argparse
import itertools
import collections
import numpy as np
import numpy.random
import scipy.stats
from scipy.stats import multivariate_normal as mvnorm # uses python2.7
import random 


# In[5]:


def get_N_gwas(traitsL):
    '''return dict of traitCtD[trait] = # non missing for that gwas
    '''
    traitCtD = {}
    # use the count file corresponding to the cohort used for gwas, which may be a subset of all available ukbb
    countFile = "/u/flashscratch/l/lgai/MultiTrait_UKBB/data_pheno/caseControlCounts.rand150k_indivs.txt"
    with open(countFile) as f:
        header = f.readline()
        # pheno case    control caseConTotal    nonmissing
        for line in f:
            lineL = line.rstrip().split()
            traitCtD[lineL[0]] = float(lineL[-1])

    N_gwas = np.zeros(len(traitsL))
    for i in range(len(traitsL)):
        N_gwas[i] = traitCtD[traitsL[i]]

    return N_gwas


# In[3]:

'''
snpIDsL, betahatsL, W_j_L_L

for i in bins:
	snpIDs, betahats, W_j_L = readInGwasData(..., ..., ..., i)
	snpIDsL.append(snpIDs)
	..
	..
'''

def readBin(binName):
    snpIDs = set()

    with open(binName) as f:
        f.readline() # skip header
        for line in f:
            lineL = line.rstrip().split()
            snpIDs.add(lineL[1]) #hardcoded

    return snpIDs

def readInGwasData(fileNameL, fieldIndices, binName, gwasFormat="myplink"):
    nTraits = len(fileNameL)

    # get consensus set of snps
    snpSetsL = []
    for t in range(nTraits):
        snpSet_t = getSnpIdsTrait(fileNameL[t], gwasFormat, fieldIndices)
        snpSetsL.append(snpSet_t)
    binSnpSet = readBin(binName)
    snpSetsL.append(binSnpSet)
    snpIDs = set.intersection(*snpSetsL)
    nSnp  = len(snpIDs)
    
    print(nSnp)
    betahats = np.empty((nSnp, nTraits))
    N = np.empty((nSnp, nTraits))

    # get summary stats
    for t in range(nTraits):
        (snpIDsL, betahats_t, N_t) = getSummaryStatsTrait(fileNameL[t], gwasFormat, snpIDs, fieldIndices)
        betahats[:,t] = betahats_t
        N[:,t] = N_t
        
    W_j_L = []
    
    
    for snp in range(len(snpIDs)):
        W_j = np.zeros((nTraits,nTraits))
        for t in range(nTraits):
            W_j[t,t] = np.sqrt(N[snp,t])
        W_j_L.append(W_j)
    
    '''
    for snp in range(len(snpIDs)):
        W_j = np.zeros((nTraits,nTraits))
        W_j[0,0] = N[snp,0]
        W_j[1,1] = 147000
        W_j[2,2] = 149000
        W_j_L.append(W_j)    
    '''
    # return consensus snps and betahats
    return snpIDsL, betahats, W_j_L



# In[8]:


def getSnpIdsTrait(fileName, gwasFormat, fieldIndices):
    snpIDs = set()
    #SE_betahats_t = []

    with open( fileName ) as f:
        f.readline() # skip header
        for line in f:
            lineL = line.rstrip().split()
            if gwasFormat=="myplink":
                snpIDs.add(lineL[fieldIndices['SNP']])
            else:
                print("ERROR: '%s' is not a recognized option for GWAS format" % gwasFormat)
                exit(1)
    return snpIDs


# In[9]:


# get betahat (and pval?) 
def getSummaryStatsTrait(fileName, gwasFormat, snpIDs_consensus, fieldIndices):
    print("Reading trait %s..." % fileName.split('/')[-1])
    snpIDsL = []
    betahats_t = []
    N_t = []
    #SE_betahats_t = []

    with open( fileName ) as f:
        f.readline() # skip header
        for line in f:
            lineL = line.rstrip().split()
            if gwasFormat=="myplink":
                if lineL[fieldIndices['SNP']] in snpIDs_consensus and lineL[fieldIndices['SNP']] not in snpIDsL:
                    snpIDsL.append(lineL[fieldIndices['SNP']]) # TODO assumes snps in same order in all files
                    betahats_t.append(float(lineL[fieldIndices['BETA']])) # may also be odds ratio.
                    N_t.append(float(lineL[fieldIndices['NMISS']]))
                    
                #TODO read header and munge the BETA/OR field, which may be in different position for different plink options
                
            else:
                print("ERROR: '%s' is not a recognized option for GWAS format" % gwasFormat)
                exit(1)

    print(len(betahats_t))
    print(len(snpIDsL))
    return (snpIDsL, betahats_t, N_t)


# In[24]:


# e.g. fields = ('SNP', 'CHR', 'BP', 'A1', 'NMISS', 'OR', 'SE', 'STAT', 'P') 
# fields must contain 'SNP' as the first entry, I've hardcoded the fields names and order in MT_run_ukbb
def parseSumStats(fileName, fields):
    '''returns of dict of snps with relevant info. assumes plink sumstat file, assumes already filtered for missingness. could adjust the column names otherwise'''

    snpD = {} # snpD[snpName] = (snp, chr, pos, a1, ..., p)
    # tuple of strings

    with open(fileName) as f:
        
        headerL = f.readline().rstrip().split()
        fieldIndicesL = [headerL.index(col) for col in fields]
        fieldIndices = {}
        for col in headerL:
            fieldIndices[col] = headerL.index(col)
        snpIndex = headerL.index('SNP')
        
        print(fieldIndices)

        for line in f:
            lineL = line.rstrip().split()
            snpName = lineL[snpIndex]
            snpD[snpName] = tuple([lineL[i] for i in fieldIndicesL])

    return snpD, headerL, fieldIndices  # also return all cols contained in original gwas file


# In[30]:

