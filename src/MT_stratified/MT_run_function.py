# Run MTAG and MT on UK Biobank data
# Output new effect sizes in UKBB 

# TODO split into multiple functions, e.g. if you only want MTAG estimate, or only want GMM, with or without writing to file
import statistics
import itertools
#import time
import random
from math import sqrt
import numpy as np
import sys

from MT_functions_cov import *
from MT_functions_gwas_bins import *
from MT_functions_gmm import *
from MT_functions_effects import *


epsilon = 1e-12 # the small non-zero value

np.set_printoptions(precision=3)




# function so you can run on diff sets of traits
# TODO Sigma_j should use N_miss to get SE for each snp, currently using same one for all snps
def run_MTAG_MT(exptName, traitL, gwasNamesL, t0, M_effective, outdir, gcov, h2, bins):

    outdir = outdir

    ########## Parse original GWAS file ############

    # get chr, pos info from original gwas in order rto write to file in PRS format
    with open(gwasNamesL[t0]) as f:
        gwasHeaderL = f.readline()
    isCaseCon = 'OR' in gwasHeaderL # whether the primary trait is case-control or not
        # tuple of column names to parse

    if isCaseCon:
        #plinkOutFields = ('SNP', 'A1', 'NMISS', 'OR', 'SE', 'STAT', 'P', 'A2', 'MAF')
        plinkOutFields = ('SNP', 'CHR', 'BP', 'A1', 'NMISS', 'OR', 'SE', 'STAT', 'P', 'A2') 
    else:
        #plinkOutFields = ('SNP', 'A1', 'NMISS', 'BETA', 'SE', 'STAT', 'P', 'A2', 'MAF')
        plinkOutFields = ('SNP', 'CHR', 'BP', 'A1', 'NMISS', 'BETA', 'SE', 'STAT', 'P', 'A2')


    snpInfoD, gwasHeaderL, fieldIndices = parseSumStats(gwasNamesL[t0], plinkOutFields)
    print('Parsed following fields from GWAS sumstats')
    print(plinkOutFields)
    
    # Read in GWAS data files
    # get just the snp names and betas (get rest of snp info later, could combine to only read through file once)
    # TODO


    snpIDsL = []
    betahatsL = []
    W_j_L_L = []


    for i in bins:
        snpIDs, betahats, W_j_L = readInGwasData(gwasNamesL, fieldIndices, gwasFormat='myplink', binName=i)
        snpIDsL.append(snpIDs)
        betahatsL.append(betahats)
        W_j_L_L.append(W_j_L)

    all_snpIDs = [snp for snpIDs in snpIDsL for snp in snpIDs]
    BETA_MTAG = []
    P_MTAG = []
    Z_MTAG = []

    # betahats is dim M x K
    for bin in range(len(snpIDsL)): 
        snpIDs = snpIDsL[bin]
        betahats = betahatsL[bin]
        W_j_L = W_j_L_L[bin]

        M = len(snpIDs)
        print("Bin 1...")
        print("M = %d snps" % M)

        ################################################


        # get genetic cov and estimation error cov from table of ldsc results
        gcovTableFile = gcov[bin]
        h2TableFile = h2[bin]
        pairInfoD, headerL_gcov, headerL_h2 = getPairInfo(gcovTableFile, h2TableFile)

        # Cov_g = # ldsc genetic covariance for trait i,j
        Cov_g = makeMatFromTable(traitL, pairInfoD, headerL_gcov, headerL_h2, "gcov")
        Cov_g = closestPSD(Cov_g)

        # Sigma_LD =  # ldsc intercept for trait i,j
        Sigma_LD = makeMatFromTable(traitL, pairInfoD, headerL_gcov, headerL_h2, "gcovIntercept")
        Sigma_LD = closestPSD(Sigma_LD)

        #print("Cov_g"); print(Cov_g)
        #print("Cov_g after closest PSD"); print(Cov_g)
        print("Sigma_LD"); print(Sigma_LD)
        print("Sigma_LD after closestPSD"); print(Sigma_LD)

        Sigma_jL = []
        W_invL = []
        K = len(traitL)

        for snpindex in range(len(snpIDs)) : 

            # W is a matrix with sqrt(N) along diag, used to convert Sigma_LD to Sigma_j
            W_j = W_j_L[snpindex]
            W_inv = np.linalg.inv(W_j)
            # overall Sigma_j, for snps with no missing genos
            Sigma_j = np.dot(np.dot(W_inv, Sigma_LD), W_inv.T)

            W_invL.append(W_inv)
            Sigma_jL.append(Sigma_j)
            #print("Sigma_j"); print(Sigma_j)


        # Msr is a matrix with sqrt(M_effective) along diag, used to convert Cov(XB_t,XB_t') to Cov(B_jt, B_jt')
        Msr = np.zeros((K,K)); np.fill_diagonal(Msr, 1.0/np.sqrt(M_effective))
    
    
        Cov_g = np.zeros((K,K))
        for i in range(len(snpIDs)):
            Cov_g = Cov_g + np.matmul(betahats[i,:],np.transpose(betahats[i,:])) - Sigma_jL[i]
   

        #Omega_1 = np.matrix([[4.17908012e-06, 3.85316097e-06, 3.44092853555e-06],[3.853160966e-06, 4.0626684619e-06, 3.3126080013e-06],[3.4409285355e-06, 3.312608001e-06, 3.72338694557e-06]])
        Omega_1 = 1.0/M * Cov_g
        #Omega_1 = np.dot(np.dot(Msr, Cov_g), Msr.T)
        print('Omega_1'); print(Omega_1)
        Omega_1 = closestPSD(Omega_1)
        print('Omega_1 after closestPSD'); print(Omega_1)

        #fname = outdir + '/' + exptName + '.Sigma_j'
        #numpy.savetxt(fname, Sigma_j, delimiter = '\t', header='\t'.join(traitL) )
        fname = outdir + '/' + exptName + '.Omega_1' + '_' + str(bin)
        numpy.savetxt(fname, Omega_1, delimiter = '\t', header='\t'.join(traitL) )




        ################# Estimate effects (MTAG) ###################
        # MTAG # TODO also depends on M_eff
        wts_mtagL = get_wts_t_mtag(t0, Omega_1, Sigma_jL)

        ########## The following parts needs to be changed for MT ##########
        newStats = np.empty(len(snpIDs))
        newZs = np.empty(len(snpIDs))
        newPvals = np.empty(len(snpIDs))
        beta_mtags_k = np.empty(len(snpIDs))

        print(betahats.shape)
        print(wts_mtagL[0].shape)
        print(len(wts_mtagL))
        print(np.mean(np.asarray(wts_mtagL), axis = 0))

        for i in range(len(wts_mtagL)):
            beta_mtags_k[i] = np.matmul(betahats[i,:], np_t(wts_mtagL[i]))
            newStats[i] = beta_mtags_k[i]
            snp = snpIDs[i]
            gwasLineL = snpInfoD[snp] # includes gwas beta and pval
            SE = gwasLineL[plinkOutFields.index('SE')] # TODO get index of SE from header]
            SE = float(SE)
            newZs[i] = newStats[i]/SE
            newPvals[i] = 2*norm.cdf(-1.0*abs(newZs[i]))
            BETA_MTAG.append(beta_mtags_k[i])
            P_MTAG.append(newPvals[i])
            Z_MTAG.append(newZs[i])


        #mtagOut = outdir + "/" + exptName + ".MTAGest" + "_" + str(bin)
        #print('Estimated mtag effects')

        

    gmmOut = outdir + "/" + exptName + ".assoc.MTAG.MT_estimates"
    with open(gmmOut, 'w') as outf:
        outf.write('\t'.join(plinkOutFields) + '\t') # write header
        outf.write('BETA_MTAG' + '\t')
        outf.write('Z_MTAG' + '\t')
        outf.write('P_MTAG' + '\n') 
        for i in range(len(all_snpIDs)):
            # replace BETA/OR, STATS, and P field with new value (or append to end or line if not already in the field)
            snpName = all_snpIDs[i]
            outf.write('\t'.join(snpInfoD[snpName]) + '\t')
            outf.write(str(BETA_MTAG[i]) + '\t')
            outf.write(str(Z_MTAG[i]) + '\t')
            outf.write(str(P_MTAG[i]) + '\n')


