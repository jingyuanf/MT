# Run MTAG and MT on UK Biobank data
# Output new effect sizes in UKBB 

# TODO split into multiple functions, e.g. if you only want MTAG estimate, or only want GMM, with or without writing to file

import statistics
import itertools
import random
from math import sqrt
import numpy as np
import sys

from MT_functions_cov import *
from MT_functions_gwas import *
from MT_functions_gmm import *
from MT_functions_effects import *


epsilon = 1e-12 # the small non-zero value

np.set_printoptions(precision=3)


# function so you can run on diff sets of traits
# TODO Sigma_j should use N_miss to get SE for each snp, currently using same one for all snps

def run_MTAG_MT(exptName, traitL, gwasNamesL, t0, M_effective, outdir, gcov, h2):

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
    snpIDs, betahats, W_j_L = readInGwasData(gwasNamesL, fieldIndices, gwasFormat='myplink')


    # betahats is dim M x K

    M = len(snpIDs)
    print("M = %d snps" % M)
    ################################################


    # get genetic cov and estimation error cov from table of ldsc results
    gcovTableFile = gcov
    h2TableFile = h2
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
    fname = outdir + '/' + exptName + '.Omega_1'
    numpy.savetxt(fname, Omega_1, delimiter = '\t', header='\t'.join(traitL) )



    # TODO take in which Omega_0 as an arg, or just try all of them

    #### multiple options for Omega_0 ####
    # a. sparse in other traits 
    # b. sparse in main trait
    # c. non-sparse but main trait has no corr with other traits
    # d. non-sprase but no correlation between any traits

    Omega_0a = np.zeros_like(Omega_1)
    Omega_0a = np.zeros((K,K))
    np.fill_diagonal(Omega_0a, epsilon)

    Omega_0a[t0,t0] = Omega_1[t0,t0]

    Omega_0b = np.copy(Omega_1) # shallow copy
    Omega_0b[t0,:] = 0; Omega_0b[:,t0] = 0
    Omega_0b[t0,t0] = epsilon

    Omega_0c = np.copy(Omega_1)
    Omega_0c[t0,:] = 0; Omega_0c[:,t0] = 0
    Omega_0c[t0,t0] = Omega_1[t0,t0]

    Omega_0d = np.zeros_like(Omega_1)
    np.fill_diagonal(Omega_0d, np.diag(Omega_1))


    Omega_0s = [Omega_0a, Omega_0b, Omega_0c, Omega_0d]

    isVerbose = False
    if isVerbose:
        print("Omega_1"); print(Omega_1)
        #print("Sigma_j"); print(Sigma_j)

        for i in range(4):
            print('Omega0' + 'abcd'[i])
            print(Omega_0s[i])



    
    # Note: str(beta_est) preserves more digits than %f

    # np arrays to store MT and GMM statistics
    newStats = np.zeros((M, 1+len(Omega_0s)))
    newPvals = np.zeros((M, 1+len(Omega_0s)))  # univar normal pvals for primary trait

    ################# Estimate effects (MTAG) ###################
    # MTAG # TODO also depends on M_eff
    wts_mtagL = get_wts_t_mtag(t0, Omega_1, Sigma_jL)
    beta_mtags_k = np.empty(len(snpIDs))
    print(betahats.shape)
    print(wts_mtagL[0].shape)
    print(len(wts_mtagL))
    print(np.mean(np.asarray(wts_mtagL), axis = 0))
    for i in range(len(wts_mtagL)):
        beta_mtags_k[i] = np.matmul(betahats[i,:], np_t(wts_mtagL[i]))

    mtagOut = outdir + "/" + exptName + ".MTAGest"
    newStats[:,0] = beta_mtags_k
    print('Estimated mtag effects')



    ########## Do MT for each Omega_0 ############

    # write pi_0 and pi_1 estimates to same file for all Omega_0s
    piOut = outdir + "/" + exptName + ".Omega_0s.pi"
    with open(piOut, 'w') as outf:
        outf.write("Omega0name\tpi0est\tpi1est\n")


    for Omega_0_index in range(4):
        Omega_0_name = "Omega_0" + "abcd"[Omega_0_index]
        Omega_0 = Omega_0s[Omega_0_index]

        ############ Estimate mixing weights (GMM) ############

        pi_0_est, pi_1_est, gamma_0s_est, gamma_1s_est = estimatePi(betahats, Omega_0, Omega_1, Sigma_jL, K, M)

        print('Estimated mixing weights for %s' % Omega_0_name)
        print(pi_0_est, pi_1_est)


        with open(piOut,'a') as outf:
            outf.write('%s\t%s\t%s\n' % (Omega_0_name, str(pi_0_est), str(pi_1_est)) )

        # p(gamma | betahat) for each snp and component
        gammaOut = outdir + "/" + exptName + "." + Omega_0_name + ".gammas"

        with open(gammaOut,'w') as outf:
            outf.write('snp_id\tgamma_0\tgamma_1\n')
            for i in range(M):
                outf.write('\t'.join([snpIDs[i], str(gamma_0s_est[i]), str(gamma_1s_est[i])]) + '\n')


        ############# Estimate effects (GMM) ###################

        # GMM

        # pre-compute some stuff that is the same for all snps
        zerov = np.zeros(K)
        mvn_0L = [multivariate_normal(mean=zerov, cov=Omega_0 + Sigma_j) for Sigma_j in Sigma_jL]
        mvn_1L = [multivariate_normal(mean=zerov, cov=Omega_1 + Sigma_j) for Sigma_j in Sigma_jL]
        AV_0L, AV_1L = get_gmm_wts(Omega_0, Omega_1, Sigma_jL)

        beta_gmms = np.zeros((M,K))

        print("")
        # TODO could probably vectorize and do in batches
        for j in range(M):
            if j % 100000 == 0:
                print("Estimating MT-Omega_0%s on snp %d..." % ("abcd"[Omega_0_index], j))
            Sigma_j = Sigma_jL[j]
            mvn_0 = mvn_0L[j]
            mvn_1 = mvn_1L[j]
            AV_0 = AV_0L[j]
            AV_1 = AV_1L[j]
            beta_gmms[j,:] = get_beta_gmm(betahats[j,:], 
                Omega_0, Omega_1, Sigma_j, 
                pi_0_est, pi_1_est,
                mvn_0, mvn_1, AV_0, AV_1)

        newStats[:,(Omega_0_index + 1)] = beta_gmms[:,t0]
        print("estimated gmm effects for ", Omega_0_name)

       # TODO report GMM estimates for other traits somewhere


    # Note: the p-value for MTAG or MT is used only for clumping to prioritize which SNPs to keep, so just use GWAS SE
    # TODO appropriate MTAG or MT model probably models log(OR) instead of OR, e.g.
    #       newZscores[i] = 1.0*log(newOR_t0)/newSD
    # note: for OR, take log to get z-score


    # compute fake p-value using each snp's gwas SE, for each method
    for i in range(M):
        snp = snpIDs[i]
        gwasLineL = snpInfoD[snp] # includes gwas beta and pval
        SE = gwasLineL[plinkOutFields.index('SE')] # TODO get index of SE from header]
        SE = float(SE)


        if isCaseCon:  # take log of OR to get z-score for each method
            print("caseCon not supported yet TODO")
            newZscores_i = 1.0*np.log(newStats[i,:])/SE
            
        else:
            newZscores_i = newStats[i,:]/SE
            
        # two-sided univariate test for zscore
        newPvals[i, :] = 2*norm.cdf(-1.0*abs(newZscores_i))


    # write each estimate (mtag and all mt versions) in plink format for clumping, as additional columns in a gwas plink file


    # TODO write something else for caseCOn header
    betaCols = ['BETA_MTAG'] + ['BETA_GMM' + 'abcd'[i] for i in range(4)]
    pvalCols = ['PVAL_MTAG'] + ['PVAL_GMM' + 'abcd'[i] for i in range(4)]


    # write estimates for MTAG and MT-Omega_0x with each M_eff
    gmmOut = outdir + "/" + exptName + ".assoc.MTAG.MT_estimates"
    with open(gmmOut, 'w') as outf:
        outf.write('\t'.join(plinkOutFields) + '\t') # write header
        outf.write('\t'.join(betaCols) + '\t')
        outf.write('\t'.join(pvalCols) + '\n') 
        for i in range(M):
            
            # replace BETA/OR, STATS, and P field with new value (or append to end or line if not already in the field)
            snpName = snpIDs[i]
            outf.write('\t'.join(snpInfoD[snpName]) + '\t')
            outf.write('\t'.join(str(val) for val in newStats[i,:]) + '\t')
            outf.write('\t'.join(str(val) for val in newPvals[i,:]) + '\n')

            if i%50000 == 0:
                print("Writing snp %d..." % i)

  

