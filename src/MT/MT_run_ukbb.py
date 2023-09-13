import traceback # for exceptions such as inverting non-singular matrix
import sys
sys.path.insert(0, '/u/scratch/l/lgai/MultiTrait_UKBB/code/')

import itertools
#import time
import random
from math import sqrt
import numpy as np
import argparse
import time
import os, re
import sys, gzip, bz2
import logging
from argparse import Namespace

# ensure it can still be run from outside the code folder

# TODO add a log file that includes info like if any singular matrices were encountered, or non-PD cov converted to nearest PD

from MT_run_function import * 
# contains a function to get the genetic cov, non-genetic cov, read in gwas data, run mtag, and run gmm on 4 different types of Omega 0's

# Replace --stat BETA with BETA_MTAG, BETA_GMMa, etc


##### Set output file locations ######


#exptName = 'armfat7ww_Meff_2e7'





#M_effective = int(1/(5*10**-8)) # based on gwas sig threshold
#M_effective = 500000.0  # assume M_eff roughly equal to num snps in study/analysis
# M_effective = 500000.0 #/4.0 # assume every 4th snp in data is roughly indept/ if 1/5e-8 indep snps out of 5 mil in human genome.
#M is either approx num indept snps in genome, number of snps in filtered SS list, not in .sorted (may contain snps not in the bfile)

# TODO get rid of this
#M_prealloc = 500000 # upper bound on number of snps file, used for preallocating array

'''
traitLL = []

traitL = ["armfat_percent", 
"bmi21001",
"height",
"height_seated",
"hipc",
"trunkfat_percent",
"weight21001",
"weight23098"]
traitLL.append(traitL)

traitL = ["armfat_percent", 
"bmi21001",
"height",
"height_seated",
"trunkfat_percent",
"weight21001"]
traitLL.append(traitL)

traitL = ["armfat_percent", "weight21001", "weight23098"]
traitLL.append(traitL)


traitL = ["diabp",
"pulse_automated",
"pulse_manual",
"sysbp"]
traitLL.append(traitL)

traitL = ["pulse_automated",
"diabp",
"pulse_manual",
"sysbp"]
traitLL.append(traitL)
'''

traitLL = [] # TODO
traitL = ["waistc",
"trunkfat_percent",
"armfat_percent"]
traitLL.append(traitL)

'''
traitL = ["pulse_automated",
"pulse_manual"]
traitLL.append(traitL)

traitL = ["armfat_percent",
"sysbp"]
traitLL.append(traitL)
'''

# TODO can mention in discussion that some traits may be overweighted, e.g. measuring arm_fat with aux traits all some variation of weight, redundant aux traits





settingsL = [(t,m) for t in traitLL for m in (int(5e5), int(1/(5*10**-8)) )]


t0 = 0 # index of main trait. could have as sep command line options for primary and aux traits


for traitL, M_effective in settingsL:
    if M_effective == 5e5:
        Meff_str = '5e5'
    elif M_effective == 2e7:
        Meff_str = '2e7'
    else:
        Meff_str = 'TODO_other'


    exptName = traitL[0] + '_' + str(len(traitL)) + ''.join([name[0] for name in traitL[1:]]) + '_Meff_' + Meff_str 
    # e.g. 'armfat_7ww_Meff_2e7'
    gwasNamesL = ['/u/scratch/f/fujy2038/BIG/ukbb_gwas_ss_plink_wSE_linear_150k_A2_MAF_rsID/' + t + ".rand150k_indivs.assoc.linear_150k_A2_MAF_150k_A2" for t in traitL] # currently only using linear traits
    print(exptName)

    try:
        run_MTAG_MT(exptName, traitL, gwasNamesL, t0, M_effective)
        
    # note, this skips it if any of Omegaa-d results in a singular matrix when computing AV_0 or AV_1
    except np.linalg.LinAlgError as err:
        print("ERROR: Skipping ", exptName)
        print("ERROR: Encountered linAlgError:\n", err)
        print("Printing traceback")
        traceback.print_exc()

    print("\n\n")


'''
Rscript PRSice.R --dir . --prsice PRSice_linux --base

Rscript ./PRSice_mar17/PRSice.R --dir ./PRSice_mar18 \
    --prsice ./PRSice_mar17/PRSice_linux \
    --base ./gwas_ss_plink_wSE/armfat_percent.rand150k_indivs.assoc.linear \
    --target ./data/filter4_validation_1k_indivs \
    --thread 1 \
    --stat BETA \
    --beta \
    --binary-target F \
    --pheno-file ./data_pheno/armfat_percent.validation_1k_indivs.pheno \
    --cov-file ./data/ukb21970_validation_1k_indivs.hicho.mt.covar \
    --out armfat_percent.gwas

'''

