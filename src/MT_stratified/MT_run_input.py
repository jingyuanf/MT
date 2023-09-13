import traceback # for exceptions such as inverting non-singular matrix
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

outdir = ''
traitL = []
traitPostfix = ''
gwasPath = ''


parser = argparse.ArgumentParser()

parser.add_argument('-g', '--gwasdir', required=True, dest='gwasdir', help='Folder path containing input gwas data')
parser.add_argument('-t', '--traitnames', required=True, dest='traitnames', help='input trait names of interest, separated by a comma. The first trait name is the default main trait.')
parser.add_argument('-sf', '--suffix', default='', help='suffix of input GWAS file names')
parser.add_argument('--gcov', required=True, help='Directory of a list of file names containing genetic covariance from LDSC results')
parser.add_argument('--gcovpath', required=True, help='Path to the files containing genetic covariance from LDSC results')
parser.add_argument('--h2', required=True, help='Directory of a list of file names containing estimation error covariance from LDSC results')
parser.add_argument('--h2path', required=True, help='Path to the files containing estimation error covariance from LDSC results')
parser.add_argument('-o', '--outdir', required=True, dest='outdir',help='output directory')
parser.add_argument('--binpath', help='Path to the bins')
parser.add_argument('--bin', help='Directory of a list of ldsc bins file names')

args = parser.parse_args()

#settingsL = [m for m in (int(5e5), int(1/(5*10**-8)))]
M_effective = 51945
t0 = 0

traitL = args.traitnames.split(',')
gwasdir = args.gwasdir
traitPostfix = args.suffix
outdir = args.outdir
gcovfile = args.gcov
h2file = args.h2
binfile = args.bin
gcovpath = args.gcovpath
h2path = args.h2path
binpath = args.binpath


'''
for M_effective in settingsL:
    if M_effective == 51945:
        Meff_str = '5e5'
    elif M_effective == 2e7:
        Meff_str = '2e7'
    else:
        Meff_str = 'TODO_other'
'''

exptName = traitL[0] + '_' + str(len(traitL)) + ''.join([name[0] for name in traitL[1:]])
gwasNamesL = [gwasdir + t + traitPostfix for t in traitL]
print(exptName)

gcov = []
h2 = []
bins = []

with open(gcovfile) as f:
    for line in f:
        file = line.rstrip()
        gcov.append(gcovpath + file)

with open(h2file) as f:
    for line in f:
        file = line.rstrip()
        h2.append(h2path + file)

with open(binfile) as f:
    for line in f:
        file = line.rstrip()
        bins.append(binpath + file)
    

try:
    run_MTAG_MT(exptName, traitL, gwasNamesL, t0, M_effective, outdir, gcov, h2, bins)
        
# note, this skips it if any of Omegaa-d results in a singular matrix when computing AV_0 or AV_1
except np.linalg.LinAlgError as err:
    print("ERROR: Skipping ", exptName)
    print("ERROR: Encountered linAlgError:\n", err)
    print("Printing traceback")
    traceback.print_exc()

    print("\n\n")





