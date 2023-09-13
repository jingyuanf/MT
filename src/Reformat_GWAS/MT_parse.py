#!/usr/bin/env python
# coding: utf-8


import itertools
import random
from math import sqrt
import numpy as np
import argparse
import time
import os, re
import sys, gzip, bz2
import logging
from argparse import Namespace

epsilon = 1e-12 # the small non-zero value

np.set_printoptions(precision=3)


fields = []
parser = argparse.ArgumentParser()
parser.add_argument('--snp', required= True, help='header name for SNPs, required')
parser.add_argument('--n', help='header name for effective sample size')
parser.add_argument('--n_val', help='integer number for effective sample size, uniformly distributed')
parser.add_argument('--chr', help = 'header name for chromosome')
parser.add_argument('--bp', help = 'header name for base position')
parser.add_argument('--a1', help = 'header name for allele 1')
parser.add_argument('--a2', help = 'header name for allele 2')
parser.add_argument('--chrbp', help = 'header name for chromosome and base position if they are combined in one column')
parser.add_argument('--chrbpsplit', default= ':', help = 'splitting symbol between chromosome and base position. Default is colon.')
parser.add_argument('--savechrbp', help = 'header name for chromosome and base position if they are combined in one column and you want to save this column, cannot be specified simultaneously with chrbp')
parser.add_argument('--no_chr_data', help = 'no chromosome data')
parser.add_argument('--beta', default= None, help='header name for estimated effect sizes')
parser.add_argument('--or', default = None, help = 'header name for odds ratio')
parser.add_argument('--se', default = None, help ='header name for standard error')
parser.add_argument('--con', required= True, help='y/n for case control')
parser.add_argument('--z', default= None, help='header name for z-score')
parser.add_argument('--op', default= '',help='any other header names of interest, separated by a comma')
parser.add_argument('--p', default = None, help = 'p-value')
parser.add_argument('--ignore', default = '', help = 'any header names to ignore, separated by a comma')

parser.add_argument('-d', '--data', required=True, dest='data_file', help='File containing data.')
parser.add_argument('-o', '--outputFile', required=True, dest='output_file',help='Output file name.')

args = parser.parse_args()

if args.con=='y':
    isCaseCon = True
elif args.con=='n':
    isCaseCon = False
else:
    print("ERROR: for case control, must be 'y' or 'n'")
    exit(1)

inputdir = args.data_file
outdir = args.output_file

snp = args.snp
n = args.n
n_val = args.n_val
chrm = args.chr
bp = args.bp
chrbp = args.chrbp
chrbpsplit = args.chrbpsplit
savechrbp = args.savechrbp
beta = args.beta
se = args.se
z = args.z
p = args.p
a1 = args.a1
a2 = args.a2
op = args.op.split(',')
ignore = args.ignore.split(',')

use_beta_se = False
use_chrbp = False

if not (beta and se) and not z:
    print("ERROR: Beta and Standard Error must be provided if Z score is not provided.")
    exit(1)
elif not z:
    use_beta_se = True
    
if not n and not n_val:
    print("ERROR: N column header or N values must be specified.")
    exit(1)

if chrbp:
    use_chrbp = True

header_map = {}

if snp:
    header_map[snp] = 'SNP'
if n:
    header_map[n] = 'NMISS'
if beta:
    header_map[beta] = 'BETA'
if z:
    header_map[z] = 'STAT'
if p:
    header_map[p] = 'P'
if a1:
    header_map[a1] = 'A1'
if a2:
    header_map[a2] = 'A2'
if savechrbp:
    header_map[savechrbp] = 'snpid'

print(fields)
print(header_map)


with open(inputdir) as f:
    headerL = f.readline().rstrip().split()
    print(headerL)

## A list of column names that could be inferred to be converted to required column names. The user may modify at their own need.

inferHeaderDict = {
    'SNP': 'SNP',
    'snp': 'SNP',
    'MARKERNAME': 'SNP',
    'markername': 'SNP',
    'RS':  'SNP',
    'RSID':  'SNP',
    'RS_NUMBER':  'SNP',
    'RS_NUMBERS':  'SNP',
    #'rsID': 'SNP',
    #'snpid':'SNP',

    # n
    'N': 'NMISS',
    'n': 'NMISS',
    'sample_size': 'NMISS',
    'ncol': 'NMISS',

    # freq
    'FREQ': 'MAF',
    'A1FREQ': 'MAF',
    'a1freq': 'MAF',
    'EAF': 'MAF',
    'eaf': 'MAF',
    'FRQ': 'MAF',
    'frq': 'MAF',
    'AF': 'MAF',
    'FRQ': 'MAF',
    'MAF': 'MAF',
    'FRQ_U': 'MAF',
    'F_U': 'MAF',
    'freq': 'MAF',

    # chr
    'CHR': 'CHR',
    'Chromosome': 'CHR',
    'chromosome': 'CHR',
    'Chr': 'CHR',
    'chr': 'CHR',

    # bpos
    'BPOS': 'BP',
    'Bpos': 'BP',
    'BP': 'BP',
    'bp': 'BP',
    'POS': 'BP',
    'Pos': 'BP',
    'pos': 'BP',
    'position': 'BP',
    'Position': 'BP',
    'bpos': 'BP',

    # a1
    'A1': 'A1',
    'ALLELE1': 'A1',
    'allele1': 'A1',
    'EFFECT_ALLELE': 'A1',
    'effect_allele': 'A1',
    'EA': 'A1',
    'ea': 'A1',
    'a1': 'A1',

    # a2
    'A2': 'A2',
    'ALLELE0': 'A2',
    'allele0': 'A2',
    'ALLELE2': 'A2',
    'allele2': 'A2',
    'OTHER_ALLELE': 'A2',
    'other_allele': 'A2',
    'OA': 'A2',
    'oa': 'A2',
    'a2': 'A2',

    # beta
    'BETA': 'BETA',
    'Beta': 'BETA',
    'EFFECT': 'BETA',
    'Effect': 'BETA',
    'effect': 'BETA',
    'b': 'BETA',
    'beta': 'BETA',

    # se
    'SE': 'SE',
    'SE_unadj': 'SE',
    'se_unadj': 'SE',
    'SE_UNADJ': 'SE',
    'se': 'SE',
    's': 'SE',  

    # z
    'Z':'STAT',
    'STAT': 'STAT',
    'Z_unadj': 'STAT',
    'z_unadj': 'STAT',
    'Z_UNADJ': 'STAT',
    'z': 'STAT',
    'Z-score':'STAT',
    'z-score':'STAT',
    'ZSCORE':'STAT',

    # pval
    'PVAL': 'P',
    'Pval': 'P',
    'P_BOLT_LMM_INF': 'P',
    'P_BOLT_LMM': 'P',
    'P': 'P',
    'p': 'P',
    'P_unadj': 'P',
    'p_unadj': 'P',
    'P_UNADJ': 'P',
    'pval': 'P',

    # info
    'INFO': 'INFO',
    'info': 'INFO',
    'RSQ': 'INFO',
    'rsq': 'INFO',
}

with open(inputdir) as f:
        
        headerL = f.readline().rstrip().split()
        fieldIndices = {}
        inferHeaderMap = {}
        for col in headerL:
            fieldIndices[col] = headerL.index(col)
        #inferHeaderMap[col] = transform_header(col)
        #snpIndex = headerL.index('SNP')
        
        if use_beta_se:
            headerL.append('STAT')
        if not n:
            headerL.append('N')
        if chrbp:
            headerL.append('CHR')
            headerL.append('BP')
            
        with open(outdir, 'w') as g:

            toWrite = ''
            for header in headerL:
                if header in ignore:
                    continue
                if header in header_map:
                    toWrite += header_map[header]
                elif header in inferHeaderDict:
                    toWrite += inferHeaderDict[header]
                else:
                    toWrite += header
                toWrite += '\t'
            toWrite += '\n'

            g.write(toWrite)
            toWrite = ''
            for line in f:
                toWrite = line.rstrip() + '\t'
                lineL = toWrite.split()
                toWrite = ''
                for i in range(len(lineL)):
                    if headerL[i] not in ignore:
                        toWrite += lineL[i] + '\t'
                if use_beta_se:
                    z = float(lineL[fieldIndices[beta]])/float(lineL[fieldIndices[se]])
                    toWrite += str(z) + '\t'
                if not n:
                    toWrite += str(n_val) + '\t'
                #write to file...
                if chrbp:
                    cb = lineL[fieldIndices[chrbp]].split(chrbpsplit)
                    for i in cb:
                        toWrite += i + '\t'
                toWrite += '\n'
                g.write(toWrite)
