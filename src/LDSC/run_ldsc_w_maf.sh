#!/bin/bash
#$ -cwd
#$ -V
#$ -N ldsc
#$ -l h_data=16G,h_rt=23:59:59,highp

source /etc/profile

module load anaconda
source activate ldsc

## BFILE is the prefix of a list of files with suffix .bed, .bim, .fam and .log.

BFILE='PATH_TO_BFILES'
OUTFILE='OUTPUT_DIRECTORY'

python ldsc/ldsc_keep_maf.py --bfile $BFILE --l2 --ld-wind-kb 2000 --out $OUTFILE
