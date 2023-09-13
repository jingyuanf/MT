#!/bin/bash
#$ -cwd
#$ -V
#$ -N gwas
#$ -l h_data=16G,h_rt=4:59:59,highp

source /etc/profile     # so module command is recognized
module load plink

# BFILE is the prefix of a list of files with suffix .bed, .bim, .fam and .log, which should be the bfiles generated during previous sub-steps.

BFILE="PATH_TO_BFILE"

plink --bfile $BFILE --freq

