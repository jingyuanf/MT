#!/bin/bash
#$ -cwd
#$ -V
#$ -N generateMAF
#$ -l h_rt=1:00:00,h_data=16G

source /etc/profile
module load plink

BFILE="/u/scratch/c/cmhuang/data_ukbb_other_traits/ukbb_150k_indiv_for_gwas/filter4_rand150k_indivs"
OUTFILE="/u/scratch/c/cmhuang/data_ukbb_other_traits/plink150k"

plink --bfile $BFILE --freq --out $OUTFILE
