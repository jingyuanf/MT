#!/bin/bash
#$ -cwd 
#$ -o /u/scratch/f/fujy2038/BIG/ukbb_output/job/
#$ -N mt_mt_jy_apr24
#$ -l h_data=4G,h_rt=01:00:00,highp
#$ -M fujy2038
#$ -m bea
#$ -t 1-4:1 

# PHENO_NAME_file is a text file containing the phenotype names with one name per line (the names should correspond to the phenotype names of the output GWAS results)

source /etc/profile     # so module command is recognized
module load anaconda
source activate ldsc

traitL="`sed -n ${SGE_TASK_ID}p PHENO_NAME_file`"

echo $SGE_TASK_ID

./munge_sumstats.py --sumstats ./PATH_TO_GWAS/${traitL[$SGE_TASK_ID-1]}.sumstats --out ./PATH_TO_GWAS_MUNGED/${traitL[$SGE_TASK_ID-1]} --N-col NMISS


