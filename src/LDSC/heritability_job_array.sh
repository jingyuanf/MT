#!/bin/bash
#$ -cwd 
#$ -o myjob/
#$ -N Heritability
#$ -l h_data=4G,h_rt=01:00:00,highp
#$ -M myname/
#$ -m bea
#$ -t 1-8:1 # Change 1:8-1 according to the number of traits that the user is analysing

# PHENO_NAME_file is a text file containing the phenotype names with one name per line (the names should correspond to the phenotype names of the output GWAS results)

source /etc/profile     # so module command is recognized
module load anaconda
source activate ldsc

PHENO1="`sed -n ${SGE_TASK_ID}p PHENO_NAME_file`" # get first phenotype line by line

echo ${PHENO1}

for i in {1..8} # loop through the number of traits to get trait 2, this number should be changed according to the actual number of traits
do

	PHENO2='`sed -n ${i}p PHENO_NAME_file`' # get second phenotype
	echo ${PHENO2}

	./ldsc.py \
 	--rg PATH_TO_TRAIT_1/${PHENO1}.sumstats.gz,PATH_TO_TRAIT_2/${PHENO2}.sumstats.gz \
 	--ref-ld /u/scratch/f/fujy2038/BIG/LDSC/ldsc_val10k_withdrawn_filtered_results/val10k \
	 --w-ld /u/scratch/f/fujy2038/BIG/LDSC/ldsc_val10k_withdrawn_filtered_results/val10k \
	 --out /u/scratch/f/fujy2038/BIG/LDSC/heritability/200k_withdraw/${PHENO1}_${PHENO2}

done

