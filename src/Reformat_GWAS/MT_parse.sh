#!/bin/bash
#$ -cwd 
#$ -o myjob/
#$ -N Parse
#$ -l h_data=4G,h_rt=01:00:00,highp
#$ -M fujy2038
#$ -m bea
#$ -t 1-8:1 

# The above is an example for running this code in UCLA Hoffman2 Job submission system, the user may change the code for different job submission systems. 


source /etc/profile     # so module command is recognized
module load python/3.7.2

# PHENO: get the phenotype name from a text file with one name per line
# the phenotype names in the lines for GWAS summary statistics are the trait names listed in pheno_name_list that correspond to GWAS summary statistics names. Example see Pheno_list_example.txt.


PHENO="`sed -n ${SGE_TASK_ID}p PATH_TO_GWAS_list/pheno_name_list.txt'"

echo $PHENO

python3 /u/scratch/f/fujy2038/BIG/BIG_code_documentation/MT/GWAS_reformat_code/MT_parse_MT.py\
 -d ${PHENO}.rand200k_indivs.assoc.linear_A2\
 -o /u/scratch/f/fujy2038/BIG/ukbb_formatted/gwas_ss_plink_wSE_w_filtered_A2/${PHENO}.rand200k_indivs.assoc.linear_A2.formatted\
 --ignore TEST,L95,U95 --snp SNP --con n --beta BETA --z STAT --n NMISS




