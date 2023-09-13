#!/bin/bash
#$ -cwd 
#$ -o /myjob
#$ -N MT_job_example
#$ -l h_data=4G,h_rt=1:00:00,highp
#$ -M MyName
#$ -m bea
#$ -t 1:8:1 # Depending on the number of traits you would like to analyze

source /etc/profile
module load R/3.6.0

PHENONAME="`sed -n ${SGE_TASK_ID}p PHENO_NAME_file`" # Takes the phenotype names one by one. File example see Pheno_list_example.txt
TRAIT_SET="`sed -n ${SGE_TASK_ID}p SET_NAME_file`" # Takes the trait set names one by one. Usually named by the initials of the traits in the set in multi-trait analysis. The total number of lines in the trait set name file should match the total number of traits. For example, armfat, trunkfat and waistc will have a set name of atw. The trait set name file should have atw repeated in 3 lines corresponding to these three traits. Example see Trait_set_name_example.txt
PHENO="`sed -n ${SGE_TASK_ID}p PHENO_FILE_NAME_file`" # correspond to phenotype file name prefixes (excluding .pheno). 

# Details of the function check generate_r2_df.R

echo ${PHENONAME}

Rscript --no-save --no-restore --verbose generate_r2_df.R path/Covariance_example.covar ${PHENONAME} path/${PHENO}.pheno ldpred_${PHENONAME}_${TRAIT_SET}_result.txt ${TRAIT_SET} output_path > DIRECTORY_TO_ROUT_FILES/${PHENONAME}.generate.df.Rout 2>&1

