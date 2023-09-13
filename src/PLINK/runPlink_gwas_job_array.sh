#!/bin/bash
#$ -cwd
#$ -V
#$ -N gwas
#$ -l h_data=16G,h_rt=36:59:59,highp
#$ -t 1:18:1 # TODO diff number of caseCon or quant traits

# note:
# -o: what directory to write any printed output or error messages to
# -t: makes it a job array. change 1:14:1 to however many traits you are running

source /etc/profile     # so module command is recognized
module load plink

# The above is an example for running this code in UCLA Hoffman2 Job submission system, the user may change the code for different job submission systems. 

# PHENO_file is a text file containing the phenotype names with one name per line (the names should correspond to the input phenotype file names). Example see Pheno_list_example.txt
# PHENO_NAME_file is a text file containing the phenotype names with one name per line (the names should correspond to the phenotype names the user would like for the output GWAS results)
# BFILE is the prefix of a list of files with suffix .bed, .bim, .fam and .log, which should be the bfiles generated during previous sub-steps.
# PHENOFILE is the directory to the phenotype files that the user would input into the plink script. The phenotype files should have names that are the same as listed in PHENO_file, and should have suffix ".pheno". Example file see: Phenotype_example.pheno. The phenotype files should have the first column as FID, second column as IID, and third column as phenotype values.
# PLINKOUTFILE is the directory to the output files that the script will produce. The output file names will be the phenotype names specified in PHENO_NAME_file followed by a suffix that the user specifies (in this case it's ".rand200k_indivs")
# COVFILE is the directory to covariance files of the individuals that the study is conducted on, the covariance includes sex and first 10 pcs in our experiment. Example file see: Covariance_example.covar
# 

# get the phenotype name from a text file with one name per line
PHENO="`sed -n ${SGE_TASK_ID}p PHENO_file`"

PHENO_NAME="`sed -n ${SGE_TASK_ID}p PHENO_NAME_file`"

BFILE="PATH_TO_BFILE"

PHENOFILE="PATH_TO_PHENOFILE/${PHENO}.pheno"

PLINKOUTFILE="PATH_TO_PLINKOUTFILE/${PHENO_NAME}.rand200k_indivs"

# currently sex/first 10 PCs
COVFILE="PATH_TO_COVFILE/Covariance_example.covar"


## Below is the script used for running PLINK that takes in one phenotype at a time and runs PLINK respectively.

echo $PHENO
date

# run plink linear or logistic regression

# linear
# (standard-beta to standardize phenotypes when estimating beta) 
# (--allow-no-sex flag to use indivs with no sex, or no sex in the .fam file)
# (--ci 0.95)  95% confidence upper and lower bound, plus SE
plink --linear standard-beta hide-covar --allow-no-sex --covar $COVFILE --bfile $BFILE --pheno $PHENOFILE --out ${PLINKOUTFILE} --ci 0.95 


echo $PHENO # makes it easy to check which trait unfinished
date







