#!/bin/bash
#$ -cwd 
#$ -V
#$ -N PLINK_filterIndivs
#$ -l h_data=16G,h_rt=23:59:59,highp

source /etc/profile     # so module command is recognized
module load plink

# The above is an example for running this code in UCLA Hoffman2 Job submission system, the user may change the code for different job submission systems. 

# BFILE_ORIG should be a folder containing the original bed, bim, fam, and log files (with prefix of the files).
# ".bed" files contain binary encoded genetic information
# ".bim" files contain SNP names and map positions
# ".fam" files contain family information
# The script will output a new set of bed, bim, fam, and log files for the filtered individuals, in the directory specified in OUTBASE.
# OUTBASE is the directory and prefix of output.
# KEEPFILE specifies the list of individuals that the PLINK is used for analysis, with the columns being Family ID (FID) and Individual ID (IID). Example see: Ind_file_example.ind

BFILE_ORIG="TODO"
KEEPFILE="TODO"
OUTBASE="TODO"

plink --bfile $BFILE_ORIG --keep $KEEPFILE --make-bed --out $OUTBASE

echo "\n\n"


