#!/bin/bash
#$ -cwd 
#$ -o /myjob
#$ -N MT_job_example
#$ -l h_data=4G,h_rt=1:00:00,highp
#$ -M MyName
#$ -m bea

source /etc/profile
module load R/3.6.0

# This script generates barplots for a specific set of analysis. For example, armfat_percent, trunkfat_percent and waistc. It generates an incremental R^2 plot as well as an incremental adjusted R^2 plot.

Rscript --no-save --no-restore --verbose barplots.R path/rds_filenames.txt path/phenotype_names.txt path/colors.txt atw output_path > DIRECTORY_TO_ROUT_FILES/atw.barplot.Rout 2>&1

