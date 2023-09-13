#!/bin/bash
#$ -cwd 
#$ -o /myjob
#$ -N MT_job_example_4bins
#$ -l h_data=16G,h_rt=47:59:59,highp
#$ -M MyName
#$ -m bea

source /etc/profile     # so module command is recognized
module load python/3.7.2

# The above is an example for running this code UCLA Hoffman2 Job submission system, the user may change the code for different job submission systems.

# Below is a for loop submission script that takes in a file containing all trait sets that need analysis, echos the command line to a separate bash file for each trait set, and submits the bash files one by one.

# DIRECTORY_OF_GWAS_FILES_FOLDER: GWAS file generated from PLINK and formatted using MT_parse_MT.py. Example see GWAS_example.rand200k.
# SUFFIX_OF_GWAS_FILES: The suffix of the input GWAS file, in this example it is rand200k.
# GENETIC_COVARIANCE_FILENAMES: File containing the list of genetic covariance file names. Genetic covariance is computed using ldsc and parsed using parse_ldsc.py. Example of genetic covariance file see Genetic_covariance_file_example.txt
# HERITABILITY_FILENAMES: File containing the list of heritability file names. Heritability is computed using ldsc and parsed using parse_ldsc.py. Example see Heritability_file_example.txt
# LDSC_FILENAMES: File containing the list of four LD Score file names (as LD scores are splitted into 4 bins). Example of LD Score files see LDSC_file_example.txt
# PATH_TO_OUTPUT_DIRECTORY: Directory for the outputting MT results
# PATH_TO_BASH_FILE: The location the user would like the bash file to be.
# TRAIT_FILE: The file containing the trait sets that the user would like to analyze. Example see Trait_file_example.txt

#traitL=("armfat_percent,trunkfat_percent,waistc" "trunkfat_percent,armfat_percent,waistc" "waistc,trunkfat_percent,armfat_percent" "pulse_automated,pulse_manual"  "pulse_manual,pulse_automated" "pulse_automated,height" "height,pulse_automated")

#traitL=("height,pulse_automated")


while read trait; do

	echo "python3 MT_run_input.py -g DIRECTORY_OF_GWAS_FILES_FOLDER -t ${trait} -sf SUFFIX_OF_GWAS_FILES --gcov GENETIC_COVARIANCE_FILENAMES --gcovpath PATH_TO_GENETIC_COVARIANCE_FILES --h2 HERITABILITY_FILENAMES --h2path PATH_TO_HERITABILITY_FILES --bin LDSC_FILENAMES --binpath PATH_TO_LDSC_FILES --o PATH_TO_OUTPUT_DIRECTORY" > PATH_TO_BASH_FILE/${trait}_job_bins.sh

qsub PATH_TO_BASH_FILE/${trait}_job_bins.sh

done < PATH_TO_TRAIT_FILE

