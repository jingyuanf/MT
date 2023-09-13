README

### relevant files ###

## Main Script
MT_run_input.py - script to run effect size estimations on traits. Use -h to check functional handles it takes for input.

## Job Submission Script
MT_run_bash.sh - example bash file for submitting jobs


## Code called in the Main Script
MT_functions_cov.py - get genetic covariance and estimation standard error correlation matrix entries from text file summarizing ldsc output
MT_functions_effects.py - functions for estimating effect sizes on single trait using multiple traits. It includes MTAG and mixture model estimates
MT_functions_gmm.py - implements functions for gaussian mixture model
MT_functions_gwas_bins.py - functions for parsing gwas summary statistics that are split into bins. Currently only parses plink files or other specific format
MT_run_function.py - function for running multiple steps of effect size estimates for a set of traits


## Input Data
Trait_file_example.txt - example text file that contains a list of trait sets for analysis. Each trait set should be in different lines and each trait in the trait set should be separated by commas.
Heritability_file_example.txt - example text file that contains the genetic heritability of each trait, generated from ldsc
Genetic_covariance_file_example.txt - example text file that contains the genetic covariance and genetic correlation between each pair of traits, generated from ldsc


