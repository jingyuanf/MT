#TODO also get all_phe_cor number, as using the sample_size_pairwise.txt
#seems to be N1*N2/N_shared

import argparse
from argparse import Namespace

parser = argparse.ArgumentParser()
parser.add_argument('--loglist', required = True, help = 'A file containing the list of the log file names generated from ldsc. Example see Log_file_names_example.txt')
parser.add_argument('--logpath', required = True, help = 'The directory containing the log files listed in argument --loglist')
parser.add_argument('--gcov', required = True, help = 'The directory and file name specified by the user for genetic covariance. Example: mypath/ldsc_gcov/ldsc_results_gcov.txt')
parser.add_argument('--h2', required = True, help = 'The directory and file name specified by the user for heritability. Example : Example: mypath/ldsc_h2/ldsc_results_h2.txt')


args = parser.parse_args()

listOfLogFiles_F = args.loglist
logFilePath = args.logpath
outName = args.gcov
outName2 = args.h2



def fileName2Trait(fileName):
    '''get my_trait from my_trait's path''
    fileName = fileName.split("/")[-1] # remove path
    return fileName.split(".")[0] # return trait nly


def get_val_and_se(line):
    '''given str "some text: 0.123 (0.567)", return ("0.123", "0.567")'''
    temp = line.split(": ")[1]  # "0.123 (0.567)"
    val = temp.split()[0]
    se = temp.split()[1].lstrip("(").rstrip(")")
    if val == "nan": 
        # e.g. Genetic Correlation: nan (nan) (h2  out of bounds) 
        val = "NA"
        se = "NA"
    return (val, se)

def get_trait_LDSC_results(lines, i):
    # Call after seeing "Heritability of phenotype X". next line is dashes, then h2 + SE, gc factor, chi2, intercept + SE, ratio

    h2, h2_se = get_val_and_se(lines[i+2])

    lambda_GC = lines[i+3].split(": ")[1]
    mean_chisq = lines[i+4].split(": ")[1]
    intercept, intercept_se = get_val_and_se(lines[i+5])
    
    if lines[i+6] == "Ratio < 0 (usually indicates GC correction).":
        ratio = "-9" # want to indicate it is negative, could also encode as nan or NA
        ratio_se = "NA"
    
    else:
        ratio, ratio_se = get_val_and_se(lines[i+6])

    return (h2, h2_se, lambda_GC, mean_chisq, intercept, intercept_se, ratio)




with open(listOfLogFiles_F) as f:
    logFileL = [l.rstrip() for l in f.readlines()]
    if logFileL[0] == "_AND_.log": # sometimes run with no traits, artifact of qsub job arr
        logFileL = logFileL[1:]

logFileL = [logFilePath + l for l in logFileL]



# write gcorr results as we read the log for each trait pair
# use a dict to store h2, etc since traits will appear in multiple pairs
traitResultsD = {} # traitResultsD[trait] = (h2, h2_se, ...)
traitPairResultsD = {} # traitResultsD[(trait1, trait2)]
trait_nSnpD = {}


# write the results for the trait pair to table
outf = open(outName,'w')
outf.write("\t".join(["trait1", "trait2", "gcov", "gcovSE", "meanZ1Z2", "gcovIntercept", "gcovInterceptSE", "gcorr", "gcorrSE", "gcorrZ", "gcorrP"]) + "\n")

# go through log file for each trait pair
for logFile in logFileL:
    print(logFile)
    with open(logFile) as logf:
        lines = [line.rstrip() for line in logf.readlines()]
    
    gotTraitName1 = False

    for i in range(len(lines)):
        line = lines[i]
        if line[:21] == "Beginning analysis at":
            # next two lines have file name and num snps for first trait
            traitName1 = fileName2Trait(lines[i+1].split()[-2]) # ignore ellipses at end
            nSnp_trait1 = lines[i+2].split()[-2]
            gotTraitName1 = True

        if line[:31] == "Reading summary statistics from" and gotTraitName1:
            traitName2 = fileName2Trait(lines[i].split()[-2])

        elif line.split()[1:] == ["SNPs", "with", "valid", "alleles."]:
            nSnp_used = line[0]

        elif line == "Heritability of phenotype 1":
            trait1_results = get_trait_LDSC_results(lines, i)

        elif line == "Heritability of phenotype 2" or line[:28] == "Heritability of phenotype 2/": # e.g. phenotype 2/2
            trait2_results = get_trait_LDSC_results(lines, i)

        elif line == "Genetic Covariance": # observed scale
            gcov, gcov_se = get_val_and_se(lines[i+2])
            mean_z1z2 = lines[i+3].split(": ")[1]
            gcov_intercept, gcov_intercept_se = get_val_and_se(lines[i+4])

        elif line == "Genetic Correlation": 
            gcorr, gcorr_se = get_val_and_se(lines[i+2])
            if gcorr == "NA":
                # if gcorr reads "nan (nan) (h2  out of bounds)"
                gcorr_zscore = "NA"
                gcorr_pval = "NA"
            else:
                gcorr_zscore = lines[i+3].split(": ")[1]
                gcorr_pval = lines[i+4].split(": ")[1]
            

        # else continue

    # update h2 and other trait-specific results. traits will be entered multiple times
    traitResultsD[traitName1] = trait1_results
    traitResultsD[traitName2] = trait2_results
    trait_nSnpD[traitName1] = nSnp_trait1

    # write gcorr results to file
    outf.write("\t".join([traitName1, traitName2, gcov, gcov_se, mean_z1z2, gcov_intercept, gcov_intercept_se, gcorr, gcorr_se, gcorr_zscore, gcorr_pval]) + "\n")

outf.close()

# some trait may never be first trait in pair
# can get num of snps manually or from  munging logs
for trait in traitResultsD.keys():
    if trait not in trait_nSnpD:
        trait_nSnpD[trait] = "NA"


with open(outName2,'w') as outf:
    outf.write("\t".join(["trait", "nSnpMunged", "h2", "h2SE", "lambda_GC", "meanChisq", "intercept", "interceptSE", "ratio"]) + "\n")
    for trait in sorted(traitResultsD.keys()):
        outf.write(trait + "\t" + trait_nSnpD[trait] + "\t" + "\t".join(traitResultsD[trait]) + "\n")



# TODO run on all pairs, record results in a big table for 1 trait, and a big table for all pairs




