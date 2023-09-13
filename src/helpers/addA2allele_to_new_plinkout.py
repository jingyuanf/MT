# adds A2 allele column to plink output (both ref and alt allele required by LDSC, LDpred)
# for future reference, plink may have an option to keep both alleles (e.g. if A ref G alt should be treated differently as G ref A alt, there's some option to handle that)

#bim file from which to take A2 alleles
bimName = "PATH_TO_BIM_FILE"
#directory containing files to convert
plinkdir = "PATH_TO_ORIGINAL_PLINK_RESULT(FOLDER)"
#directory to output converted files to 
outdir = "PATH_TO_OUTPUT_DIRECTORY(FOLDER)"
#list of files to convert, example see pheno_name_list.txt
listFile = "TEXT_FILE_CONTAINING_A_LIST_OF_PHENOTYPE_NAMES_AS_PREFIX_OF_FILENAMES"

bimD = {} # bimD[(chro, pos)] = A2 allele
with open(bimName) as bimf:
    for line in bimf:
        lineL = line.rstrip().split()
        #hardcoded for plinks
        chro = lineL[0]
        pos = lineL[3]
        A2allele = lineL[-1]
        bimD[(chro,pos)] = A2allele

# get list of gwas files that need conversion
with open(listFile) as f:
    fileL = [x.rstrip() for x in f.readlines()]

# assumes alleles in plinkf are subset of those in bimf
for fname in fileL:
    plinkName = "%s/%s.rand200k_indivs.assoc.linear" % (plinkdir, fname)
    outName = "%s/%s.rand200k_indivs.assoc.linear_A2" % (outdir, fname)
    with open(plinkName) as plinkf:
        with open(outName, 'w') as outf:
            header = plinkf.readline() # skip header
            outf.write(header.rstrip() + "\tA2\n")
            for line in plinkf:
                lineL = line.rstrip().split()
                chro = lineL[0]
                pos = lineL[2]
                outline = line.rstrip() + "\t" + bimD[(chro,pos)] + "\n"
                outf.write(outline)

    print('added A2 to %s' % plinkName)


