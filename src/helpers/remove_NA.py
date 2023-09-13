# removes individuals with missing phenotype value (required for ldpred)

phenodir = "/u/scratch/c/cmhuang/ukbb_new/pheno/filter4"
outdir = "/u/scratch/c/cmhuang/ukbb_new/pheno_noNA"
listFile = "/u/scratch/c/cmhuang/ukbb_new/pheno_list_linear.txt"

# get list of phenos files that need conversion
with open(listFile) as f:
    fileL = [x.rstrip() for x in f.readlines()]

# assumes alleles in plinkf are subset of those in bimf
for fname in fileL:
    phenoName = "%s/%s" % (phenodir, fname)
    outName = "%s/%s" % (outdir, fname)
    with open(phenoName) as phenof:
        with open(outName, 'w') as outf:
            for line in phenof:
                lineL = line.rstrip().split()
                pheno = lineL[2]
                if pheno != "NA":
                    outline = line
                    outf.write(outline)

    print('removed NA from %s' % phenoName)


