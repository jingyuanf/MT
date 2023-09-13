# maps chr:pos:ref_allele to rsID and adds column to plink sumstats file

variantsFile = "/u/scratch/c/cmhuang/variants.tsv"
plinkdir = "/u/scratch/c/cmhuang/data_ukbb_other_traits/ukbb_gwas_ss_plink_wSE_linear_150k_MAF"
outdir = "/u/scratch/c/cmhuang/data_ukbb_other_traits/ukbb_gwas_ss_plink_wSE_linear_150k_A2_MAF_rsID"
listFile = "/u/scratch/c/cmhuang/data_ukbb_other_traits/ukbb_gwas_ss_plink_wSE_linear_150k_MAF/plink_assoc_150k_MAF_list.txt"

dict = {} # dict[(chr:pos, A1)] = rsID
with open(variantsFile) as varf:
    for line in varf:
        lineL = line.rstrip().split()
        chro = lineL[1]
        pos = lineL[2]
	chropos = chro + ":" + pos
        A1allele = lineL[4]
	rsID = lineL[5]
        dict[(chropos,A1allele)] = rsID

# get list of gwas files that need conversion
with open(listFile) as f:
    fileL = [x.rstrip() for x in f.readlines()]

# assumes alleles in plinkf are subset of those in bimf
for fname in fileL:
    plinkName = "%s/%s" % (plinkdir, fname)
    outName = "%s/%s_rsID" % (outdir, fname)
    with open(plinkName) as plinkf:
        with open(outName, 'w') as outf:
            header = plinkf.readline() # skip header
            outf.write(header.rstrip() + "\trsID\n")
            for line in plinkf:
                lineL = line.rstrip().split()
                chropos = lineL[1]
                A1 = lineL[3]
		if (chropos,A1) in dict:
                	outline = line.rstrip() + "\t" + dict[(chropos,A1)] + "\n"
                	outf.write(outline)

    print('added rsID to %s' % plinkName)


