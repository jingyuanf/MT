# maps chr:pos:ref_allele to rsID and adds column to plink bim file

variantsFile = "/u/scratch/c/cmhuang/variants.tsv"
filterdir = "/u/scratch/c/cmhuang/ukbb_new"
outdir = "/u/scratch/c/cmhuang/mtag_snp_list"

dict = {} # dict[(chr:pos, A1)] = rsID
with open(variantsFile) as varf:
    for line in varf:
        lineL = line.rstrip().split()
        chro = lineL[1]
        pos = lineL[2]
        chropos = chro + ":" + pos
        A1allele = lineL[4]
        rsID = lineL[5]
        dict[chropos, A1allele] = rsID


# assumes alleles in plinkf are subset of those in bimf

fname = "filter4_rand200k_indivs.bim"
filterName = "/u/scratch/c/cmhuang/ukbb_new/withdrawn/filter4_rand200k_indivs.bim"
outName = "%s/%s_rsID" % (outdir, fname)
with open(filterName) as plinkf:
    with open(outName, 'w') as outf:
        for line in plinkf:
            lineL = line.rstrip().split()
            chropos = lineL[1]
            A1 = lineL[4]
            if (chropos,A1) in dict:
                outline = line.rstrip() + "\t" + dict[(chropos,A1)] + "\n"
               	outf.write(outline)

    print('added rsID to %s' % filterName)


