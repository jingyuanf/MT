from sets import Set

# given plink sumstat and filter including chr:bp, removes SNPs not in filter

variantsFile = "/u/scratch/c/cmhuang/ukbb_new/gwas_ss_plink_wSE_A2_filtered/armfat_percent.rand200k_indivs.assoc.linear_filtered"
sumstatdir = "/u/scratch/c/cmhuang/ukbb_new/gwas_ss_plink_wSE_A2"
outdir = "/u/scratch/c/cmhuang/ukbb_new/gwas_ss_plink_wSE_A2_filtered"
#sumstatdir = "/u/scratch/c/cmhuang/data_ukbb_other_traits/ukbb_gwas_ss_plink_wSE_linear"
#outdir = "/u/scratch/c/cmhuang/data_ukbb_other_traits/150k_filtered"
listFile = "/u/scratch/c/cmhuang/ukbb_new/gwas_ss_plink_wSE/gwas_list.txt"
#listFile = "/u/scratch/c/cmhuang/ukbb_new/gwas_name_list.txt"

snp_set = Set()
with open(variantsFile) as varf:
    header = varf.readline() #skip header
    for line in varf:
        lineL = line.rstrip().split()
        snp = lineL[1] #+ ':' + lineL[3]
        snp_set.add(snp)

# get list of gwas files that need conversion
with open(listFile) as f:
    fileL = [x.rstrip() for x in f.readlines()]

for fname in fileL:
    sumstatName = "%s/%s_A2" % (sumstatdir, fname)
    #sumstatName = "%s/%s.rand150k_indivs.assoc.linear" % (sumstatdir, fname)
    outName = "%s/%s.rand150k_indivs.assoc.linear_filtered" % (outdir, fname)
    with open(sumstatName) as sumstatf:
        with open(outName, 'w') as outf:
            header = sumstatf.readline().rstrip()
            outf.write(header + '\n')
            for line in sumstatf:
                lineL = line.rstrip().split() 
                snpID = lineL[1] #+ ":" + lineL[3]
                if snpID in snp_set:    
                    #print snpID
                    outline = line.rstrip() + '\n'
                    outf.write(outline)

    print('filtered %s' % sumstatName)


