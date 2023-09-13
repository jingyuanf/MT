from sets import Set

# maps chr:pos:ref_allele to rsID and adds column to sumstats file

bimFile = "/u/scratch/c/cmhuang/ukbb_new/gwas_ss_plink_wSE_A2_filtered/armfat_percent.rand200k_indivs.assoc.linear_filtered"
outdir = "/u/scratch/c/cmhuang/mtag_snp_list/"

bim_set = Set()
with open(bimFile) as varf:
    header = varf.readline() #skip header
    for line in varf:
        lineL = line.rstrip().split()
        bim_set.add(lineL[1])

outName = outdir + "filter_snps"
  
with open(outName, 'w') as outf:
    for snp in bim_set:
        outf.write(snp + '\n')
  
print('found intersection')


