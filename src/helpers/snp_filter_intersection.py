from sets import Set

# maps chr:pos:ref_allele to rsID and adds column to sumstats file

mtagFile = "/u/scratch/c/cmhuang/mtag_snp_list/Main.snps"
filter4File = "/u/scratch/c/cmhuang/mtag_snp_list/filter4_rand200k_indivs.bim_rsID"
outdir = "/u/scratch/c/cmhuang/ukbb_new/"

mtag_set = Set()
with open(mtagFile) as varf:
    header = varf.readline() #skip header
    for line in varf:
        lineL = line.rstrip()
        mtag_set.add(lineL)

filter_set = Set()
filter_dict = {}
with open(filter4File) as varf:
    for line in varf:
        lineL = line.rstrip().split()
        filter_set.add(lineL[6])
        filter_dict[lineL[6]] = lineL[1]


# get list of gwas files that need conversion
inter = filter_set.intersection(mtag_set)


outName = outdir + "filter_intersection"
  
with open(outName, 'w') as outf:
    for snp in inter:
        outf.write(filter_dict[snp] + '\t' + snp + '\n')
  
print('found intersection')


