#from sets import Set

# maps chr:pos:ref_allele to rsID and adds column to sumstats file

indivFile = "/u/home/f/fujy2038/project-zarlab/BIG/ukbb_new_indiv/w33127_20200820.csv"
filter4File = "/u/home/f/fujy2038/project-zarlab/MT_github/Data/PLINK/SNP_filter/filter4.fam"
outdir = "/u/home/f/fujy2038/project-zarlab/BIG/ukbb_new_indiv/"

indivL = []
with open(indivFile) as indivf:
    for line in indivf:
        lineL = line.rstrip()
        indivL.append(lineL)

filterL = []
with open(filter4File) as varf:
    for line in varf:
        indiv = line.rstrip().split()[0:2]
        if indiv[1] not in indivL:
            filterL.append(tuple(indiv))

outName = outdir + "filter4_withdrawn.txt"
  
with open(outName, 'w') as outf:
    for ind in filterL:
        outf.write('\t'.join(ind) + '\n')
  
print('withdrew participants')


