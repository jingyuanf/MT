from sets import Set

# given sumstat and filter, removes snps from sumstat not in filter

variantsFile = "/u/scratch/c/cmhuang/BIG/data_mtag/mtag_snp_list/Main.snps"
sumstatdir = "/u/scratch/c/cmhuang/BIG/data_mtag"
outdir = "/u/scratch/c/cmhuang/BIG/data_mtag/filtered"
listFile = "/u/scratch/c/cmhuang/BIG/data_mtag/sumstat_list.txt"

#for snp filter in format:
'''
	header
	rsID
	rsID
	rsID
	...
'''
snp_set = Set()
with open(variantsFile) as varf:
    header = varf.readline() #skip header
    for line in varf:
        lineL = line.rstrip()
        snp_set.add(lineL)

# get list of gwas files that need conversion
with open(listFile) as f:
    fileL = [x.rstrip() for x in f.readlines()]

# assumes alleles in plinkf are subset of those in bimf
for fname in fileL:
    sumstatName = "%s/%s" % (sumstatdir, fname)
    outName = "%s/%s" % (outdir, fname)
    with open(sumstatName) as sumstatf:
        with open(outName, 'w') as outf:
            header = sumstatf.readline().rstrip() # skip header
            headerL = header.split()
            fieldIndices = {}
            for col in headerL:
                fieldIndices[col] = headerL.index(col)

            outf.write(header + '\n')
            for line in sumstatf:
                lineL = line.rstrip().split() 
                snpID = lineL[0]
                chr = lineL[fieldIndices['chr']]
                #print(chr)
                bp = lineL[fieldIndices['bpos']]
		if snpID in snp_set and not (chr == 8 and bp > 7962590 and bp < 11962591):
                    #print snpID
                    outline = line.rstrip() + '\n'
                    outf.write(outline)

    print('filtered %s' % sumstatName)


