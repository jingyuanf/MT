# splits variant column (chr:bpos:A2:A1)

nealedir = "/u/scratch/c/cmhuang/UKB_Neale_round1_sep2017"
outdir = "/u/scratch/c/cmhuang/UKB_Neale_round1_sep2017/new_reformatted"
listFile = "/u/scratch/c/cmhuang/UKB_Neale_round1_sep2017/neale_list.txt"

# get list of gwas files that need conversion
with open(listFile) as f:
    fileL = [x.rstrip() for x in f.readlines()]

# assumes alleles in plinkf are subset of those in bimf
for fname in fileL:
    nealeName = "%s/%s" % (nealedir, fname)
    outName = "%s/reformatted_%s" % (outdir, fname)
    with open(nealeName) as nealef:
        with open(outName, 'w') as outf:
            header = nealef.readline().rstrip() # skip header
            #add chr, bp, a1, a2, and MAF columns to neale group
            outf.write(header + "\tchr\tbpos\ta2\ta1\tMAF\n")
            for line in nealef:
                lineL = line.rstrip().split()
                #split variant column to get desired columns
                chroposa1a2 = lineL[0].split(':')
                outline = line.rstrip()
                for item in chroposa1a2:
                    outline += "\t" + item
                #MAF = allele count / (2 * number of nonmissing individuals)
                outline += '\t' + str(float(lineL[3])/(2 * float(lineL[2]))) + '\n'
                outf.write(outline)

    print('reformatted %s' % nealeName)


