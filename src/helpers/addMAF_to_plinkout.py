# adds MAF column to plink output (needed for mtag)

#.frq file are SNP frequency files generated from plink in step 1.4
frqName = "DIRECTORY_OF_FRQ_FILE"
plinkdir = "DIRECTORY_OF_PLINK_OUTPUT_FILE_WITHOUT_MAF_ADDED"
outdir = "OUTPUT_DIRECTORY"
listFile = "TEXT_FILE_INCLUDING_A_LIST_OF_PLINK_OUTPUT_FILE_NAMES_WITHOUT_MAF"

frqD = {} # frq[(chro, pos)] = MAF
with open(frqName) as frqf:
    for line in frqf:
        lineL = line.rstrip().split()
        chro = lineL[0]
        SNP = lineL[1]
        frq = lineL[4]
        frqD[(chro,SNP)] = frq

# get list of gwas files that need conversion
with open(listFile) as f:
    fileL = [x.rstrip() for x in f.readlines()]

for fname in fileL:
    plinkName = "%s/%s" % (plinkdir, fname)
    outName = "%s/%s_MAF" % (outdir, fname)
    with open(plinkName) as plinkf:
        with open(outName, 'w') as outf:
            header = plinkf.readline() # skip header
            outf.write(header.rstrip() + "\tMAF\n")
            for line in plinkf:
                lineL = line.rstrip().split()
                chro = lineL[0]
                SNP = lineL[1]
                outline = line.rstrip() + "\t" + frqD[(chro,SNP)] + "\n"
                outf.write(outline)

    print('added frq to %s' % plinkName)


