

# maps chr:pos:ref_allele to rsID and adds column to sumstats file
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--maf', help = 'Minor Allele Frequency splitting threshold, default is 0.1. Must be numeric.', default = '0.1')
par
parser.add_argument('--ldsc', required= True, help = 'LD Score splitting threshold, required. Must be numeric.')
parser.add_argument('-s', '--sumstatdir', required=True, help = 'The directory holding ldscore files with MAF and L2. Example ldscore files see LDSC_file_example.txt.')
parser.add_argument('-o', '--outdir', required=True, help = 'The directory of SNP assigned into bins that this code generates
')
parser.add_argument('--suffix', help = 'File containing the list of output suffix for the 4 bins. Default suffix are _maf_0-0.1_ldsc_0-50, _maf_0-0.1_ldsc_50-100, _maf_0.1-0.5_ldsc_0-50 and _maf_0.1-0.5_ldsc_50-100. These correspond to the four bins separated by MAF with a threshold of 0.1, and by LD Score with a threshold of 50% of all LD Scores.')
parser.add_argument('-f', '--filelist', required=True, help = 'File containing the list of ldscore file names in sumstatdir. Each file name is separated by a line. Example ldscore files see LDSC_file_example.txt.')

args = parser.parse_args()

sumstatdir = args.sumstatdir
outdir = args.outdir
listFile = args.filelist

maf_thrd = float(args.maf)
ldsc_thrd = float(args.ldsc)

suffix = args.suffix


# get list of gwas files that need conversion
with open(listFile) as f:
    fileL = [x.rstrip() for x in f.readlines()]

if not suffix:
    out = ["_maf_0-0.1_ldsc_0-50", "_maf_0-0.1_ldsc_50-100", "_maf_0.1-0.5_ldsc_0-50", "_maf_0.1-0.5_ldsc_50-100"]
else:
    with open(suffix) as f:
        out = [x.rstrip() for x in f.readlines()]

# assumes alleles in plinkf are subset of those in bimf
for fname in fileL:
    sumstatName = "%s/%s" % (sumstatdir, fname)
    outNameL = ["%s/%s" % (outdir, fname+out[i]) for i in range(len(out))]
    

    with open(sumstatName) as sumstatf:
        header = sumstatf.readline().rstrip()
        print(header)
        headerL = header.split()
        fieldIndices = {}
        for col in headerL:
            fieldIndices[col] = headerL.index(col)
        print(fieldIndices)

        with open(outNameL[0], 'w') as a, open(outNameL[1], 'w') as b, open(outNameL[2], 'w') as c, open(outNameL[3], 'w') as d:

            outputFileL = [a,b,c,d]

            for outputFile in outputFileL:
                outputFile.write(header + '\n')

            for line in sumstatf:
                lineL = line.rstrip().split()
                outline = line.rstrip() + '\n'
                maf = float(lineL[fieldIndices['MAF']])
                ldsc = float(lineL[fieldIndices['L2']])

                if maf < maf_thrd:
                    if ldsc < ldsc_thrd: 
                # Currently the ldsc and maf thresholds are hardcoded here from the log file generated by ldsc
                        outputFileL[0].write(outline)
                    else:
                        outputFileL[1].write(outline)
                else:
                    if ldsc < ldsc_thrd:
                        outputFileL[2].write(outline)
                    else:
                        outputFileL[3].write(outline)


    print('filtered %s' % sumstatName)


