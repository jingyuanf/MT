import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gwas', required = True, help = 'Folder holding all of the GWAS files')
parser.add_argument('-s', '--suffix', required = True, help = 'The suffix of the GWAS files')
parser.add_argument('-t', '--trait', required = True, help = 'The names of GWAS files (usually trait names)')
parser.add_argument('--binpath', required = True, help = 'The path holding all of the SNP bin files generated by stratify_snps_4bins.py')
parser.add_argument('--binfile', required = True, help = 'File containing the list of names of bin files generated by stratify_snps_4bins.py. Each name is separated in a line.')
parser.add_argument('--outdir', required = True, help = 'The directory of all output')
parser.add_argument('--outfile', required = True, help = 'File containing the suffix names of all 4 outputs')

args = parser.parse_args()

gwaspath = args.gwas
gwasfilesub = args.suffix
traitfile = args.trait
binpath = args.binpath
binfile = args.binfile
outdir = args.outdir
outfile = args.outfile

with open(traitfile) as f:
    traitL = [x.rstrip() for x in f.readlines()]

with open(binfile) as f:
    binsfiles = [x.rstrip() for x in f.readlines()]

with open(outfile) as f:
    outputBinName = [x.rstrip() for x in f.readlines()]

gwasfileL = []
allsnps = []

for i in traitL:
    gwasfileL.append(gwaspath + i + gwasfilesub)

gwasDict = {}

for i in range(len(gwasfileL)):
    traitDict = {}
    traitSnp = set()
    gwasfile = gwasfileL[i]
    with open(gwasfile) as f:
        gwasheaderL = f.readline().rstrip().split()
        for line in f:
            lineL = line.rstrip().split()
            snp = lineL[gwasheaderL.index('SNP')]
            traitDict[snp] = lineL
            traitSnp.add(snp)
        gwasDict[traitL[i]] = traitDict
        allsnps.append(traitSnp)
        


bins = []

for i in binsfiles:
    bins.append(binpath+i)


outputFilesL = []
for i in traitL:
    for j in outputBinName:
        outputFilesL.append(i+j)

for i in range(len(bins)):
    snpids = set()
    binsfile = bins[i]
    with open(binsfile) as f:
        headerL = f.readline().rstrip().split()
        for line in f:
            lineL = line.rstrip().split()
            snpids.add(lineL[headerL.index('SNP')])
    for j in range(len(traitL)):
        snp_intersect = snpids.intersection(allsnps[j])
        with open(outdir+traitL[j]+outputBinName[i], 'w') as g:
            g.write("\t".join(header for header in gwasheaderL) + '\n')
            for k in snp_intersect:
                g.write("\t".join(gwasDict[traitL[j]][k]) + '\n')
            
    
        
        



