#!/bin/python

from sys import argv

def readPredictFile(gtMatrix, predictFile):
    samples = []
    with open(predictFile) as pFile:
        for line in pFile:
            if line[0] == '#':
                if line[1] != '#' and len(samples) == 0:
                    sampleMap = {}
                    i = 0
                    for s in line.split()[9:]:
                        sampleMap[s] = i      ##Create a map sample:index
                        i += 1
                continue
            genotypes = line.split()[9:]
            gts = []
            for s in genotypes:
                s = s.split(":", 1)[0]        ## Read genotypes of samples
                if s == "1/1":                ## Replace genotypes with numbers
                    gts.append(2)
                elif s == "0/1":# or s == "1/0":
                    gts.append(1)
                elif s == "0/0":
                    gts.append(0)
                else:
                    gts.append(-1)          ## not genotyped
            gtMatrix.append(gts)
    return sampleMap

def loadPedFile(pedFile):
    trios = []
    with open(pedFile) as pFile:
        for line in pFile:
            trio = line.split()
            trios.append(trio)
    return trios

# ------------------------------------------------------------------------------

if len(argv) != 3:
    print 'Usage: ' + argv[0] + ' <predicted.vcf> <ped.txt>'
    exit(1)
predictFile = argv[1]
pedFile = argv[2]
# Store a dictonnary of vcf records corresponding to the predictions.
# Keys are chrom:pos:PD and if non-unique, _<i> is appended
gtMatrix = []
sampleMap = readPredictFile(gtMatrix, predictFile)
trios = loadPedFile(pedFile)

##Counts of error types##
correct = 0
ReReHe = 0    ## 0/0 + 0/0 -> 0/1
ReReHo = 0    ## 0/0 + 0/0 -> 1/1
ReHeHo = 0    ## 0/0 + 0/1 -> 1/1
ReHoRe = 0    ## 0/0 + 1/1 -> 0/0
ReHoHo = 0    ## 0/0 + 1/1 -> 1/1
HeHoRe = 0    ## 0/1 + 1/1 -> 0/0
HoHoRe = 0    ## 1/1 + 1/1 -> 0/0
HoHoHe = 0    ## 1/1 + 1/1 -> 0/1

for var in gtMatrix:                        ## for each variant...
    for trio in trios:                      ## check each trio
        if trio[0] in sampleMap:
            p1index = sampleMap[trio[0]]    ## get the index of each trio member
        else:
            print ('Warning: Sample ' + trio[0] + ' is not part of VCF file. Skipping trio ') + "-".join(trio) + '.'
            trios.remove(trio)
            continue
        if trio[1] in sampleMap:
            p2index = sampleMap[trio[1]]
        else:
            print 'Warning: Sample ' + trio[1] + ' is not part of VCF file. Skipping trio ' + "-".join(trio) + '.'
            trios.remove(trio)
            continue
        if trio[2] in sampleMap:
            f1index = sampleMap[trio[2]]
        else:
            print 'Warning: Sample ' + trio[2] + ' is not part of VCF file. Skipping trio ' + "-".join(trio) + '.'
            trios.remove(trio)
            continue
        p1 = var[p1index]               ## get the genotype predictions
        p2 = var[p2index]
        f1 = var[f1index]
        if p1 == 0:
            if p2 == 0:
                if f1 == 0:
                    correct += 1
                elif f1 == 1:
                    ReReHe += 1
                elif f1 == 2:
                    ReReHo += 1
            elif p2 == 1:
                if f1 == 0 or f1 == 1:
                     correct += 1
                elif f1 == 2:
                    ReHeHo += 1
            elif p2 == 2:
                if f1 == 0:
                    ReHoRe += 1
                elif f1 == 1:
                    correct += 1
                elif f1 == 2:
                    ReHoHo += 1
        elif p1 == 1:
            if p2 == 0:
                if f1 == 0 or f1 == 1:
                    correct += 1
                elif f1 == 2:
                    ReHeHo +=1
            elif p2 == 1 and f1 != -1:
                correct += 1
            elif p2 == 2:
                if f1 == 0:
                    HeHoRe += 1
                elif f1 == 1 or f1 == 2:
                    correct += 1
        elif p1 == 2:
            if p2 == 0:
                if f1 == 0:
                    ReHoRe += 1
                elif f1 == 1:
                    correct += 1
                elif f1 == 2:
                    ReHoHo += 1
            elif p2 == 1:
                if f1 == 0:
                    HeHoRe += 1
                elif f1 == 1 or f1 == 2:
                    correct += 1
            elif p2 == 2:
                if f1 == 0:
                    HoHoRe += 1
                elif f1 == 1:
                    HoHoHe += 1
                elif f1 == 2:
                    correct += 1

##print the results
tot = float(correct + ReReHe + ReReHo + ReHeHo + ReHoRe + ReHoHo + HeHoRe + HoHoRe + HoHoHe)
err = float(tot - correct)
print 'Total     : ' + str(int(tot))
print 'Consistent: ' + str(correct) + '\t(' + "%.4f" % (correct / tot * 100) + '%)'
print 'Errors    : ' + str(int(err))+ '\t(' + "%.4f" % (err     / tot * 100) + '%)'
print 'Re+Re->He : ' + str(ReReHe)  + '\t(' + "%.4f" % (ReReHe  / tot * 100) + '%)'
print 'Re+Re->Ho : ' + str(ReReHo)  + '\t(' + "%.4f" % (ReReHo  / tot * 100) + '%)'
print 'Re+He->Ho : ' + str(ReHeHo)  + '\t(' + "%.4f" % (ReHeHo  / tot * 100) + '%)'
print 'Re+Ho->Re : ' + str(ReHoRe)  + '\t(' + "%.4f" % (ReHoRe  / tot * 100) + '%)'
print 'Re+Ho->Ho : ' + str(ReHoHo)  + '\t(' + "%.4f" % (ReHoHo  / tot * 100) + '%)'
print 'He+Ho->Re : ' + str(HeHoRe)  + '\t(' + "%.4f" % (HeHoRe  / tot * 100) + '%)'
print 'Ho+Ho->Re : ' + str(HoHoRe)  + '\t(' + "%.4f" % (HoHoRe  / tot * 100) + '%)'
print 'Ho+Ho->He : ' + str(HoHoHe)  + '\t(' + "%.4f" % (HoHoHe  / tot * 100) + '%)'

##Meaning of error abreviations:
## Re - Reference (0/0)
## He - Heterozygous (0/1) or (1/0)
## Ho - Homozygous variant (1/1)
## The first two genotypes belong to the parents, the third one to the offspring.
