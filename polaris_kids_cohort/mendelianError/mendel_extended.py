#!/bin/python

from sys import argv

def readPredictFile(gtMatrix, predictFile, minQ, tool):
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
                if tool == "popdel" or tool == "delly":
		    gq = int(s.split(":", 3)[2])
                elif tool == "lumpy":
                    gq = s.split(":", 3)[1]
                    if gq != ".":
                        gq = int(gq)
                    else:
                        gq = 0
                s = s.split(":", 1)[0]        ## Read genotypes of samples
                if s == "1/1" and gq >= minQ:                ## Replace genotypes with numbers
                    gts.append(2)
                elif s == "0/1" and gq >= minQ:
                    gts.append(1)
                elif s == "0/0" and gq >= minQ:
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

if len(argv) != 5:
    print 'Usage: ' + argv[0] + ' <predicted.vcf> <ped.txt> <minQ> <popdel|delly|lumpy>'
    exit(1)
predictFile = argv[1]
pedFile = argv[2]
minQ = int(argv[3])
tool = argv[4]
if tool not in ["popdel","delly","lumpy"]:
    print ("Invalid tool \'" + tool + "\'. Must be one of \'popdel\', \'delly\', \'lumpy\'.")
    exit(1)
# Store a dictonnary of vcf records corresponding to the predictions.
# Keys are chrom:pos:PD and if non-unique, _<i> is appended
gtMatrix = []
sampleMap = readPredictFile(gtMatrix, predictFile, minQ, tool)
trios = loadPedFile(pedFile)

## Counts
correct = 0
## Consistent
ReReRe = 0      ## 0/0 + 0/0 -> 0/0
ReHeRe = 0      ## 0/0 + 0/1 -> 0/0
HeHeRe = 0      ## 0/1 + 0/1 -> 0/0
ReHeHe = 0      ## 0/0 + 0/1 -> 0/1
ReHoHe = 0      ## 0/0 + 1/1 -> 0/1
HeHeHe = 0      ## 0/1 + 0/1 -> 0/1
HeHoHe = 0      ## 0/1 + 1/1 -> 0/1
HeHeHo = 0      ## 0/1 + 0/1 -> 1/1
HeHoHo = 0      ## 0/1 + 1/1 -> 1/1
HoHoHo = 0      ## 1/1 + 1/1 -> 1/1
## Inconsistent
incorrect = 0
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
                    ReReRe += 1
                elif f1 == 1:
                    incorrect += 1
                    ReReHe += 1
                elif f1 == 2:
                    incorrect += 1
                    ReReHo += 1
            elif p2 == 1:
                if f1 == 0:
                     correct += 1
                     ReHeRe += 1
                elif f1 == 1:
                     correct += 1
                     ReHeHe += 1
                elif f1 == 2:
                    incorrect += 1
                    ReHeHo += 1
            elif p2 == 2:
                if f1 == 0:
                    incorrect += 1
                    ReHoRe += 1
                elif f1 == 1:
                    correct += 1
                    ReHoHe += 1
                elif f1 == 2:
                    incorrect += 1
                    ReHoHo += 1
        elif p1 == 1:
            if p2 == 0:
                if f1 == 0:
                    correct += 1
                    ReHeRe += 1
                elif f1 == 1:
                    correct += 1
                    ReHeHe += 1
                elif f1 == 2:
                    incorrect += 1
                    ReHeHo +=1
            elif p2 == 1:
                if f1 == 0:
                    correct += 1
                    HeHeRe += 1
                elif f1 == 1:
                    correct += 1
                    HeHeHe += 1
                elif f1 == 2:
                    correct += 1
                    HeHeHo += 1
            elif p2 == 2:
                if f1 == 0:
                    incorrect += 1
                    HeHoRe += 1
                elif f1 == 1:
                    correct += 1
                    HeHoHe += 1
                elif f1 == 2:
                    correct += 1
                    HeHoHo += 1
        elif p1 == 2:
            if p2 == 0:
                if f1 == 0:
                    incorrect += 1
                    ReHoRe += 1
                elif f1 == 1:
                    correct += 1
                    ReHoHe += 1
                elif f1 == 2:
                    incorrect += 1
                    ReHoHo += 1
            elif p2 == 1:
                if f1 == 0:
                    incorrect += 1
                    HeHoRe += 1
                elif f1 == 1:
                    correct += 1
                    HeHoHe += 1
                elif f1 == 2:
                    correct += 1
                    HeHoHo += 1
            elif p2 == 2:
                if f1 == 0:
                    incorrect += 1
                    HoHoRe += 1
                elif f1 == 1:
                    incorrect += 1
                    HoHoHe += 1
                elif f1 == 2:
                    correct += 1
                    HoHoHo += 1

##print the results
tot = float(correct + incorrect)
#print 'Total:\t' + str(int(tot))
if tot == 0:
    print ('Consistent:\t0')
    print ('Errors    :\t0')
else:
    print 'Consistent:\t' + str(correct) + '\t(' + "%.4f" % (correct / tot * 100) + '%)'
    print 'Consistent_w/o_all_Ref:\t' + str(correct - ReReRe) + '\t(' + "%.4f" % ((correct - ReReRe) / (tot - ReReRe) * 100) + '%)'
    print 'Re+Re->Re:\t' + str(ReReRe)  + '\t(' + "%.4f" % (ReReRe  / tot * 100) + '%)'
    print 'Re+He->Re:\t' + str(ReHeRe)  + '\t(' + "%.4f" % (ReHeRe  / tot * 100) + '%)'
    print 'He+He->Re:\t' + str(HeHeRe)  + '\t(' + "%.4f" % (HeHeRe  / tot * 100) + '%)'
    print 'Re+He->He:\t' + str(ReHeHe)  + '\t(' + "%.4f" % (ReHeHe  / tot * 100) + '%)'
    print 'Re+Ho->He:\t' + str(ReHoHe)  + '\t(' + "%.4f" % (ReHoHe  / tot * 100) + '%)'
    print 'He+He->He:\t' + str(HeHeHe)  + '\t(' + "%.4f" % (HeHeHe  / tot * 100) + '%)'
    print 'He+He->Ho:\t' + str(HeHeHo)  + '\t(' + "%.4f" % (HeHeHo  / tot * 100) + '%)'
    print 'He+Ho->He:\t' + str(HeHoHe)  + '\t(' + "%.4f" % (HeHoHe  / tot * 100) + '%)'
    print 'He+Ho->Ho:\t' + str(HeHoHo)  + '\t(' + "%.4f" % (HeHoHo  / tot * 100) + '%)'
    print 'Ho+Ho->Ho:\t' + str(HoHoHo)  + '\t(' + "%.4f" % (HoHoHo  / tot * 100) + '%)'

    print 'Errors:\t' + str(int(incorrect))+ '\t(' + "%.4f" % (incorrect / tot * 100) + '%)'
    print 'Re+Re->He:\t' + str(ReReHe)  + '\t(' + "%.4f" % (ReReHe  / tot * 100) + '%)'
    print 'Re+Re->Ho:\t' + str(ReReHo)  + '\t(' + "%.4f" % (ReReHo  / tot * 100) + '%)'
    print 'Re+He->Ho:\t' + str(ReHeHo)  + '\t(' + "%.4f" % (ReHeHo  / tot * 100) + '%)'
    print 'Re+Ho->Re:\t' + str(ReHoRe)  + '\t(' + "%.4f" % (ReHoRe  / tot * 100) + '%)'
    print 'Re+Ho->Ho:\t' + str(ReHoHo)  + '\t(' + "%.4f" % (ReHoHo  / tot * 100) + '%)'
    print 'He+Ho->Re:\t' + str(HeHoRe)  + '\t(' + "%.4f" % (HeHoRe  / tot * 100) + '%)'
    print 'Ho+Ho->Re:\t' + str(HoHoRe)  + '\t(' + "%.4f" % (HoHoRe  / tot * 100) + '%)'
    print 'Ho+Ho->He:\t' + str(HoHoHe)  + '\t(' + "%.4f" % (HoHoHe  / tot * 100) + '%)'




##Meaning of error abbreviations:
## Re - Reference (0/0)
## He - Heterozygous (0/1) or (1/0)
## Ho - Homozygous variant (1/1)
## The first two genotypes belong to the parents, the third one to the offspring.
