#!/usr/bin/python

import sys
from sys import argv

if len(argv) != 3:
    print("Usage: " + argv[0] + " <PED file> <VCF file>")
    exit(1)

pedfilename = argv[1]
vcffilename = argv[2]

# ------------------------------------------------------------------------------

def getIndexInFormat(formatString, what):
    return formatString.split(":").index(what)

def getGenotype(field, gqIdx):
    fields = field.split(':')

    gt = -1
    gq = -1
    if fields[0] == '0/0':
        gt = 0
    elif fields[0] == '0/1' or fields[0] == '1/0':
        gt = 1
    elif fields[0] == '1/1':
        gt = 2
    gq = fields[gqIdx]
    if (gq != "."):
        gq = min(int(gq), 255)

    return gt, gq

def isGenotyped(gt1, gt2, gt3):
    if gt1[0] == -1 or gt2[0] == -1 or gt3[0] == -1:
        return False
    return True

def isHighQ(gt1, gt2, gt3, minQ):
    if gt1[1] < minQ or gt2[1] < minQ or gt3[1] < minQ:
        return False
    return True

def addToColumn(allCounts, minQ, column):
    for i in range(0, minQ+1):
        allCounts[i][column] += 1

def addToMatrix(allCounts, p1, p2, c, minQ):
    if p1 == 0 and p2 == 0 and c == 0:
        addToColumn(allCounts, minQ, 0)
    elif p1 == 0 and p2 == 0 and c == 1:
        addToColumn(allCounts, minQ, 1)
    elif p1 == 0 and p2 == 0 and c == 2:
        addToColumn(allCounts, minQ, 2)
    elif p1 == 0 and p2 == 1 and c == 0:
        addToColumn(allCounts, minQ, 3)
    elif p1 == 0 and p2 == 1 and c == 1:
        addToColumn(allCounts, minQ, 4)
    elif p1 == 0 and p2 == 1 and c == 2:
        addToColumn(allCounts, minQ, 5)
    elif p1 == 0 and p2 == 2 and c == 0:
        addToColumn(allCounts, minQ, 6)
    elif p1 == 0 and p2 == 2 and c == 1:
        addToColumn(allCounts, minQ, 7)
    elif p1 == 0 and p2 == 2 and c == 2:
        addToColumn(allCounts, minQ, 8)
    elif p1 == 1 and p2 == 0 and c == 0:
        addToColumn(allCounts, minQ, 3)
    elif p1 == 1 and p2 == 0 and c == 1:
        addToColumn(allCounts, minQ, 4)
    elif p1 == 1 and p2 == 0 and c == 2:
        addToColumn(allCounts, minQ, 5)
    elif p1 == 1 and p2 == 1 and c == 0:
        addToColumn(allCounts, minQ, 9)
    elif p1 == 1 and p2 == 1 and c == 1:
        addToColumn(allCounts, minQ, 10)
    elif p1 == 1 and p2 == 1 and c == 2:
        addToColumn(allCounts, minQ, 11)
    elif p1 == 1 and p2 == 2 and c == 0:
        addToColumn(allCounts, minQ, 12)
    elif p1 == 1 and p2 == 2 and c == 1:
        addToColumn(allCounts, minQ, 13)
    elif p1 == 1 and p2 == 2 and c == 2:
        addToColumn(allCounts, minQ, 14)
    elif p1 == 2 and p2 == 0 and c == 0:
        addToColumn(allCounts, minQ, 6)
    elif p1 == 2 and p2 == 0 and c == 1:
        addToColumn(allCounts, minQ, 7)
    elif p1 == 2 and p2 == 0 and c == 2:
        addToColumn(allCounts, minQ, 8)
    elif p1 == 2 and p2 == 1 and c == 0:
        addToColumn(allCounts, minQ, 12)
    elif p1 == 2 and p2 == 1 and c == 1:
        addToColumn(allCounts, minQ, 13)
    elif p1 == 2 and p2 == 1 and c == 2:
        addToColumn(allCounts, minQ, 14)
    elif p1 == 2 and p2 == 2 and c == 0:
        addToColumn(allCounts, minQ, 15)
    elif p1 == 2 and p2 == 2 and c == 1:
        addToColumn(allCounts, minQ, 16)
    elif p1 == 2 and p2 == 2 and c == 2:
        addToColumn(allCounts, minQ, 17)

# ------------------------------------------------------------------------------

trios = []
with open(pedfilename, 'r') as pedfile:
    for line in pedfile:

        line = line.split()
        if len(line) != 3:
            print("Unexpected line in pedigree file: " + '\t'.join(line))
            exit(1)

        trios += [line]

header = []
genotypes = {}
allCounts = [[[0 for i in range(18) ] for j in range(256)] for k in range(50)]
anyCounts = [[0 for i in range(18) ] for j in range(256)]

with open(vcffilename, 'r') as vcffile:
    for line in vcffile:

        # Skip VCF header
        if line[0] == '#' and line[1] == '#':
            continue

        line = line.split()

        # Get sample names
        if line[0][0] == '#':
            header = line
            header[0] = header[0][1:]
            continue

        # Only consider variants on autosomes
        if line[0] not in ['chr' + str(x) for x in range(1,23)]:
            continue

        # Get Index of the genotype quality
        gqIdx = getIndexInFormat(line[8], "GQ")

        # Read genotypes of the variant
        allele_count = 0
        for i in range(9, len(header)):
            gt, gq = getGenotype(line[i], gqIdx)
            genotypes[header[i]] = (gt, gq)

        # Count trios carrying the variant
        numTrios = 0
        for trio in trios:

            parent1 = genotypes[trio[0]]
            parent2 = genotypes[trio[1]]
            child = genotypes[trio[2]]

            if not isGenotyped(parent1, parent2, child):
                continue

            if (parent1[0] + parent2[0] + child[0]) != 0:
                numTrios += 1

        # Add counts to matrix
        for trio in trios:

            parent1 = genotypes[trio[0]]
            parent2 = genotypes[trio[1]]
            child = genotypes[trio[2]]

            if not isGenotyped(parent1, parent2, child):
                continue

            minQ = min(min(parent1[1], parent2[1]), child[1])
            ## Uncomment for printing de novo candidates:
            ##if child[0] == 1 and parent1[0] + parent2[0] == 0 and minQ > 50:
            ##    print(line[0] + ":" + line[1] + " is a candidate de novo in trio " + trio[0] + "," + trio[1] + "," + trio[2] + " with min(GQ)=" + str(minQ) + ". " + str(parent1[1]) + "  " + str(parent2[1]) + "  " + str(child[1]) )
            ##############
            addToMatrix(anyCounts, parent1[0], parent2[0], child[0], minQ)
            addToMatrix(allCounts[numTrios], parent1[0], parent2[0], child[0], minQ)

#OUPUT
print('\t'.join(['NTrio', 'GQ', '0-0-0', '0-0-1', '0-0-2', '0-1-0', '0-1-1', '0-1-2', '0-2-0', '0-2-1', '0-2-2', '1-1-0', '1-1-1', '1-1-2', '1-2-0', '1-2-1', '1-2-2', '2-2-0', '2-2-1', '2-2-2', 'TransmissionRate']))
minQ = 0
for row in anyCounts:
    counts = [minQ] + row
    notTransmitted = counts[4] + counts[14]
    transmitted = counts[5] + counts[15]
    if transmitted + notTransmitted == 0:
        transmission = float("nan")
    else:
        transmission = float(transmitted) / (transmitted + notTransmitted) * 100
    print('\t'.join(str (x) for x in ['NA', minQ] + row) + '\t' + str(transmission))
    minQ += 1

for i in range(1,len(allCounts)):
    minQ = 0
    for row in allCounts[i]:
        counts = [minQ] + row
        notTransmitted = counts[4] + counts[14]
        transmitted = counts[5] + counts[15]
        if transmitted + notTransmitted == 0:
            transmission = float("nan")
        else:
            transmission = float(transmitted) / (transmitted + notTransmitted) * 100
        print('\t'.join(str (x) for x in [i, minQ] + row) + '\t' + str(transmission))
        minQ += 1
