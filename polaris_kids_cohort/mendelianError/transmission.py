#!/usr/bin/python

import sys
from sys import argv

if len(argv) != 4:
    print("Usage: " + argv[0] + "<PED file> <VCF file> <popdel|delly|lumpy>")
    exit(1)

pedfilename = argv[1]
vcffilename = argv[2]
tool = argv[3]

# ------------------------------------------------------------------------------

def getGenotype(field, tool):
    fields = field.split(':')

    gt  = -1
    if fields[0] == '0/0':
        gt = 0
    elif fields[0] == '0/1':
        gt = 1
    elif fields[0] == '1/1':
        gt = 2

    gq = -1
    if tool == "popdel":
        pls = [int(x) for x in fields[1].split(',')]
        pls.sort()
        gq = pls[1]
    elif tool == "lumpy":
        if gt != -1:
            gq = int(fields[1])
    elif tool == "delly":
        gq = min(int(fields[2]), 255)
    else:
        print("Unknown tool: " + tool)
        exit(1)

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
allCounts = [[0 for i in range(18) ] for j in range(256)]
singleCounts = [[0 for i in range(18) ] for j in range(256)]

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

        # Read genotypes of the variant
        allele_count = 0
        for i in range(9, len(header)):
            gt, gq = getGenotype(line[i], tool)
            genotypes[header[i]] = (gt, gq)

        p1, p2, c, gq = [None]*4
        singleTrio = True

        for trio in trios:

            parent1 = genotypes[trio[0]]
            parent2 = genotypes[trio[1]]
            child = genotypes[trio[2]]

            if not isGenotyped(parent1, parent2, child):
                continue

            minQ = min(min(parent1[1], parent2[1]), child[1])
            addToMatrix(allCounts, parent1[0], parent2[0], child[0], minQ)

            # Check if this is the only trio that carries this variant
            if (parent1[0] + parent2[0] + child[0]) != 0:
                if p1 != None:
                    singleTrio = False
                else:
                    p1 = parent1[0]
                    p2 = parent2[0]
                    c = child[0]
                    gq = minQ

        if singleTrio and p1 != None:
            addToMatrix(singleCounts, p1, p2, c, gq)


#OUPUT
print('\t'.join(['AnyTrio', '0-0-0', '0-0-1', '0-0-2', '0-1-0', '0-1-1', '0-1-2', '0-2-0', '0-2-1', '0-2-2', '1-1-0', '1-1-1', '1-1-2', '1-2-0', '1-2-1', '1-2-2', '2-2-0', '2-2-1', '2-2-2']))
minQ = 0
for row in allCounts:
    print('\t'.join(str (x) for x in [minQ] + row))
    minQ += 1

print('\t'.join(['OneTrio', '0-0-0', '0-0-1', '0-0-2', '0-1-0', '0-1-1', '0-1-2', '0-2-0', '0-2-1', '0-2-2', '1-1-0', '1-1-1', '1-1-2', '1-2-0', '1-2-1', '1-2-2', '2-2-0', '2-2-1', '2-2-2']))
minQ = 0
for row in singleCounts:
    print('\t'.join(str (x) for x in [minQ] + row))
    minQ += 1
