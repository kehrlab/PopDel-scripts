#!/usr/bin/python

from sys import argv

if len(argv) != 4:
    print("Usage: " + argv[0] + " <Min GT likelihood> <PED file> <VCF file>")
    exit(1)

minPL = int(argv[1])
pedfilename = argv[2]
vcffilename = argv[3]

# ------------------------------------------------------------------------------

def getGenotype(field, minPL):
    fields = field.split(':')

    gt  = -1
    if fields[0] == '0/0':
        gt = 0
    elif fields[0] == '0/1':
        gt = 1
    elif fields[0] == '1/1':
        gt = 2

    pls = [int(x) for x in fields[1].split(',')]
    pls.sort()

    if pls[1] < minPL:
        return gt, -1

    return gt, gt

def isGenotyped(gt1, gt2, gt3):
    if gt1 == -1 or gt2 == -1 or gt3 == -1:
        return False
    return True

def isConsistent(p1, p2, c):
    if c == 0 and (p1 == 2 or p2 == 2):
        return False
    if c == 2 and (p1 == 0 or p2 == 0):
        return False
    if c == 1 and p1 == 0 and p2 == 0:
        return False
    if c == 1 and p1 == 2 and p2 == 2:
        return False
    return True

def isInformative(p1, p2, c):
    if p1 + p2 == 1 and c != 2:
        return True
    if p1 + p2 == 3 and c != 0:
        return True
    return False

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
with open(vcffilename, 'r') as vcffile:
    for line in vcffile:

        if line[0] == '#' and line[1] == '#':
            continue

        line = line.split()

        if line[0][0] == '#':
            header = line
            header[0] = header[0][1:]
            continue

        if line[0] == 'chrX':
            continue

        af = float(line[7].split(";")[3][3:])

        allele_count = 0
        for i in range(9, len(header)):
            alleles, gt = getGenotype(line[i], minPL)
            genotypes[header[i]] = gt
            if alleles != -1:
                allele_count += alleles

        tr_avg = "NA"
        tr_std = "NA"
        num_informative_trios = 0
        num_consistent_trios = 0
        num_gt_trios = 0

        trs = []
        for trio in trios:
            parent1 = genotypes[trio[0]]
            parent2 = genotypes[trio[1]]
            child = genotypes[trio[2]]

            if not isGenotyped(parent1, parent2, child):
                continue
            num_gt_trios += 1

            #if allele_count == 1 and child == 1:
            #    print(line[0] + ":" + line[1] + " is a candidate de novo in trio " + trio[0] + "," + trio[1] + "," + trio[2] + ".")

            if not isConsistent(parent1, parent2, child):
                continue
            num_consistent_trios += 1

            if not isInformative(parent1, parent2, child):
                continue
            num_informative_trios += 1

            gtSum = parent1 + parent2 + child
            if gtSum == 2 or gtSum == 5:
                trs += [1.0]
            elif gtSum == 1 or gtSum ==4:
                trs += [0.0]
            else:
                print("Genotypes are not consistent or " + line[0] + ":" + line[1] + " in trio " + trio[0] + "," + trio[1] + "," + trio[2] + ".")
                exit(1)

        if len(trs) != num_informative_trios:
            print("Inconsistent number of trios for " + line[0] + ":" + line[1])
            exit(1)

        if num_informative_trios > 0:
            tr_avg = sum(trs) / len(trs)
            tr_std = (sum((x - tr_avg)**2 for x in trs) / len(trs))**0.5

        print('\t'.join([str(x) for x in [line[0], line[1], tr_avg, tr_std, num_informative_trios, num_consistent_trios, num_gt_trios, str(af)]]))
