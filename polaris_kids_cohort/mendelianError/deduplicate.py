#!/bin/python
from sys import argv

def overlap(startA, startB, endA, endB): ## Does NOT check for same chromosome!
    lenA = endA - startA
    lenB = endB - startB
    l = max(startA, startB)
    r = min(endA, endB)
    i = r - l
    if i >= 0.5 * lenA and i >= 0.5 * lenB: ## >=50% reciprocal overlap
        return True
    else:
        return False  

def readPredictFile(calls, predictFile):
    with open(predictFile) as pFile:
        for line in pFile:
            if line[0] == '#':
                print(line.rstrip())
            else:
                calls.append(line.rstrip())

def deDup(calls, tool):
    drop = set()
    e = len(calls) - 2
    for i, a in enumerate(calls[:e], start=1):
        chromA = a.split()[0]
        startA = int(a.split()[1])
        if tool == "lumpy":
            endA = int(a.split()[7].split(";", 3)[2].strip("END="))
        elif tool == "delly":
            endA = int(a.split()[7].split(";", 5)[4].strip("END="))
        elif tool == "popdel":
             endA = startA - int(a.split()[7].split(";", 2)[1].strip("SVLEN="))
        for j, b in enumerate(calls[i:], start = i):
            chromB = b.split()[0]
            if chromB != chromA:
                break
            startB = int(b.split()[1])
            if startB > startA + 5000:  ## max del lenght is 10000. After 5000 there can be no reciprocal overlap of at least 50%.
                break
            if tool == "lumpy":
                endB = int(b.split()[7].split(";", 3)[2].strip("END="))
            elif tool == "delly":
                endB = int(b.split()[7].split(";", 5)[4].strip("END="))
            elif tool == "popdel":
                endB = startB - int(b.split()[7].split(";", 2)[1].strip("SVLEN="))
            if overlap(startA, startB, endA, endB):
                drop.add(j)
    for k, c in enumerate(calls):
        if k not in drop:
            print(c)
# ------------------------------------------------------------------------------

if len(argv) != 2:
    print('Usage: ' + argv[0] + ' <predicted.[popdel|delly|lumpy].vcf>')
    exit(1)

vcf = argv[1]
if "lumpy" in vcf:
    tool = "lumpy"
elif "delly" in vcf:
    tool = "delly"
elif "popdel" in vcf:
    tool = "popdel"
else:
    print("Error: Could not determine tool frome vcf-filename. Filename must contain toolname \'popdel\' or \'lumpy\' as a substring. Terminating")
    exit(1)

calls = []
readPredictFile(calls, vcf)
deDup(calls, tool)


            

    











