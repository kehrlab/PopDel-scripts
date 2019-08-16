#!/bin/python

''' Assumes both input VCF files to be sorted first by chromosome and second by position.

Legend for counts:
-------------------------------
         S 0/0   S 0/1   S 1/1     S .. Simulated genotype
-------------------------------
 P 0/0     -      FN1     FN2      P .. Predicted genotype
 P 0/1    FP1     TP1     FN3
 P 1/1    FP2     FP3     TP2
-------------------------------

'''

from sys import argv
import re
from intervaltree import Interval, IntervalTree

class VcfRecord:
    def __init__(self, fields):
        self.chrom = fields[0]
        self.pos = int(fields[1])
        self.end = self.pos
        self.id = fields[2]
        self.alt = fields[4] ##only necessary for GRIDSS
        self.found = []

        qual = fields[5]
        if qual == '.':
            self.qual = None
        else:
            self.qual = round(float(qual))

        if ':' in self.alt: ##Only necessary for GRIDSS
            if len(re.split(':|\[|\]', self.alt)) >= 3: ## Manta can have this smaller than 2
                self.end = int(re.split(':|\[|\]', self.alt)[2])

        info = fields[7].split(';')
        self.info = {}
        for i in info:
            i = i.split('=')
            if len(i) == 1:
                i += [True]
            self.info[i[0]] = i[1]

        gts = fields[9:]
        self.gts = []
        for gt in gts:
            gt = gt.split(':')
            self.gts += [gt[0]]

    def str(self):
        infoFields = []
        for i in self.info:
            if self.info[i] == True:
                infoFields += [i]
            else:
                infoFields += [i + '=' + self.info[i]]
        info = ';'.join(infoFields)
        return '\t'.join([self.chrom, str(self.pos), self.id, 'N', '<DEL>', str(self.qual), '.', info, 'GT'] + self.gts)

def readSimFile(itrees, records, simFile):
    samples = []
    with open(simFile, 'r') as sFile:
        for line in sFile:
            if line[0] == '#':
                if line[1] != '#' and len(samples) == 0:
                    samples = line.split()[9:]
                continue

            rec = VcfRecord(line.split())

            if rec.chrom not in itrees:
                itrees[rec.chrom] = IntervalTree()

            itrees[rec.chrom][rec.pos - pDelta:rec.pos + pDelta] = rec.id
            records[rec.id] = rec
    return samples

def openVcfWrite(filename, headerfile):
    outfile = open(filename, 'w')
    with open(headerfile, 'r') as infile:
        for line in infile:
            if line[0] == '#':
                outfile.write(line)
    return outfile

def readPredictFile(records, predictFile):
    with open(predictFile) as pFile:
        for line in pFile:
            if line[0] == '#':
                continue

            rec = VcfRecord(line.split())
            key = rec.chrom + ':' + str(rec.pos) + ':PD'

            k = key
            i = 1
            while key in records:
                key = k + '_' + str(i)
                i += 1

            rec.id = key
            records[key] = rec

def addFalsePositive(counts, gts):
    for i in range(0, len(gts)):
        if gts[i] == '0/1':
            counts[i]['FP1'] += 1
        elif gts[i] == '1/1':
            counts[i]['FP2'] += 1

def addFalseNegative(counts, gts):
    for i in range(0, len(gts)):
        if gts[i] == '0/1':
            counts[i]['FN1'] += 1
        elif gts[i] == '1/1':
            counts[i]['FN2'] += 1

def addTruePositive(counts, sGts, records, ids):
    # Select the matching prediction with the highest likelihood ratio.
    maxLR = 0
    if 'LR' in records[ids[0]].info:
        which = 0
        for id in ids:
            if records[id].info['LR'] > maxLR:
                pGts = records[id].gts
                maxLR = records[id].info['LR']
                which = id
    else:
        pGts = records[ids[0]].gts
    rec = records[which]
    for i in range(0, len(sGts)):
        if sGts[i] == '0/0':
            if pGts[i] == '0/1':
                counts[i]['FP1'] += 1
            if pGts[i] == '1/1':
                counts[i]['FP2'] += 1
        elif sGts[i] == '0/1':
            if pGts[i] == '0/1':
                counts[i]['TP1'] += 1
            if pGts[i] == '1/1':
                counts[i]['FP3'] += 1
        elif sGts[i] == '1/1':
            if pGts[i] == '0/1':
                counts[i]['FN3'] += 1
            if pGts[i] == '1/1':
                counts[i]['TP2'] += 1

def printPerSampleCounts(counts, samples):
    total = {'TP1':0, 'TP2':0, 'FP1':0, 'FP2':0, 'FP3':0, 'FN1':0, 'FN2':0, 'FN3':0}
    print '\t'.join(['#SAMPLE'] + total.keys())
    print '##'

    for i in range(0, len(samples)):
        print '\t'.join([samples[i]] + [str(x) for x in counts[i].values()])

    for key in total:
        for c in counts:
            total[key] += c[key]
    print '##'
    print '\t'.join(['TOTAL'] + [str(x) for x in total.values()])


# ------------------------------------------------------------------------------

if len(argv) != 6:
    print 'Usage: ' + argv[0] + ' <simulated.vcf> <predicted.vcf> <pos delta> <length delta> <prefix for output files>'
    exit(1)

simFile = argv[1]
predictFile = argv[2]
pDelta = int(argv[3])
lDelta = int(argv[4])
prefix = argv[5]

# Store one interval tree per chromosome.
# Each interval tree stores one interval of size 2*pDelta per simulated variant around the variants true position.
itrees = {}
sRecords = {}
samples = readSimFile(itrees, sRecords, simFile)

# Store a dictonnary of vcf records corresponding to the predictions.
# Keys are chrom:pos:PD and if non-unique, _<i> is appended
pRecords = {}
readPredictFile(pRecords, predictFile)

# Open output files for true positives (TP), false positives (FP) and false negatives (FN).
tpSimFile = openVcfWrite(prefix + '.TP.sim.vcf', simFile)
tpPredFile = openVcfWrite(prefix + '.TP.pred.vcf', predictFile)
fpFile = openVcfWrite(prefix + '.FP.vcf', predictFile)
fnFile = openVcfWrite(prefix + '.FN.vcf', simFile)

FN = 0
FP = 0
TP = 0
counts = []
for i in range(0, len(samples)):
    counts += [{'TP1':0, 'TP2':0, 'FP1':0, 'FP2':0, 'FP3':0, 'FN1':0, 'FN2':0, 'FN3':0}]

# Assign predicted records to simulated variants and count false positives.
for r in pRecords:
    rec = pRecords[r]

    if not rec.chrom in itrees: ## Check for FP on wrong chromosome.
        sRecs = {}
    else:
        sRecs = itrees[rec.chrom][rec.pos]
    if len(sRecs) == 0:
        FP += 1
        addFalsePositive(counts, rec.gts)
        fpFile.write(rec.str() + '\n')
    else:
        for id in sRecs:
            sLen = int(sRecords[id[2]].info['SVLEN'])

            plen = None
            if 'SVLEN' in rec.info: ## PopDel & Lumpy
                pLen = -int(rec.info['SVLEN'])
            elif 'END' in rec.info: ## Delly
                pLen = int(rec.info['END']) - rec.pos
            else:                   ## GRIDSS
                pLen = abs(rec.end - rec.pos)
            if pLen < sLen + lDelta and pLen > sLen - lDelta:
                sRecords[id[2]].found += [r]
                tpPredFile.write(str(id) + '\t' + rec.str() + '\n')
            else:
                FP += 1
                addFalsePositive(counts, rec.gts)
                fpFile.write(rec.str() + '\n')

# Count true positives and false negatives.
for r in sRecords:
    rec = sRecords[r]
    if len(rec.found) == 0:
        FN += 1
        addFalseNegative(counts, rec.gts)
        fnFile.write(rec.str() + '\n')
    else:
        TP += 1 # If len(rec.found) is larger than 1, this causes TP + FP < Total predicted
        addTruePositive(counts, rec.gts, pRecords, rec.found)
        tpSimFile.write(rec.str() + '\n')

tpSimFile.close()
tpPredFile.close()
fpFile.close()
fnFile.close()

print '## TP: ' + str(TP)
print '## FP: ' + str(FP)
print '## FN: ' + str(FN)
print '##'
printPerSampleCounts(counts, samples)
