#!/bin/python

from sys import argv
from datetime import datetime

if len(argv) < 2:
    print 'Usage: ' + argv[0] + ' <sample dir 1> [... <sample dir N>]'
    exit(1)

samplepaths = argv[1:]

def status(msg):
    print '[{:%Y-%m-%d %H:%M:%S}]'.format(datetime.now()),
    print msg
    

def loadHaplotype(deletions, hap1):
    with open(hap1, 'r') as infile:
        l = 0
        for line in infile:
            l += 1
            if line[0] == '#':
                continue
            line = line.strip()
            
            if line in deletions:
                deletions[line] = '1/1'
            else:
                deletions[line] = '0/1'

def addSample(vcf, sample, deletions):
    for d in deletions:
        if d in vcf:
            vcf[d][sample] = deletions[d]
        else:
            vcf[d] = {sample:deletions[d]}

def printVcfHeader(samples):
    print '##fileformat=VCFv4.3'
    print '##fileDate=[{:%Y%m%d}]'.format(datetime.now())
    print '##source=simulation2vcf.py'
    print '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">'
    print '##INFO=<ID=AFS,Number=A,Type=Float,Description="Allele frequency input to simulation.">'
    print '##INFO=<ID=SVLEN,Number=1,Type=String,Description="Difference in length between REF and ALT alleles">'
    print '##Format=<ID=GT,Number=1,Type=String,Description="Genotype">'
    print '\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + samples)
            
# ------------------------------------------------------------------------------
            
# Load two haplotype files per sample.
vcf = {}
samples = []
for samplepath in samplepaths:

    # Get the haplotype file names.
    if samplepath[-1] == '/':
        samplepath = samplepath[:-1]
    sample = samplepath.split('/')[-1]
    samples += [sample]
    prefix = samplepath + '/' + sample
    hap1 = prefix + '-hap1.deletions.txt'
    hap2 = prefix + '-hap2.deletions.txt'
    
    # Load the files.
    deletions = {}
    #status('Loading ' + hap1 + ' and ' + hap2 + ' ...')
    loadHaplotype(deletions, hap1)
    loadHaplotype(deletions, hap2)
    
    # Add the sample to the final vcf records.
    addSample(vcf, sample, deletions)

# Print the vcf file.
records = []
for v in vcf:
    record = v.split()
    record[3] = record[3].split('.')
    record[3] = record[3][0] + '.' + record[3][1][0:4]
    count = 0
    for s in samples:
        if s in vcf[v].keys():
            record += [vcf[v][s]]
            if vcf[v][s] == '0/1':
                count += 1
            if vcf[v][s] == '1/1':
                count += 2
        else:
            record += ['0/0']
    record = record[:4] + [str(count/2.0/len(samples))] + record[4:]
    records += [record]

records.sort(key=lambda r: int(r[1])) 

printVcfHeader(samples)
for line in records:
    print '\t'.join([line[0], line[1], line[0] + ':' + line[1], 'N', '<DEL>', '.', '.', 'SVLEN=' + line[2] + ';AF=' + line[4] + ';AFS=' + line[3], 'GT'] + line[5:])