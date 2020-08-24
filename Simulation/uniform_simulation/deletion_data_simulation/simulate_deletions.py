#!/usr/bin/python

from sys import argv
from datetime import datetime
from random import seed,randint
from numpy.random import seed as npseed
from numpy.random import uniform # seed,lognormal, exponential
from intervaltree import Interval, IntervalTree

if len(argv) != 5:
    print('Usage: ' + argv[0] + ' <seed> <numvar> <in.fa> <out.txt>')
    exit(1)

inputseed = int(argv[1])
numvar = int(argv[2])
inputfile = argv[3]
outputfile = argv[4]

def status(msg):
    print('[{:%Y-%m-%d %H:%M:%S}]'.format(datetime.now())),
    print(msg)


### Load the sequence from the input file.
status("Loading input file...")

id = None
seq = ''

with open(inputfile, 'r') as infile:
    for line in infile:
        if line[0] == '>' and id == None:
           id = line[1:].strip()
        elif line[0] == '>':
            print('ERROR: Expecting only a single fasta record in input file.')
            exit(1)
        else:
            seq += line.strip()

seqlen = len(seq)


### Collect all intervals of +/- 1000 bp around Ns
status('Collecting intervals around Ns...')

itree = IntervalTree()
first = 0
last = 0
for i in range(0, seqlen):
    if seq[i] not in ['A', 'C', 'G', 'T']:
        if i - 1000 > last:
            if last > 0:
                itree[first:last] = 'N'
            first = i - 1000
        last = i + 1000

        
### Generate a set of deletion variants.
status('Generating deletion variants...')

variants = []
seed(inputseed)
npseed(inputseed)

for i in range(0, numvar):
    # Sample an allele frequency.
    af = uniform(0.0, 1.0) # exponential(.2)
#    while af > 1:
#        af = exponential(.2)

    # Sample a deletion length.
    dlen = int(uniform(100,10000)) # int(lognormal(0, 1) * 1500)
#    while dlen > 10000:
#        dlen = int(lognormal(0, 1) * 1500)

    # Sample a position
    pos = randint(1000, seqlen - 1000)
    while len(itree[pos:pos+dlen]) != 0:
        pos = randint(1000, seqlen - 1000)

    # Add the variant to the list of variants and to the interval tree of forbidden positions.
    variants += [[id, str(pos), str(dlen), str(af)]]
    itree[pos - 1000:pos + dlen + 1000] = 'D'


### Write the variants to the output file.
status('Writing variants to output file...')

variants.sort(key=lambda var: int(var[1]))

with open(outputfile, 'w') as outfile:
    outfile.write('#' + '\t'.join(['CHROM', 'POS', 'LEN', 'AF']) + '\n')
    for var in variants:
        outfile.write('\t'.join(var) + '\n')
