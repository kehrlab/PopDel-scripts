#!/bin/python

from sys import argv
from random import randint, seed
from datetime import datetime

if len(argv) != 6:
    print('Usage: ' + argv[0] + ' <seed> <in.txt> <out.txt> <in.fa> <out.fa>')
    exit(1)

inputseed = int(argv[1])
variantfile = argv[2]
selectedvariantfile = argv[3]
inputfile = argv[4]
outputfile = argv[5]


def status(msg):
    print('[{:%Y-%m-%d %H:%M:%S}]'.format(datetime.now())),
    print(msg)


### Load the sequence from the input file.
status('Loading the input sequence...')

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


### Load the simulated variants.
status('Loading the variants...')

variants = []

with open(variantfile, 'r') as varfile:
    for line in varfile:
        if line[0] != '#':
            variants += [line.split()]


### Select a set of variants according to allele frequencies.
status('Selecting a set of variants according to allele frequencies...')

vars = []
seed(inputseed)
for var in variants:
    if randint(0, 1000) < float(var[3]) * 1000:
        vars += [var]


### Write the selected variants to the output file.
status('Writing selected variants to output file...')

vars.sort(key=lambda var: int(var[1]))
with open(selectedvariantfile, 'w') as outfile:
    outfile.write('#' + '\t'.join(['CHROM', 'POS', 'LEN', 'AF']) + '\n')
    for var in vars:
        outfile.write('\t'.join(var) + '\n')


### Insert the variants into the sequence = delete the correponding ranges.
status('Inserting variants into the input sequence...')

vars.sort(key=lambda var: -int(var[1]))
for var in vars:
    seq = seq[:int(var[1])] + seq[int(var[1])+int(var[2]):]


### Output the modified sequence.
status('Writing the modified sequence to output file...')

with open(outputfile, 'w') as outfile:
    outfile.write('>' + id + '|Seed=' + str(inputseed) + '\n')
    for i in range(0, len(seq), 60):
        end = min(i+60, len(seq))
        outfile.write(seq[i:end] + '\n')
