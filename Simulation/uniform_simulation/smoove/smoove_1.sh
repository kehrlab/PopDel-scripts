#!/bin/bash

bamBase="bam/"
REF="data/reference/GRCh38/genome.fa"
BED="scripts/smoove/exclude.bed"
OUT="results-smoove/"
TIME="/usr/bin/time -ao metrics/time/smoove_call.time -f "
EX='chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr22,chrX,chrY,chrM,chrEBV,~.*_.*'   ##ignore everything except chr21 (negative lookahead does not seem to work for the regex)

for i in {0001..1000}
do
   ${TIME} "${i}\t%e\t%U\t%S\t%M" scripts/smoove/smoove call --outdir ${OUT} --name S0${i} --fasta ${REF} --excludechroms ${EX} --exclude ${BED} --processes 1 --genotype ${bamBase}/S0${i}/S0${i}.bam
done
