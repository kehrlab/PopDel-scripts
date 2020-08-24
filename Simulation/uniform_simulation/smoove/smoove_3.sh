#!/bin/bash

bamBase="bam"
REF="data/reference/GRCh38/genome.fa"
OUT="results-smoove/"
TIME="/usr/bin/time -ao metrics/time/smoove_genotype.time -f "


for batch in {0002..0009} {0010..0100..10} {0200..1000..100}
do
   for sample in $( eval echo {0001..$batch} )
   do
      ${TIME} "${batch}\t${sample}\t%e\t%U\t%S\t%M" scripts/smoove/smoove genotype -d -x -p 1 --name S0${sample}joint_batch0${batch} --outdir ${OUT}genotyped/  --fasta ${REF} --vcf ${OUT}merged_S0${batch}.sites.vcf.gz ${bamBase}/S0${sample}/S0${sample}.bam
   done
done
