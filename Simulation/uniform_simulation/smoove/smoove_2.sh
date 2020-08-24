#!/bin/bash

REF="data/reference/GRCh38/genome.fa"
OUT="results-smoove/"
TIME="/usr/bin/time -ao metrics/time/smoove_merge.time -f "
for batch in {0002..0009} {0010..0100..10} {0200..1000..100}
do
   vcfs=""
   for sample in $( eval echo {0001..$batch} )
   do
      vcfs="${vcfs} results-smoove/S0${sample}-smoove.genotyped.vcf.gz"
   done
   ${TIME} "${batch}\t%e\t%U\t%S\t%M" scripts/smoove/smoove merge --outdir ${OUT} --name merged_S0${batch} --fasta ${REF} $vcfs
done
