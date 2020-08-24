#!/bin/bash

OUT="results-smoove/"
TIME="/usr/bin/time -ao metrics/time/smoove_paste.time -f "

for batch in {0002..0009} {0010..0100..10} {0200..1000..100}
do
   ${TIME} "${batch}\t%e\t%U\t%S\t%M" scripts/smoove/smoove paste -o ${OUT}pasted --name cohort_0${batch} ${OUT}genotyped/S*joint_batch0${batch}-smoove.genotyped.vcf.gz
done
