#!/bin/bash
VCF="calls/popdel/polaris150.rnd.popdel.500-10000.NoCentromeres.vcf"
OUT="popdel.perSampleVariants.noCent.GTfilter.txt"

## Replace 27 with your desired minimum GQ

## Count all variants
for i in {10..159..1} ## Adapt upper limit to the number of samples (9 + samplenum)
do
grep -v '^##' $VCF | head -n1 | cut -f $i > tmp
grep -v '^#' $VCF | cut -f $i | awk -F '[:,]' '{if ( ($1 == "0/1" && $4-$3 >= 27 && $2-$3 >= 27) || ($1 == "1/1" && $2-$4 >= 30 && $3-$4 >= 27) ) count++} END {print count}' | paste tmp - >> $OUT
done
rm tmp

OUT="popdel.perSampleVariants.noCent.GTfilter.hetHom.txt"
## Count HomVar and HetVar
for i in {10..159..1} ## Adapt number of samples (9 + samplenum)
do
grep -v '^##' $VCF | head -n1 | cut -f $i > tmp
grep -v '^#' $VCF | cut -f $i | awk -F '[:,]' '{if ($1 == "1/1" && $2-$4 >= 30 && $3-$4 >= 27) hom++; else if ($1 == "0/1" && $4-$3 >= 27 && $2-$3 >= 27) het++} END {print het, hom}' | paste tmp - >> $OUT
done
rm tmp
## Count HetVar
