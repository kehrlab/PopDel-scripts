#!/bin/bash
VCF="calls/popdel/polaris150.rnd.popdel.500-10000.NoCentromeres.vcf"
OUT="calculations/pca/popdel.NoCentromeres.perSampleVariants.GTfilter.txt"

## Replace '26' with your desired minimum GQ

## Count all variants
for i in {10..159..1} ## Adapt upper limit to the number of samples (9 + samplenum)
do
grep -v '^##' $VCF | head -n1 | cut -f $i > tmp
grep -v '^#' $VCF | cut -f $i | awk -F ':' '{if ( ($1 == "0/1" && $3 >= 26) || ($1 == "1/1" && $3 >= 26) ) count++} END {print count}' | paste tmp - >> $OUT
done
rm tmp

OUT="calculations/pca/popdel.NoCentromeres.GTfilter.hetHom.txt"
## Count HomVar and HetVar
for i in {10..159..1} ## Adapt number of samples (9 + samplenum)
do
grep -v '^##' $VCF | head -n1 | cut -f $i > tmp
grep -v '^#' $VCF | cut -f $i | awk -F ':' '{OFS="\t"} {if ($1 == "1/1" && $3 >= 26) hom++; else if ($1 == "0/1" && $3 >= 26) het++} END {print het, hom}' | paste tmp - >> $OUT
done
rm tmp
## Count HetVar
