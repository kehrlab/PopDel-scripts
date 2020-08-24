#!/bin/bash
set -euxo pipefail

function vcf2bed_truth
{
   grep ^${2} ${1} | sed -e 's/;/\t/g' -e 's;SVLEN=;;g' | awk -v OFS='\t' '{{print$1,$2,$2+$8}}' | sort -u -k1,1 -k2,2n -k3,3 > ${1}.bed
}

function vcf2bed_popdel
{
   grep -v "^#" ${1} | sed -e 's/;/\t/g' -e 's;SVLEN=-;;g' | awk -v OFS='\t' '{print$1,$2,$2+$9}' | sort -u -k1,1 -k2,2n -k3,3 > ${1}.bed
}

function vcf2bed_delly
{
   bcftools view ${1} | grep  -v "^#" | grep "SVTYPE=DEL" | sed 's/;/\t/g' | cut -f 1,2,12 | sed 's/END=//g' |  sort -u -k1,1 -k2,2n -k3,3n > ${1}.bed # && $7=="PASS" # Would mainly filter TP on single sample.
}

function vcf2bed_lumpy
{
   zcat ${1}.gz | grep  -v "^#" |  grep "SVTYPE=DEL" | sed 's/;/\t/g' | cut -f 1,2,10 | sed 's/END=//g' |  sort -u -k1,1 -k2,2n -k3,3n > ${1}.bed
}

function vcf2bed_manta
{
   zcat ${1}.gz | grep  -v "^#" |  grep "SVTYPE=DEL"| cut -f 1,2,8 | sed 's/;/\t/g' | cut -f 1-3 | sed 's/END=//g' | sort -u -k1,1 -k2,2n -k3,3n > ${1}.bed # && $7=="PASS" # Would only filter one TP and nothing else (at least in uniform simulation).
}

function vcf2bed
{
   if [ "$1" = "popdel" ]
   then
      vcf2bed_popdel ${2}
   else
      if [ "$1" = "delly" ]
      then
         vcf2bed_delly ${2}
      else
         if [ "$1" = "lumpy" ]
         then
            vcf2bed_lumpy ${2}
         else
            if [ "$1" = "manta" ]
            then
               vcf2bed_manta ${2}
            else
               echo "ERROR! Could not determine tool from file \'${2}\'. Terminating."
               exit 1
            fi
         fi
      fi
   fi
}


function evaluate
{
   tool=$1
   chrom=$2
   chr=${chrom#"chr"}
   callPath=WG_G1k_simulation/calls
   outdir=WG_G1k_simulation/eval/${tool}
   out=${outdir}/${chrom}.${tool}.bed.counts
   if test -f "$out"
   then
     echo "Output file $out already exists. Terminating."
     exit 1
  fi

   mkdir -p ${outdir}
   touch -a ${out}
   for n in {1..10} {20..100..10} {200..500..100}
   do
      if [ "$tool" == "delly" ]
      then
          callVCF=${callPath}/${tool}/${chrom}.${n}.${tool}.bcf
      else
          callVCF=${callPath}/${tool}/${chrom}.${n}.${tool}.vcf
      fi
      calls=${callVCF}.bed
      if [ ! -f $calls ]
      then
         echo "Creating BED-file of $tool calls for $n samples."
         vcf2bed $tool $callVCF
      fi
      zeros="0"
      if [ "$n" -lt "1000" ]
      then
        zeros="00"
      fi
      if [ "$n" -lt "100" ]
      then
        zeros="000"
      fi
      if [ "$n" -lt "10" ]
      then
        zeros="0000"
      fi
      truthVCF=WG_G1k_simulation/variants/${chrom}/deletions.S00001-S${zeros}${n}.vcf
      truth=${truthVCF}.bed
      if [ ! -f $truth ]
      then
         echo "Creating BED-file of truth set for $n samples."
         vcf2bed_truth $truthVCF $chr
      fi
      echo $n > ${out}.tmp
      bedtools intersect -r -f 0.5 -wa -u -a ${truth} -b ${calls} | wc -l | cut -d ' ' -f1 >> ${out}.tmp
      bedtools intersect -r -f 0.5 -wa -v -b ${truth} -a ${calls} | wc -l | cut -d ' ' -f1 >> ${out}.tmp
      bedtools intersect -r -f 0.5 -wa -v -a ${truth} -b ${calls} | wc -l | cut -d ' ' -f1 >> ${out}.tmp
      paste $out ${out}.tmp | sed -e 's;^\t;;g' > ${out}.tmp2
      mv ${out}.tmp2 $out
      rm ${out}.tmp
   done
}

function evaluate_gt
{
  tool=$1
  chrom=$2
  chr=${chrom#"chr"}
  callPath=WG_G1k_simulation/calls
  outdir=WG_G1k_simulation/eval/${tool}/GT
  out=${outdir}/${chrom}.${tool}.py.counts
  outGT=${outdir}/${chrom}.${tool}.GT.counts
  comp=scripts/compare_results.py
  if test -f "$out"
  then
    echo "Output file $out already exists. Terminating."
    exit 1
 fi

  mkdir -p ${outdir}/tmp
  touch -a ${out}

  for n in {1..10} {20..100..10} {200..500..100}
  do
    callVCF=${callPath}/${tool}/${chrom}.${n}.${tool}.vcf
    if [ "$tool" == "delly" ]
    then
      bcftools view ${callPath}/${tool}/${chrom}.${n}.${tool}.bcf > $callVCF
    fi
    if [ "$tool" == "lumpy" ] || [ "$tool" == "manta" ]
    then
      zcat ${callPath}/${tool}/${chrom}.${n}.${tool}.vcf.gz > $callVCF
    fi
    zeros="0"
    if [ "$n" -lt "1000" ]
    then
      zeros="00"
    fi
    if [ "$n" -lt "100" ]
    then
      zeros="000"
    fi
    if [ "$n" -lt "10" ]
    then
      zeros="0000"
    fi
    truthVCF=WG_G1k_simulation/variants/${chrom}/deletions.S00001-S${zeros}${n}.vcf
    echo $n > ${out}.tmp
    python2 $comp $truthVCF $callVCF 300 150 $outdir/tmp > tmp
    head -n3 tmp | sed -e 's;## [A-Z][A-Z]: ;;g' >> ${out}.tmp
    paste $out ${out}.tmp | sed -e 's;^\t;;g' > ${out}.tmp2
    mv ${out}.tmp2 ${out}
    rm ${out}.tmp
    tail -n 1 tmp | sed -e "s;TOTAL;${n};g" >> ${outGT}
    rm tmp
    rm $outdir/tmp.*
    if [ "$tool" == "delly" ] || [ "$tool" == "lumpy" ] || [ "$tool" == "manta" ]
    then
      rm $callVCF
    fi
  done
}

### Creation of truth files:
chromosomes="chr17 chr18 chr19 chr20 chr21 chr22"
tools="popdel delly lumpy manta"
vardir=WG_G1k_simulation/variants
for chrom in $chromosomes
do
  if [ ! -f "${vardir}/${chrom}/deletions.S00001-S00500.vcf" ]
  then
    python2 scripts/simulation2vcf.py ${vardir}/${chrom}/S* > ${vardir}/${chrom}/deletions.S00001-S00500.vcf
  fi
  for i in {001..009} {010..090..10} {100..400..100}
  do
    if [ ! -f "${vardir}/${chrom}/deletions.S00001-S00${i}.vcf" ]
    then
      let j=10#$i+9; less ${vardir}/${chrom}/deletions.S00001-S00500.vcf | cut -f 1-${j} | grep -e "^#" -e "0/1" -e "1/1" > ${vardir}/${chrom}/deletions.S00001-S00${i}.vcf ## The "10#" is not a comment, but forces the interpreation as a base 10 number.
    fi
  done
done

## evaluate the G1k_callsets without genotypes
for tool in $tools
do
  for chrom in $chromosomes
  do
    evaluate $tool $chrom
  done
done

## evaluate the G1k_callsets with genotypes
for tool in $tools
do
  for chrom in $chromosomes
  do
    evaluate_gt $tool $chrom
  done
done
