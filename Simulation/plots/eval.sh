#!/bin/bash

comp=compare_results.py

function evl
{
   out=TP-FP-FN/${1}/${1}.py.counts
   out2=TP-FP-FN/${1}/${1}.GTcounts
   touch -a ${out}
   for n in 0001 {0010..0100..10} {0200..1000..100}
   do
      truth=truth/deletions.S00001-S0${n}.vcf
      echo ${n} >> ${out}.tmp
      python2 $comp $truth calls/${1}/${n}sim.${1}.vcf 300 150 TP-FP-FN/${1}/${1}.${n} > tmp
      head -n3 tmp | sed -e 's;## [A-Z][A-Z]: ;;g' >> ${out}.tmp
      paste $out ${out}.tmp | sed -e 's;^\t;;g' > ${out}.tmp2
      mv ${out}.tmp2 ${out}
      rm ${out}.tmp
      tail -n 1 tmp | sed -e "s;TOTAL;${n};g" >> ${out2}
      rm tmp
   done
}

function evalNoGT
{
   out=TP-FP-FN/${1}/${1}.py.counts
   touch -a ${out}
   for n in 0001 {0010..0100..10} {0200..0400..100}
   do
      truth=truth/deletions.S00001-S0${n}.vcf
      echo ${n} >> ${out}.tmp
      python2 $comp $truth calls/${1}/${n}sim.${1}.deDup.vcf 300 150 TP-FP-FN/${1}/${1}.${n} > tmp
      head -n3 tmp | sed -e 's;## [A-Z][A-Z]: ;;g' >> ${out}.tmp
      paste $out ${out}.tmp | sed -e 's;^\t;;g' > ${out}.tmp2
      mv ${out}.tmp2 ${out}
      rm ${out}.tmp
      rm tmp
   done
}

evl popdel
evl delly
evl lumpy
evl manta
evalNoGT gridss ##NO GENOTYPES
