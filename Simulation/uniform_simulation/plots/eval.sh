#!/bin/bash
set -euxo pipefail

comp=/home/sebastian/data/uniform_simulation/scripts/compare_results.py

function evaluate
{
   tmpdir=/home/sebastian/data/uniform_simulation/TP-FP-FN/${1}
   out=/home/sebastian/data/uniform_simulation/eval/${1}/${1}.py.counts
   out2=/home/sebastian/data/uniform_simulation/eval/${1}/${1}.GT.counts
   if [ -e $out ]
   then
     echo "Output file $out already exists. Terminating."
     exit 1
   fi
   if [ -e $out2 ]
   then
     echo "Output file $out2 already exists. Terminating."
     exit 1
   fi

   [ -e ${tmpdir}/col.tmp ] && rm ${tmpdir}/col.tmp
   [ -e ${tmpdir}/tpfpfn.tmp ] && rm ${tmpdir}/tpfpfn.tmp
   [ -e ${tmpdir}/gt.tmp ] && rm ${tmpdir}/gt.tmp
   mkdir -p $tmpdir
   touch -a ${tmpdir}/tpfpfn.tmp

   for n in {0001..0010} {0020..0100..10} {0200..1000..100}
   do
      truth=/home/sebastian/data/uniform_simulation/truth/deletions.S00001-S0${n}.vcf
      if [ "$1" == "popdel" ]
      then
        calls=/home/sebastian/data/uniform_simulation/calls/${1}/S00001-S0${n}.${1}.vcf
      fi
      if [ "$1" == "delly_n" ]
      then
        calls=${tmpdir}/S00001-S0${n}.${1}.vcf
        bcftools view /home/sebastian/data/uniform_simulation/calls/${1}/0${n}sim.${1}.bcf > $calls
      fi
      if [ "$1" == "smoove" ]
      then
        calls=${tmpdir}/S00001-S0${n}.${1}.vcf
        zcat /home/sebastian/data/uniform_simulation/calls/${1}/${n}sim.${1}.vcf.gz > $calls
      fi
      if [ "$1" == "manta" ]
      then
        calls=/home/sebastian/data/uniform_simulation/calls/${1}/${n}sim.${1}.vcf
      fi
      python2 $comp $truth $calls 300 150 ${tmpdir}/${1}.${n} > ${tmpdir}/tmp
      echo ${n} > ${tmpdir}/col.tmp
      head -n3 ${tmpdir}/tmp | sed -e 's;## [A-Z][A-Z]: ;;g' >> ${tmpdir}/col.tmp
      paste ${tmpdir}/tpfpfn.tmp ${tmpdir}/col.tmp | sed -e 's;^\t;;g' > ${tmpdir}/o4.tmp && mv ${tmpdir}/o4.tmp ${tmpdir}/tpfpfn.tmp
      rm ${tmpdir}/col.tmp
      tail -n 1 ${tmpdir}/tmp | sed -e "s;TOTAL;${n};g" >> ${tmpdir}/gt.tmp
      rm ${tmpdir}/tmp
   done
   mv ${tmpdir}/tpfpfn.tmp ${out}
   mv ${tmpdir}/gt.tmp ${out2}
   rm ${tmpdir}/*
   wait
}

function evalNoGT
{
   out=/vol/local/data/uniformSimulation/TP-FP-FN/${1}/${1}.py.counts
   touch -a ${out}
   for n in {0001..010} {0020..0100..10} {0200..0400..100}
   do
      truth=/vol/local/data/uniformSimulation/truth/deletions.S00001-S0${n}.vcf
      echo ${n} >> ${out}.tmp
      python2 $comp $truth /home/sebastian/data/uniform_simulation/calls/${1}/${n}sim.${1}.deDup.vcf 300 150 ${tmpdir}/${1}.${n}  > tmp
      head -n3 tmp | sed -e 's;## [A-Z][A-Z]: ;;g' >> ${out}.tmp
      paste $out ${out}.tmp | sed -e 's;^\t;;g' > ${out}.tmp2
      mv ${out}.tmp2 ${out}
      rm ${out}.tmp
      rm tmp
   done
}

evaluate popdel
evaluate delly_n
evaluate smoove
evaluate manta
evalNoGT gridss ##NO GENOTYPES
