#!/bin/bash


function vcf2bed_truth
{
   grep ^"chr21" ${1} | sed -e 's/;/\t/g' -e 's;SVLEN=;;g' | awk -v OFS='\t' '{{print$1,$2,$2+$8}}' | sort -u -k1,1 -k2,2n -k3,3 > ${1}.bed
}

function vcf2bed_popdel
{
   grep -v "^#" ${1} | sed -e 's/;/\t/g' -e 's;SVLEN=-;;g' | awk -v OFS='\t' '{print$1,$2,$2+$9}' | sort -u -k1,1 -k2,2n -k3,3 > ${1}.bed
}

function vcf2bed_delly
{
   grep  -v "^#" ${1} | awk '{if ($5=="<DEL>") print $0}' | sed 's/;/\t/g' | cut -f 1,2,12 | sed 's/END=//g' |  sort -u -k1,1 -k2,2n -k3,3n > ${1}.bed # && $7=="PASS" # Would mainly filter TP on single sample.
}

function vcf2bed_lumpy
{
   grep  -v "^#" ${1} | awk '{if ($5=="<DEL>") print $0}' | sed 's/;/\t/g' | cut -f 1,2,10 | sed 's/END=//g' |  sort -u -k1,1 -k2,2n -k3,3n > ${1}.bed
}

function vcf2bed_manta
{
   grep  -v "^#" ${1} | awk '{if ($5=="<DEL>") print $0}' | cut -f 1,2,8 | sed 's/;/\t/g' | cut -f 1-3 | sed 's/END=//g' | sort -u -k1,1 -k2,2n -k3,3n > ${1}.bed # && $7=="PASS" # Would only filter one TP and nothing else.
}

function vcf2bed_gridss
{
   grep  -v "^#" ${1} | grep -v "LOW_QUAL" | awk -F '[_\t]' '$4 ~ /h/{next} 1' | cut -f 1-5 | sed -e 's/[ACGT]*//g' -e 's/:/\t/g' -e 's/\[//g' -e 's/\]//g' | cut -f 1,2,6 |  sort -u -k1,1 -k2,2n -k3,3n > ${1}.bed
}


function vcf2bed
{
   if [ "$1" = "popdel" ]
   then
      vcf2bed_popdel ${2}
   else
      if [ "$1" = "delly" ] || [ "$1" = "delly_n" ]
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
               if [ "$1" = "gridss" ]
               then
                  vcf2bed_gridss ${2}
               else
                  echo "ERROR! Could not determine tool!"
               fi
            fi
         fi
      fi
   fi
}

function eval
{
   outdir=${1}
   out=${outdir}/${1}.bed.counts
   mkdir -p ${outdir}
   touch -a ${out}
   for n in 0001 0010 0020 0030 0040 0050 0060 0070 0080 0090 0100 0200 0300 0400 0500 # 0600 0700 0800 0900 1000
   do
      callVCF=/vol/local/data/uniformSimulation/calls/${1}/${n}sim.${1}.vcf
      calls=${callVCF}.bed
      if [ ! -f ${calls} ]
      then
         echo "Creating BED-file of ${1} calls for ${n} samples."
         vcf2bed ${1} ${callVCF}
      fi
      truthVCF=truth/deletions.S00001-S0${n}.vcf
      truth=${truthVCF}.bed
      if [ ! -f ${truth} ]
      then
         echo "Creating BED-file of truth set for ${n} samples."
         vcf2bed_truth ${truthVCF}
      fi
      echo ${n} > ${out}.tmp
      bedtools intersect -r -f 0.5 -wa -u -a ${truth} -b ${calls} | wc -l | cut -d ' ' -f1 >> ${out}.tmp
      bedtools intersect -r -f 0.5 -wa -v -b ${truth} -a ${calls} | wc -l | cut -d ' ' -f1 >> ${out}.tmp
      bedtools intersect -r -f 0.5 -wa -v -a ${truth} -b ${calls} | wc -l | cut -d ' ' -f1 >> ${out}.tmp
      paste $out ${out}.tmp | sed -e 's;^\t;;g' > ${out}.tmp2
      mv ${out}.tmp2 ${out}
      rm ${out}.tmp
   done
}


eval popdel
eval delly
eval delly_n
eval lumpy
eval manta
eval gridss
