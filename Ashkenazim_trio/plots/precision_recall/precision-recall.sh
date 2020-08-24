#!/bin/bash
set -euxo pipefail

#### This script gets the TP and FP of the tools PopDel, Delly, Lumppy/Smoove and Manta for the HG002 ground truth in the high-confidence regions and filters them for varying
#### genotype quality thresholds.

truth=scripts/bedtools/truth/HG002_SVs_Tier1_v0.6.vcf.gz
highConf=scripts/bedtools/truth/HG002_SVs_Tier1_v0.6.bed

zcat $truth | grep -v ^# | grep "SVTYPE=DEL" | grep -v "^X" | grep -v "^Y"| awk -F '[\t]' '{f=$10; split(f,gt,":"); if (gt[1] == "0/1" || gt[1] == "1/1") print $0}' | sed 's/;/\t/g' | cut -f 1,2,36 | sed 's;END=;;g' | awk '{if($3 - $2 >= 500 && $3 - $2 <= 10000){print $0}}' | sort -u -n -k1,1 -k2,2 -k3,3 > eval/truth.500-10000.bed
bedtools intersect -wa -u -f 1 -a eval/truth.500-10000.bed -b $highConf > eval/truth_highConf.500-10000.bed
wc -l eval/truth_highConf.500-10000.bed  ##get total truth count


############## POPDEL ##############

tool=popdel
#VCF=calls/${tool}/HG002.${tool}.DelOnly.500-10000.GRCh37.vcf         ## Choose one and adapt OUT name
VCF=calls/${tool}/HG002_joint.${tool}.HQ.DelOnly.500-10000.GRCh37.vcf

grep ^# $VCF > header
bedtools intersect -wa -u -f 1 -a $VCF -b $highConf | cat header - > highConfDel.popdel.vcf

bedtools intersect -wa -u -r -f 0.5 -a highConfDel.popdel.vcf -b eval/truth_highConf.500-10000.bed > TP.popdel.vcf
bedtools intersect -v -r -f 0.5 -a highConfDel.popdel.vcf -b eval/truth_highConf.500-10000.bed > FP.popdel.vcf

OUT=eval/${tool}_PR.tsv
echo -e "QUAL\tFP\tTP" > $OUT
for q in {0..255}
do
   echo -e "$q\t$(awk -v Q="$q" -F '[:;=\t]' '{if ($34>=Q) {print $0}}' FP.popdel.vcf | wc -l)\t$(awk -v Q="$q" -F '[:;=\t]' '{if ($34>=Q) {print $0}}' TP.popdel.vcf | wc -l)" >> $OUT
done


###################### DELLY #############

tool=delly
#VCF=calls/${tool}/HG002_joint.${tool}.HQ.DelOnly.500-10000.GRCh37.vcf ## Choose one and adapt OUT name
VCF=calls/${tool}/HG002.${tool}.HQ.DelOnly.500-10000.GRCh37.vcf

grep ^# $VCF > header
bedtools intersect -wa -u -f 1 -a $VCF -b $highConf | cat header - > highConfDel.vcf

bedtools intersect -wa -u -r -f 0.5 -a highConfDel.vcf -b eval/truth_highConf.500-10000.bed > TP.vcf
bedtools intersect -v -r -f 0.5 -a highConfDel.vcf -b eval/truth_highConf.500-10000.bed > FP.vcf

OUT=eval/${tool}_PR.HG002_single.tsv
echo -e "QUAL\tFP\tTP" > $OUT
for q in {0..255}
do
   echo -e "$q\t$(cut -f10 FP.vcf | awk -v Q="$q" -F '[:;=\t]' '{if ($3>=Q) {print $0}}' | wc -l)\t$(cut -f10 TP.vcf | awk -v Q="$q" -F '[:;=\t]' '{if ($3>=Q) {print $0}}' | wc -l)" >> $OUT
done


######### Lumpy/Smoove ########

tool=smoove
#VCF=calls/${tool}/HG002_joint.${tool}.HQ.DelOnly.500-10000.GRCh37.vcf ## Choose one and adapt OUT name
VCF=calls/${tool}/HG002.${tool}.HQ.DelOnly.500-10000.GRCh37.vcf

grep ^# $VCF > header
bedtools intersect -wa -u -f 1 -a $VCF -b $highConf | cat header - > highConfDel.vcf

bedtools intersect -wa -u -r -f 0.5 -a highConfDel.vcf -b eval/truth_highConf.500-10000.bed > TP.vcf
bedtools intersect -v -r -f 0.5 -a highConfDel.vcf -b eval/truth_highConf.500-10000.bed > FP.vcf

OUT=eval/${tool}_PR.HG002_single.tsv
echo -e "QUAL\tFP\tTP" > $OUT
for q in {0..200}
do
   echo -e "$q\t$(cut -f10 FP.vcf | awk -v Q="$q" -F ':' '{if ($2>=Q) {print $0}}' | wc -l)\t$(cut -f10 TP.vcf | awk -v Q="$q" -F ':' '{if ($2>=Q) {print $0}}' | wc -l)" >> $OUT
done


############# Manta #################
tool=manta
#VCF=calls/${tool}/HG002_joint.${tool}.HQ.DelOnly.500-10000.GRCh37.vcf ## Choose one and adapt OUT name
VCF=calls/${tool}/HG002.${tool}.HQ.DelOnly.500-10000.GRCh37.vcf

grep ^# $VCF > header
bedtools intersect -wa -u -f 1 -a $VCF -b $highConf | cat header - > highConfDel.vcf

bedtools intersect -wa -u -r -f 0.5 -a highConfDel.vcf -b eval/truth_highConf.500-10000.bed > TP.vcf
bedtools intersect -v -r -f 0.5 -a highConfDel.vcf -b eval/truth_highConf.500-10000.bed > FP.vcf

OUT=eval/${tool}_PR.HG002_single.tsv
echo -e "QUAL\tFP\tTP" > $OUT
for q in {0..999}
do
   echo -e "$q\t$(cut -f10 FP.vcf | awk -v Q="$q" -F ':' '{if ($3>=Q) {print $0}}' | wc -l)\t$(cut -f10 TP.vcf | awk -v Q="$q" -F ':' '{if ($3>=Q) {print $0}}' | wc -l)" >> $OUT
done

##### Clean up #####

rm header eval/truth.500-10000.bed eval/truth_highConf.500-10000.bed highConfDel.vcf TP.vcf FP.vcf
