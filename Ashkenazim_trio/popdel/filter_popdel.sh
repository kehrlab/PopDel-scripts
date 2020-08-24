TOOL=popdel
SAMPLE=all #HG002  # one of 'HG002' or 'all'
DIR=calls/${TOOL}
VCF=${DIR}/${SAMPLE}.${TOOL}.GRCh37.vcf
SIZE=${DIR}/${SAMPLE}.${TOOL}.GRCh37.500-10000.vcf

awk -F '[\t;=]' '{if ($10 <= -500 && $10 >= -10000 || $1 ~ /^#/ ) print $0}' $VCF > $SIZE

## Don't need this if SAMPLE == HG002
HG002Only=${DIR}/HG002_joint.${TOOL}.HQ.DelOnly.GRCh37.vcf
HG002Only_SIZE=${DIR}/HG002_joint.${TOOL}.HQ.DelOnly.500-10000.GRCh37.vcf
grep ^## $VCF > t
grep -v ^## $VCF | cut -f1-10 | grep -e "^#" -e "1/1" -e "1/0" -e "0/1" > tt
cat t tt > $HG002Only
rm t tt
awk -F '[\t;=]' '{if ($10 <= -500 && $10 >= -10000 || $1 ~ /^#/ ) print $0}' $HG002Only > $HG002Only_SIZE
