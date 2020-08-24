TOOL=manta
SAMPLE=HG002  # one of 'HG002' or 'all'
DIR=calls/${TOOL}
VCF=${DIR}/${SAMPLE}.${TOOL}.GRCh37.vcf.gz
FILTERED=${DIR}/${SAMPLE}.${TOOL}.HQ.DelOnly.GRCh37.vcf
SIZE=${DIR}/${SAMPLE}.${TOOL}.HQ.DelOnly.GRCh37.500-10000.vcf

zcat $VCF | grep -e "SVTYPE=DEL" -e "^#" | grep -e "PASS" -e "^#" | grep -v "^X" | grep -v "^Y" | grep -v "^h" | grep -v "^G" | grep -v "^N" > $FILTERED
awk -F '[\t;=]' '{if ($9-$2 >= 500 && $9-$2 <= 10000 || $1 ~ /^#/ ) print $0}' $FILTERED > $SIZE

## Don't need this if SAMPLE == HG002
HG002Only=${DIR}/HG002_joint.${TOOL}.HQ.DelOnly.GRCh37.vcf
HG002Only_SIZE=${DIR}/HG002_joint.${TOOL}.HQ.DelOnly.500-10000.GRCh37.vcf
grep ^## $FILTERED > t
grep -v ^## $FILTERED | cut -f1-10 | grep -e "^#" -e "1/1" -e "1/0" -e "0/1" > tt
cat t tt > $HG002Only
rm t tt
awk -F '[\t;=]' '{if ($9-$2 >= 500 && $9-$2 <= 10000 || $1 ~ /^#/ ) print $0}' $HG002Only > $HG002Only_SIZE
