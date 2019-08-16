tool=popdel
#tool=delly
#tool=lumpy

OUT=${tool}_Mendel_Q.txt
#OUT=${tool}_Mendel_Q.deDup.txt                                                 ## Without duplicates (apply deduplicate.py on VCF file.)
echo "minQ total consistent" > $OUT
VCF=calls/${tool}/polaris_kids.${tool}.500-10000.NoCentromeres.vcf
#VCF=calls/${tool}/polaris_kids.${tool}.500-10000.NoCentromeres.deDup.vcf       ## Without duplicates (apply deduplicate.py on VCF file.)


echo "minQ consistent relative" > $OUT
for Q in {0..255..5}    ## Lumpy's GQ only goes up to 200
do
   python2 scripts/mendelianError/mendel_extended.py $VCF scripts/mendelianError/kids.ped $Q ${tool} | head -n2 | tail -n1 | cut -f2- | sed -e 's;(;;g' -e 's;%);;g' | paste <(echo $Q) - | sed -e 's;\t; ;g' >> $OUT
done

##For inlcuding the 0-0-0 sites in trios:
#echo "minQ total consistent" > $OUT
#for Q in {0..255..5}
#do
#   python2 mendel.py  $VCF kids.ped $Q ${tool} | head -n2 | sed -e 's; ;;g' | awk -v n="$Q" -F '[:(]' 'NR%2==1{t=$2;next} {print n, t, $2}' >> $OUT
#done




