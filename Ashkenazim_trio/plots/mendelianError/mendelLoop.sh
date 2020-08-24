tool=popdel
#tool=delly
#tool=smoove

OUT=${tool}_Mendel_Q.txt
VCF=calls/${tool}/all.${tool}.HQ.DelOnly.GRCh37.500-10000.vcf

echo "minQ consistent relative" > $OUT
for Q in {0..255..1}    ## Lumpy/Smoove's GQ only goes up to 200
do
   python2 scripts/mendelianError/mendel_extended.py $VCF scripts/mendelianError/Ashkenazim.ped $Q ${tool} | head -n2 | tail -n1 | cut -f2- | sed -e 's;(;;g' -e 's;%);;g' | paste <(echo $Q) - | sed -e 's;\t; ;g' >> $OUT
done
