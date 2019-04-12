function filter
{
   grep -v '#' calls/gridss/${1}sim.gridss_orig.vcf | grep -v "LOW_QUAL" | awk -F '[_\t]' '$4 ~ /h/{next} 1' > calls/gridss/${1}sim.gridss.vcf
}

for i in {0010..0100..10} {0200..0500..100}
do
   filter $i
done
