##DEL only
function filter
{
    zcat calls/smoove/${1}sim.smoove.vcf.gz | awk '{if ($1 ~ /^#/ || $5=="<DEL>") print $0}'  > calls/smoove/${1}sim.smoove.vcf
}

for i in {0001..0009} {0010..0100..10} {0200..1000..100}
do
   filter $i
done
