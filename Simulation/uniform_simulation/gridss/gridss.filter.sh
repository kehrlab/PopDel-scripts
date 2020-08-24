function filter
{
   Rscript --vanilla scripts/annotate.R calls/gridss/sim${1}.gridss.vcf > calls/gridss/sim${1}.gridss.annotated.vcf
}

for i in {0010..0100..10} {0200..0500..100}
do
   filter $i
done
