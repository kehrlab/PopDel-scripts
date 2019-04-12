for n in {0010..0100..10} 200
do
   echo -n "BAM=\"" >> t
   for i in $( eval echo {00001..$n} )
   do
      echo -n "bam/S${i}/S${i}.bam," >> t
   done
   sed -i 's;,$;\";' t
   echo -e "\nREF=\"reference/GRCh38/genome.fa\"" >> t
   cat t >> lumpy/${n}.svtyper.sh
   rm t
   echo -e "\nsvtyper -B \$BAM -T \$REF -i calls/lumpy/0${n}sim.lumpy.vcf -o calls/lumpy/0${n}sim.lumpysvtyper.vcf" >> lumpy/${n}.svtyper.sh
done
chmod 700 lumpy/*.svtyper.sh
