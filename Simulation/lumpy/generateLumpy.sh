for n in {0010..0100..10} 200
do
   echo -n "BAM=\"" >> t
   for i in $( eval echo {00001..$n} )
   do
      echo -n "bam/S${i}/S${i}.bam," >> t
   done
   sed -i 's;,$;\";' t
   cat t >> scripts/lumpy/${n}.lumpy.sh
   rm t
   echo -e "\n\nlumpyexpress -B \$BAM -o lumpy/0${n}sim.lumpy.vcf" >> lumpy/${n}.lumpy.sh
done
