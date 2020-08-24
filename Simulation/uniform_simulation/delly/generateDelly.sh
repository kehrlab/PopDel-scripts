for n in {0010..0100..10} {200..1000..100}
do
   echo -n "BAM=\"" >> t
   for i in $( eval echo {00001..$n} )
   do
      echo -n "bam/S${i}/S${i}.bam " >> t
   done
   echo -n "\"" >> t
   cat t >> delly/${n}.delly.sh
   rm t
   echo -e "\nREF=\"reference/GRCh38/genome.fa\"" >> delly/${n}.delly.sh
   echo "OUT=\"calls/delly/${n}sim.delly\"" >> delly/${n}.delly.sh
   echo -e "\ndelly call -n -o \${OUT} -g \${REF} \${BAM}" >> delly/${n}.delly.sh
done
chmod 774 delly/*.delly.sh
