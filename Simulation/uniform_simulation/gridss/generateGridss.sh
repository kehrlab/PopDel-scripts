for n in {0010..0100..10} {0200..500..100}
do
   echo -n "BAM=\"" >> t
   for i in $( eval echo {00001..$n} )
   do
      echo -n "INPUT=bam/S${i}/S${i}.bam " >> t
   done
   echo -n "\"" >> t
   cat t >> scripts/gridss/${n}.gridss.sh
   rm t
   echo -e "\nASM=\"ASSEMBLY=tmp/${n}.assembly.bam\"" >> scripts/gridss/${n}.gridss.sh
   echo "WD=\"WORKING_DIR=tmp\"" >> scripts/gridss/${n}.gridss.sh
   echo "TD=\"TMP_DIR=gridsstmp\"" >> scripts/gridss/${n}.gridss.sh
   echo "REF=\"REFERENCE_SEQUENCE=reference/GRCh38/genome.fa\""  >> scripts/gridss/${n}.gridss.sh
   echo "OUT=\"OUTPUT=calls/gridss/${n}sim.gridss.vcf\"" >> scripts/gridss/${n}.gridss.sh

   echo -e "\njava -jar -Xmx8g gridss/gridss-1.8.1-gridss-jar-with-dependencies.jar \${BAM} \${ASM} \${WD} \${TD} \${OUT} \${REF} WORKER_THREADS=1" >> scripts/gridss/${n}.gridss.sh

   echo -e "\nrm -rf gridsstmp/*" >> scripts/gridss/${n}.gridss.sh
   echo "rm -rf tmp/*" >> scripts/gridss/${n}.gridss.sh
done
chmod 700 scripts/gridss/*.gridss.sh
