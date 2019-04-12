for n in {0010..100..10} 200
do
   echo -n "BAMS=\"" >> t
   for i in $( eval echo {0001..$n} )
   do
      echo -n "--bam bam/S0${i}/S0${i}.bam " >> t
   done
   echo -n "\"" >> t
   echo -e "\nREF=\"--referenceFasta reference/GRCh38/genome.fa\"" >> t
   echo "WD=\"--runDir manta${n}\"" >> t
   echo -e "\npython manta/bin/configManta.py \$BAMS \$REF \$WD --region chr21"  >> t
   echo -e "\npython manta${n}/runWorkflow.py -m local -j 1" >> t
   cat t >> scripts/manta/${n}.manta.sh
   rm t
done
