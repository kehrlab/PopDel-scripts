#$ -S /bin/bash

#$ -P medium
#$ -pe smp 1
#$ -l h_rt=72:00:00
#$ -l h_vmem=16g
# Keep current environment variables, keeps virtualenv
#$ -V
#$ -j y
#$ -r yes
# Use current working directory as working directory of the job.  This
# requires you to cd/popd into the directory where the snake_job.sh lies.
#$ -cwd
#$ -o logs/lumpy

set -x

module purge

# Enforce existence of TMPDIR
export TMPDIR=tmp
mkdir -p $TMPDIR

# Activate bash cmd printing, debug info
set -x
>&2 hostname
>&2 date

source activate lumpy

bam="bam/NA12878_illu_platinum_GRCh38.bam"
vcfNoGt="calls/lumpy/platinum.lumpy.noGT.vcf"
vcf="calls/lumpy/platinum.lumpy.vcf"
ref="reference/hs38.fa"

lumpyexpress -B ${bam} -o ${vcfNoGt} > logs/lumpy/platinum.lumpy.log
svtyper -B ${bam} -T ${ref} -i ${vcfNoGt} -o ${vcf}  > logs/svtyper/platinum.svtyper.log

source deactivate

# Print date after finishing, for good measure
>&2 date
>&2 echo "All done. Have a nice day."
