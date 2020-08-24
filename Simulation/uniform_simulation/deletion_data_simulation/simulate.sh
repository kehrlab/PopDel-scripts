#!/bin/bash

#$ -pe smp 4
#$ -l h_vmem=8g
#$ -l h_rt=12:00:00
#$ -cwd
#$ -P medium
#$ -o sge_logs
#$ -e sge_logs

set -e -o pipefail

source ~/.bashrc

#SAMPLE=$1
#RANDOM_SEED=$2

CHR21=reference/hg38/Chromosomes/chr21/chr21.fa
GENOME=reference/hg38/BWA_index/genome.fa
VARIANTS=simulated/2018-04-05_HiSeq-deletions_uniform/deletions.txt

ART=bin/art_illumina
SIMPY=deletion_data_simulation/simulate_deletion_haplotype.py

echo "mkdir -p ${SAMPLE} && cd ${SAMPLE}"
mkdir -p ${SAMPLE} && cd ${SAMPLE}
echo "mkdir -p logs"
mkdir -p logs

SEED2=$(( ${RANDOM_SEED} * 2 ))
SEED1=$(( ${SEED2} - 1 ))

PREFIX1=${SAMPLE}-hap1
PREFIX2=${SAMPLE}-hap2
SIMREF1=${PREFIX1}.chr21.fa
SIMREF2=${PREFIX2}.chr21.fa
RG1=@RG\\tID:${SAMPLE}\\tLB:${SAMPLE}\\tSM:${SAMPLE}\\tPL:ART-HS25
RG2=@RG\\tID:${SAMPLE}\\tLB:${SAMPLE}\\tSM:${SAMPLE}\\tPL:ART-HS25


### Create the haplotype reference sequence.

echo "python ${SIMPY} ${SEED1} ${VARIANTS} ${PREFIX1}.deletions.txt ${CHR21} ${SIMREF1} 1> logs/make-hap1.out 2> logs/make-hap1.err"
python ${SIMPY} ${SEED1} ${VARIANTS} ${PREFIX1}.deletions.txt ${CHR21} ${SIMREF1} 1> logs/make-hap1.out 2> logs/make-hap1.err &
echo "python ${SIMPY} ${SEED2} ${VARIANTS} ${PREFIX2}.deletions.txt ${CHR21} ${SIMREF2} 1> logs/make-hap2.out 2> logs/make-hap2.err"
python ${SIMPY} ${SEED2} ${VARIANTS} ${PREFIX2}.deletions.txt ${CHR21} ${SIMREF2} 1> logs/make-hap2.out 2> logs/make-hap2.err &
wait


### Run the simulation.

echo "${ART} -ss HS25 -i ${SIMREF1} -p -l 150 -f 15 -m 300 -s 50 -rs ${SEED1} -o ${PREFIX1}. 1> logs/sim-hap1.out 2> logs/sim-hap1.err"
${ART} -ss HS25 -i ${SIMREF1} -p -l 150 -f 15 -m 300 -s 50 -rs ${SEED1} -o ${PREFIX1}. 1> logs/sim-hap1.out 2> logs/sim-hap1.err &
echo "${ART} -ss HS25 -i ${SIMREF2} -p -l 150 -f 15 -m 300 -s 50 -rs ${SEED2} -o ${PREFIX2}. 1> logs/sim-hap2.out 2> logs/sim-hap2.err"
${ART} -ss HS25 -i ${SIMREF2} -p -l 150 -f 15 -m 300 -s 50 -rs ${SEED2} -o ${PREFIX2}. 1> logs/sim-hap2.out 2> logs/sim-hap2.err &
wait

echo "rm ${PREFIX1}.1.aln ${PREFIX1}.2.aln ${PREFIX2}.1.aln ${PREFIX2}.2.aln"
rm ${PREFIX1}.1.aln ${PREFIX1}.2.aln ${PREFIX2}.1.aln ${PREFIX2}.2.aln

echo "${SIMREF1} ${SIMREF2}"
rm ${SIMREF1} ${SIMREF2}


### Align the simulated reads to the reference genome.

echo "bwa mem -t 4 -M -R ${RG1} ${GENOME} ${PREFIX1}.1.fq ${PREFIX1}.2.fq 2> logs/bwa-hap1.err | samtools view -Sb - | samtools sort -@ 4 -o ${PREFIX1}.aligned.sorted.bam"
bwa mem -t 4 -M -R ${RG1} ${GENOME} ${PREFIX1}.1.fq ${PREFIX1}.2.fq 2> logs/bwa-hap1.err | samtools view -Sb - | samtools sort -@ 4 -o ${PREFIX1}.aligned.sorted.bam
echo "bwa mem -t 4 -M -R ${RG2} ${GENOME} ${PREFIX2}.1.fq ${PREFIX2}.2.fq 2> logs/bwa-hap2.err | samtools view -Sb - | samtools sort -@ 4 -o ${PREFIX2}.aligned.sorted.bam"
bwa mem -t 4 -M -R ${RG2} ${GENOME} ${PREFIX2}.1.fq ${PREFIX2}.2.fq 2> logs/bwa-hap2.err | samtools view -Sb - | samtools sort -@ 4 -o ${PREFIX2}.aligned.sorted.bam

echo "rm ${PREFIX1}.1.fq ${PREFIX1}.2.fq ${PREFIX2}.1.fq ${PREFIX2}.2.fq"
rm ${PREFIX1}.1.fq ${PREFIX1}.2.fq ${PREFIX2}.1.fq ${PREFIX2}.2.fq

echo "samtools merge ${SAMPLE}.bam ${PREFIX1}.aligned.sorted.bam ${PREFIX2}.aligned.sorted.bam"
samtools merge -c ${SAMPLE}.bam ${PREFIX1}.aligned.sorted.bam ${PREFIX2}.aligned.sorted.bam

echo "rm ${PREFIX1}.aligned.sorted.bam ${PREFIX2}.aligned.sorted.bam"
rm ${PREFIX1}.aligned.sorted.bam ${PREFIX2}.aligned.sorted.bam &
echo "samtools index ${SAMPLE}.bam"
samtools index ${SAMPLE}.bam
wait
