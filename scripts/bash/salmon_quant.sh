#!/usr/bin/env bash

gffread -w output.fa -g genome.fa genome.gtf

idList=$1

if [ ! -f ${idList} ]
then
echo "id file does not exist, exiting"
exit 1
fi

for id in `cat ${idList}`
do
  if [ ! -f ../aligned/${id}_Aligned.toTranscriptome.out.bam ]
  then
  echo "bam file for ${id} do not exist, exiting"
  exit 1
else
  echo "bam file for ${id} found"
  fi
done

conda activate salmon

for id in `cat ${idList}`
do
echo "Processing sample ${id}"
salmon quant -l ISR \
         -t ~/Desktop/genomes/AplCal3.0/AplCal_gffread_RNA.fa \
         -a ../aligned/${id}_Aligned.toTranscriptome.out.bam \
         -p 8 \
         -o ../quants/${id}_quant \
         --seqBias \
         --gcBias
done
