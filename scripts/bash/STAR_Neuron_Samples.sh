#!/usr/bin/env bash

idList=$1

if [ ! -f ${idList} ]
then
echo "id file does not exist, exiting"
exit 1
fi

for id in `cat ${idList}`
do
  if [ ! -f ../trimmed_reads/${id}_1.fastq.gz -o ! -f ../trimmed_reads/${id}_2.fastq.gz ]
  then
  echo "trimmed read files for ${id} do not exist, exiting"
  exit 1
else
  echo "trimmed read files for ${id} found"
  fi
done


ulimit -n 5000

for id in `cat ${idList}`
do

  echo "processing ${id}"

  star --runThreadN 8 \
  --genomeDir ~/Desktop/genomes/AplCal3.0/AplCal3_STAR_Index \
  --readFilesIn ../trimmed_reads/${id}_1.fastq.gz \
  ../trimmed_reads/${id}_2.fastq.gz \
  --outFileNamePrefix ../aligned/${id}_ \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 5000000000 \
  --outReadsUnmapped Fastx \
  --quantMode TranscriptomeSAM GeneCounts\
  --twopassMode Basic \
  --readFilesCommand gunzip -c

  if [ -f ../aligned/${id}_Log.final.out ]
  then
    rm ../trimmed_reads/${id}_1.fastq.gz
    rm ../trimmed_reads/${id}_2.fastq.gz
  fi

  #rm ../aligned/${id}_Aligned.sortedByCoord.out.bam
  mv ../aligned/${id}_Unmapped.out.mate1 ../unmapped_reads
  mv ../aligned/${id}_Unmapped.out.mate2 ../unmapped_reads
  mv ../aligned/${id}_Log.final.out ../log_files
  mv ../aligned/${id}_Log.out ../log_files
  mv ../aligned/${id}_Log.progress.out ../log_files
  mv ../aligned/${id}_ReadsPerGene.out.tab ../log_files
  mv ../aligned/${id}_SJ.out.tab ../log_files

done
