#!/usr/bin/env bash

idList=$1

if [ ! -f ${idList} ]
then
echo "id file does not exist, exiting"
exit 1
fi

ulimit -n 5000

for id in `cat ${idList}`
do
  if [ ! -f ../unmapped_reads/${id}_Unmapped.out.mate1.gz -o ! -f ../unmapped_reads/${id}_Unmapped.out.mate2.gz ]
  then
  echo "read files for ${id} do not exist, exiting"
  exit 1
else
  echo "read files for ${id} found"
  fi
done

for id in `cat ${idList}`
do

  echo "processing ${id}"

  star --runThreadN 4 \
  --genomeDir ../genomes/ANVIndex \
  --readFilesIn ../unmapped_reads/${id}_Unmapped.out.mate1.gz \
  ../unmapped_reads/${id}_Unmapped.out.mate2.gz \
  --outFileNamePrefix ../anv_aligned/${id}_ANV_ \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 5000000000 \
  --quantMode GeneCounts \
  --readFilesCommand gunzip -c \
  --outSAMstrandField intronMotif \
  --outSAMattributes NH HI AS nM NM MD jM jI XS MC ch


  echo "${id} aligned, output in ~/projects/ANV/anv_aligned/${id}_ANV"

  echo "removing unampped reads from ${id}"

  samtools view -b -F 4 ../anv_aligned/${id}_ANV_Aligned.sortedByCoord.out.bam  > \
  ../anv_aligned/${id}_ANV_Aligned.sortedByCoord.onlyMapped.bam

  if [ -f ../anv_aligned/${id}_ANV_Log.final.out ]
  then
    #rm ../anv_aligned/${id}_ANV_Aligned.sortedByCoord.out.bam
    mv ../anv_aligned/${id}_ANV_Log.final.out ../log_files
    mv ../anv_aligned/${id}_ANV_Log.out ../log_files
    mv ../anv_aligned/${id}_ANV_Log.progress.out ../log_files
    mv ../anv_aligned/${id}_ANV_ReadsPerGene.out.tab ../log_files
    mv ../anv_aligned/${id}_ANV_SJ.out.tab ../log_files
  fi

done
