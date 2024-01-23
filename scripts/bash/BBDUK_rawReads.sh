#!/usr/bin/env zsh

idList=$1
BBDUK=~/Programs/bbmap

if [ ! -f ${idList} ]
then
echo "id file does not exist, exiting"
exit 1
fi


for id in `cat ${idList}`
do
  echo "searching for ${id} read files in ../raw_reads/"
  if [ -f ../raw_reads/${id}_1.fastq.gz ] && [ -f ../raw_reads/${id}_2.fastq.gz ]
  then
    echo "read files for ${id} found!"
  else
    echo "read files for ${id} not found, exiting"
    exit 1
  fi
done

for id in `cat ${idList}`
do

  ${BBDUK}/bbduk.sh -Xmx1g \
  in1=../raw_reads/${id}_1.fastq.gz \
  in2=../raw_reads/${id}_2.fastq.gz \
  out1=../trimmed_reads/${id}_1.fastq.gz \
  out2=../trimmed_reads/${id}_2.fastq.gz \
  ref=${BBDUK}/resources/adapters.fa \
  ktrim=r \
  k=23 \
  mink=8 \
  hdist=1 \
  tpe \
  tbo \
  qtrim=lr \
  trimq=10 \
  minlen=50 \
  maq=10 \
  threads=auto

#  if [ -f ../trimmed_reads/${id}_1.fastq.gz ] && [ -f ../trimmed_reads/${id}_2.fastq.gz ]
#  then
#    rm ../raw_reads/${id}_1.fastq.gz
#    rm ../raw_reads/${id}_2.fastq.gz
#  fi

done
