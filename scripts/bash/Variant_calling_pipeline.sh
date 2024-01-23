#! /usr/bin/env bash

SraAccList=$1
BAMdir="../../data/bams"
REF="../../data/modified_ANV_genome.fasta"

[ ! -d ${BAMdir}/GATK ] && mkdir ${BAMdir}/GATK

#make indecies for ref fasta
gatk CreateSequenceDictionary \
-R ${REF}

samtools faidx ${REF}

samples=`cat ${SraAccList}`

for samp in ${samples}
do

[[ ! -f "${BAMdir}/${samp}_Aligned.sortedByCoord.out.bam" ]] && (printf "BAM file for "${samp}" not found!" && exit )


echo "marking duplicates with picard ${samp}"
gatk MarkDuplicates \
-I ${BAMdir}/${samp}_Aligned.sortedByCoord.out.bam \
-M ${BAMdir}/GATK/${samp}_metrics.txt \
-O ${BAMdir}/GATK/${samp}_Dedupped.bam \

echo "adding read groups with picard ${samp}"
gatk AddOrReplaceReadGroups \
-I ${BAMdir}/GATK/${samp}_Dedupped.bam \
-O ${BAMdir}/GATK/${samp}_RG.bam \
-LB Lib1 \
-PL Illumina \
-PU unit1 \
-SM ${samp}

echo "building bam index for ${samp}"
gatk BuildBamIndex \
-I ${BAMdir}/GATK/${samp}_RG.bam \

echo "splitting N Cigar reads for ${samp}"
gatk SplitNCigarReads \
      -R ${REF} \
      -I ${BAMdir}/GATK/${samp}_RG.bam \
      -O ${BAMdir}/GATK/${samp}_split.bam

echo "recallibrating reads for ${samp}"
gatk BaseRecalibrator \
-I ${BAMdir}/GATK/${samp}_split.bam \
-R ${REF} \
-O ${BAMdir}/GATK/${samp}_recal_data.table \
--known-sites ../../data/annotation_data/BED12_modified_ANV_genome_for_var.bed

echo "Applying base recalibration for ${samp}"
gatk ApplyBQSR \
  -R ${REF} \
  -I ${BAMdir}/GATK/${samp}_split.bam \
  --bqsr-recal-file ${BAMdir}/GATK/${samp}_recal_data.table \
  -O ${BAMdir}/GATK/${samp}_recal.bam

echo "marking variants with Mutect2 for ${samp}"
gatk Mutect2 \
  -R ${REF} \
  -I ${BAMdir}/GATK/${samp}_recal.bam \
  -O ${BAMdir}/GATK/${samp}_variants.vcf.gz \
  --mitochondria-mode

  echo "filtering mutec calls ${samp}"
  gatk FilterMutectCalls \
  -R ${REF} \
  -V ${BAMdir}/GATK/${samp}_variants.vcf.gz \
  -O ${BAMdir}/GATK/${samp}_variants_filtered_mutec.vcf.gz

#  echo "filtering variants for ${samp}"
#  gatk VariantFiltration \
#    -R ${REF} \
#    -V ${BAMdir}/GATK/${samp}_variants_filtered_mutec.vcf.gz \
#    -O ${BAMdir}/GATK/${samp}_variants_filtered.vcf.gz \
#   --filter-expression "AB < 0.2 || MQ0 > 50" \
#   --filter-name "my_filters"

done
