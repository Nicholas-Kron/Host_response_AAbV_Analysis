#!/usr/bin/env bash

## ANV Genome = 36420
## min(14, log2(GenomeLength)/2 - 1)
## log2(36420)/2 - 1 = 15/2 - 1 = 7.5 - 1 = 6.5

if ![-d ../genomes/ANVIndex]
then
mkdir ../genomes/ANVIndex
fi

star --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ../genomes/ANVIndex \
--genomeFastaFiles ../genomes/ANV_genome.fasta \
--sjdbGTFfile ../genomes/ANV_genome.gff3 \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfeatureExon CDS \
--sjdbOverhang 99 \
--genomeSAindexNbases 6
