#!/bin/bash

cd data/sequence/alignment

# (step done in galaxy bcs too much here) 
# -- bwa-mem2 from galaxy takes reference and region8 -> bam file in data/sequence/alignment
bwa-mem2 mem -t 4 -R '@RG\tID:foo\tSM:bar' ../reference/GCF_018294505.1_genomic.fa ../../region8.fa > region8_aln.sam
samtools sort -@ 4 -o region8_aln.bam region8_aln.sam

samtools view region8_aln.bam #just to check the bam file
bedtools bamtobed -i region8_aln.bam > region8_aln.bed