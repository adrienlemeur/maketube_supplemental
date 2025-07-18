#!/usr/bin/env bash

mkdir -p VCF FILTERED_VCF

#variant calling
ls BAM/*.bam | parallel -j 10 'export sample=$(basename {} .bam) ; echo $sample ; freebayes -f "REF/H37Rv.fasta" "BAM/${sample}.bam" --vcf "VCF/${sample}.vcf" --ploidy 2 -C 4 -m 30 --report-all-haplotype-alleles && bgzip VCF/${sample}.vcf'

#older experience
#mkdir -p VCF_MPILEUP
#ls BAM/*.bam | grep -v "nopg" | \
#	parallel -j 10 '''export sample=$(basename {} .bam) ; echo $sample ; \
#	bcftools mpileup -f REF/H37Rv.fasta BAM/${sample}.bam -ABq0 -Q0 -a DP,AD -Ob | \
#	bcftools call --ploidy 2 -mg 10 -Oz -o VCF_MPILEUP/${sample}.vcf.gz'''
