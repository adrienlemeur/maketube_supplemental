#!/usr/bin/env bash

all_strains=$(find ./ -type f \( -iname \*.fasta -o -iname \*.fa \))

mkdir -p VCF
mkdir -p F/F1 F/F2 F/F3 PSEUDOVCF
mkdir -p /home/alemeur/SCRATCH/

#rm -rf VCF_MPILEUP
mkdir -p VCF_MPILEUP VCF_MPILEUP/F0 VCF_MPILEUP/F1 VCF_MPILEUP/F2 VCF_MPILEUP/F3

ls BAM_30/*.bam | grep -v "nopg" | \
	parallel -j 10 '''export sample=$(basename {} .bam) ; echo $sample ; \
	bcftools mpileup -f REF/H37Rv.fasta BAM_30/${sample}.bam -ABq0 -Q0 -a DP,AD -Ob | \
	bcftools call --ploidy 2 -mg 10 -Oz -o VCF_MPILEUP/${sample}.vcf.gz'''


for one_strain in $all_strains
do
	sample=$(basename $one_strain ".fasta" | sed 's/\.fa//g')

	>&2 echo "doing $one_strain"

	if [ ! -f "VCF_MPILEUP/F3/${sample}.vcf.gz" ]
	then
		bcftools view "VCF_MPILEUP/F0/${sample}.normed.vcf.gz" -i "MAX(INFO/AF)>0.80 && DP4[3]+DP4[4]>=5 && QUAL>20" -Oz -o "VCF_MPILEUP/F3/${sample}.vcf.gz"
	fi
done
