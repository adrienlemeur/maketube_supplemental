#!/usr/bin/env bash

all_strains=$(find ./ -type f \( -iname \*.fasta -o -iname \*.fa \) | grep -v maketube_strains)
all_strains_maketube=$(find ./ -type f \( -iname \*.fasta -o -iname \*.fa \) | grep pop1 | grep H1 )

rm -rf VCF F
mkdir -p VCF /home/alemeur/SCRATCH/

mkdir -p F/F1 F/F2 F/F3 PSEUDOVCF
mkdir -p /home/alemeur/SCRATCH/

ls BAM_30/*.bam | parallel -j 10 'export sample=$(basename {} .bam) ; echo $sample ; freebayes -f "REF/H37Rv.fasta" "BAM_30/${sample}.bam" --vcf "VCF/${sample}.vcf" --ploidy 2 -C 4 -m 30 --report-all-haplotype-alleles'

echo -e "sample\tFILTER\tDP\tRO\tQR" > info_variant_calling_metrics.tsv
for one_strain in $all_strains $all_strains_maketube
do
	sample=$(basename $one_strain ".fasta" | sed 's/\.fa//g')

	>&2 echo "doing $one_strain"

	#mv -f ~/SCRATCH/tmp.vcf.gz "VCF/${sample}.vcf.gz"

	if [ ! -f "VCF/${sample}.vcf.gz" ]
	then
		bgzip "VCF/${sample}.vcf"
		bcftools norm -a -m- "VCF/${sample}.vcf.gz" -Oz -o /home/alemeur/SCRATCH/tmp.vcf.gz
		mv -f /home/alemeur/SCRATCH/tmp.vcf.gz "VCF/${sample}.vcf.gz"

		bcftools index -f "VCF/${sample}.vcf.gz"
	fi

	if [ ! -f "F/F3/${sample}.vcf.gz" ]
	then
		bcftools view "VCF/${sample}.vcf.gz" -i "MAX(INFO/AF)>0.80 && MAX(FORMAT/AD)>=5 && QUAL>20" -Oz -o "F/F3/${sample}.vcf.gz"
	fi

done

