#!/usr/bin/env bash

mkdir -p ~/SCRATCH

REF="REF/H37Rv.fasta"

while read LINE
do
	sample=$(echo $LINE | cut -f1 -d ' ')
	source=$(echo $LINE | cut -f2 -d ' ')
	reference=$(echo $LINE | cut -f3 -d ' ')
	EQUIVALENCE=$(echo $LINE | cut -f4 -d ' ')
	BED=$(echo $LINE | cut -f5 -d ' ')
	PSEUDOVCF=$(echo $LINE | cut -f6 -d ' ')
	PAF_VCF=$(echo $LINE | cut -f7 -d ' ')

	FASTA=$(echo $LINE | cut -f15 -d ' ' | sed -E 's/.$//')
	fasta_length=$(cat $FASTA | grep -E -o "G|C|T|A|N" | wc -l)

	#cat $BED | grep "duplicated" > ~/SCRATCH/duplicated_region.bed
	cat $BED | grep "dupli" > ~/SCRATCH/dupli_region.bed
	echo | awk -v fasta_lgth=$fasta_length '{print "H37Rv\t"fasta_lgth-150000"\t"fasta_lgth"\tduplicata_region"}' >> ~/SCRATCH/dupli_region.bed

	bcftools norm -a -m- $reference | bcftools view -v snps -Oz > ~/SCRATCH/reference.vcf.gz 2>> /dev/null
	bcftools norm -a -m- $PAF_VCF | bcftools view -v snps -Oz > ~/SCRATCH/PAF_VCF.vcf.gz 2>> /dev/null
	bcftools norm -a -m- $PSEUDOVCF | bcftools view -v snps -Oz > ~/SCRATCH/PSEUDOVCF.vcf.gz 2>> /dev/null

	#missing some files for samples SV8_pop1_H10, SV6_pop1_H10, SV2_pop1_H10 so I removed them
	vcf2metrics.py -i ~/SCRATCH/PSEUDOVCF.vcf.gz --sample $sample --reference ~/SCRATCH/reference.vcf.gz --subtract ~/SCRATCH/dupli_region.bed --backtrack $EQUIVALENCE --add_col $sample $source "wo_duplicated" "nucmer" 2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/PAF_VCF.vcf.gz --sample $sample --reference ~/SCRATCH/reference.vcf.gz --subtract ~/SCRATCH/dupli_region.bed --backtrack $EQUIVALENCE --add_col $sample $source "wo_duplicated" "minimap2" 2>> /dev/null

	vcf2metrics.py -i ~/SCRATCH/PSEUDOVCF.vcf.gz --sample $sample --reference ~/SCRATCH/reference.vcf.gz --backtrack $EQUIVALENCE --add_col $sample $source "w_duplicated" "nucmer" 2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/PAF_VCF.vcf.gz --sample $sample --reference ~/SCRATCH/reference.vcf.gz --backtrack $EQUIVALENCE --add_col $sample $source "w_duplicated" "minimap2" 2>> /dev/null

done < <(grep -v '#' a_database_to_rule_them_all.tsv | grep "maketube_strains" | grep -e "H1\." -e "H10\." )

while read LINE
do
	sample=$(echo $LINE | cut -f1 -d ' ')
	source=$(echo $LINE | cut -f2 -d ' ')
	reference=$(echo $LINE | cut -f3 -d ' ')
	EQUIVALENCE=$(echo $LINE | cut -f4 -d ' ')
	BED=$(echo $LINE | cut -f5 -d ' ')
	PSEUDOVCF=$(echo $LINE | cut -f6 -d ' ')
	PAF_VCF=$(echo $LINE | cut -f7 -d ' ')

	bcftools norm -a -m- $reference -Oz -o ~/SCRATCH/reference.vcf.gz 2>> /dev/null

	bcftools norm -a -m- $PAF_VCF -Oz -o ~/SCRATCH/PAF_VCF.vcf.gz 2>> /dev/null
	bcftools norm -a -m- $PSEUDOVCF -Oz -o ~/SCRATCH/PSEUDOVCF.vcf.gz 2>> /dev/null

	vcf2metrics.py --sample $sample -i ~/SCRATCH/PSEUDOVCF.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --add_col $sample $source "no_duplicated" "nucmer"
	vcf2metrics.py --sample $sample -i ~/SCRATCH/PAF_VCF.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --add_col $sample $source "no_duplicated" "minimap2"

done < <(grep -v '#' a_database_to_rule_them_all.tsv | grep "snpmutator" )


