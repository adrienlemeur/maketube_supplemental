	#!/usr/bin/env bash

mkdir -p ~/SCRATCH

REF="REF/H37Rv.fasta"


while read sample_line
do
	sample=$(echo $sample_line | cut -f1 -d ' ')
	sample_file=$(echo $sample_line | cut -f3 -d ' ')
	sample_lineage=$(echo $sample_line | cut -f5 -d ' ')
	sample_source=$(echo $sample_line | cut -f6 -d ' ')

	fasta_length=$(cat $sample_file | grep -E -o "G|C|T|A|N" | wc -l)

	EQUIVALENCE=$(echo $sample_line | cut -f7 -d ' ')

	reference=$(echo $sample_line | cut -f9 -d ' ')

	minimap2=$(echo $sample_line | cut -f10 -d ' ')
	nucmer=$(echo $sample_line | cut -f11 -d ' ')
	strain=$(echo $sample_line | cut -f12 -d ' ')

	bcftools norm -a -m- $reference | bcftools view -v snps -Oz > ~/SCRATCH/reference.vcf.gz
	bcftools norm -a -m- $minimap2 | bcftools view -v snps -Oz > ~/SCRATCH/PAF_VCF.vcf.gz 2>> /dev/null
	bcftools norm -a -m- $nucmer | bcftools view -v snps -Oz > ~/SCRATCH/PSEUDOVCF.vcf.gz 2>> /dev/null

	vcf2metrics.py --sample $sample -i ~/SCRATCH/PSEUDOVCF.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --add_col $sample $sample_source "w_duplicated" "nucmer" 2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/PAF_VCF.vcf.gz --sample $sample --reference ~/SCRATCH/reference.vcf.gz --add_col $sample $sample_source "w_duplicated" "minimap2" 2>> /dev/null

done < <( grep -v '#' "natural_strains_info.txt" | grep "snpmutator")

while read sample_line
do
	sample=$(echo $sample_line | cut -f1 -d ' ')
	sample_file=$(echo $sample_line | cut -f3 -d ' ')
	sample_lineage=$(echo $sample_line | cut -f5 -d ' ')
	sample_source=$(echo $sample_line | cut -f6 -d ' ')

	fasta_length=$(cat $sample_file | grep -E -o "G|C|T|A|N" | wc -l)

	EQUIVALENCE=$(echo $sample_line | cut -f7 -d ' ')

	reference=$(echo $sample_line | cut -f9 -d ' ')

	minimap2=$(echo $sample_line | cut -f10 -d ' ')
	nucmer=$(echo $sample_line | cut -f11 -d ' ')
	strain=$(echo $sample_line | cut -f12 -d ' ')

	echo | awk -v fasta_lgth=$fasta_length '{print "H37Rv\t"fasta_lgth-150000"\t"fasta_lgth"\tduplicata_region"}' >> ~/SCRATCH/dupli_region.bed

	bcftools norm -a -m- $reference | bcftools view -v snps -s $strain -Ov | grep -v "	0:" > ~/SCRATCH/reference.vcf 2>> /dev/null && bgzip -f ~/SCRATCH/reference.vcf
	bcftools norm -a -m- $minimap2 | bcftools view -v snps -Oz > ~/SCRATCH/PAF_VCF.vcf.gz 2>> /dev/null
	bcftools norm -a -m- $nucmer | bcftools view -v snps -Oz > ~/SCRATCH/PSEUDOVCF.vcf.gz 2>> /dev/null

	vcf2metrics.py --sample $sample -i ~/SCRATCH/PSEUDOVCF.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --subtract ~/SCRATCH/dupli_region.bed --backtrack $EQUIVALENCE --add_col $sample $sample_source "wo_duplicated" "nucmer" 2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/PAF_VCF.vcf.gz --sample $sample --reference ~/SCRATCH/reference.vcf.gz --subtract ~/SCRATCH/dupli_region.bed --backtrack $EQUIVALENCE --add_col $sample $sample_source "wo_duplicated" "minimap2" 2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/PSEUDOVCF.vcf.gz --sample $sample --reference ~/SCRATCH/reference.vcf.gz --backtrack $EQUIVALENCE --add_col $sample $sample_source "w_duplicated" "nucmer" 2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/PAF_VCF.vcf.gz --sample $sample --reference ~/SCRATCH/reference.vcf.gz --backtrack $EQUIVALENCE --add_col $sample $sample_source "w_duplicated" "minimap2" 2>> /dev/null

done < <( grep -v '#' "natural_strains_info.txt" | grep "maketube")

