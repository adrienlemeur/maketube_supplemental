#!/usr/bin/env bash

#Align the sample genome to the reference using nucmer and minimap2, calling the variant and converting the result to a vcf

#IMPORTANT : ALL2VCF to convert .snps to .vcf
#https://github.com/MatteoSchiavinato/all2vcf

all_strains=$(find ./ARTIFICIAL_STRAINS -type f \( -iname \*.fasta -o -iname \*.fa \))

REF="REF/H37Rv.fasta"

#rm -rf DELTAS PAF PAF_VCF PSEUDOVCF
mkdir -p DELTAS PAF PAF_VCF PSEUDOVCF PAF_VCF_SHORTER_ALIGNMENT ~/SCRATCH

for one_strain in $all_strains
do
	sample=$(basename $one_strain ".fasta" | sed 's/\.fa//g')

	if [ ! -f "PAF/${sample}.paf" ]
	then
		minimap2 -cx asm10 --cs "REF/H37Rv.fasta" $one_strain > "PAF/${sample}.paf"
	fi

	if [ ! -f "PAF_VCF/$sample.vcf.gz" ]
	then
		sort -k6,6 -k8,8n "PAF/${sample}.paf" | paftools.js call -f "REF/H37Rv.fasta" -s $sample - | \
		bcftools view -i '(ILEN >= -100 && ILEN <= 100) || TYPE!="INDEL"' > ~/SCRATCH/tmp.vcf && bgzip -f ~/SCRATCH/tmp.vcf
		bcftools norm -a -m- ~/SCRATCH/tmp.vcf.gz -Oz > "PAF_VCF/$sample.vcf.gz"
	fi

	#sort -k6,6 -k8,8n "PAF/${sample}.paf" | paftools.js call -L1000 -f "REF/H37Rv.fasta" -s $sample - | \
	#bcftools view -i '(ILEN >= -100 && ILEN <= 100) || TYPE!="INDEL"' > ~/SCRATCH/tmp.vcf && bgzip -f ~/SCRATCH/tmp.vcf
	#bcftools norm -a -m- ~/SCRATCH/tmp.vcf.gz -Oz > "PAF_VCF_SHORTER_ALIGNMENT/$sample.vcf.gz"

	#according to the nucmer documentation
	>&2 echo "doing $one_strain"
	if [ ! -f "DELTAS/$sample.delta" ]
	then
		nucmer -p "DELTAS/"$sample $REF $one_strain
	fi

	if [ ! -f "DELTAS/$sample.snps" ]
	then
		delta-filter -r -q "DELTAS/$sample.delta" > "DELTAS/$sample.delta2"
		show-snps "DELTAS/$sample.delta2" -T > "DELTAS/$sample.snps"
	fi

	if [ ! -f "PSEUDOVCF/$sample.vcf.gz" ]
	then
		all2vcf mummer --snps "DELTAS/$sample.snps" --input-header --output-header --reference "REF/H37Rv.fasta" > ~/SCRATCH/tmp.vcf && bgzip -f ~/SCRATCH/tmp.vcf
		bcftools norm -a -m- ~/SCRATCH/tmp.vcf.gz -Oz > "PSEUDOVCF/$sample.vcf.gz"
	fi


done




