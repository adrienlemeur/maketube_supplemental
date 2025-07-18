#!/usr/bin/env bash

#Align the sample genome to the reference using nucmer and minimap2, calling the variant and converting the result to a vcf

#IMPORTANT : ALL2VCF to convert .snps to .vcf
#https://github.com/MatteoSchiavinato/all2vcf

all_strains=$(find ./ -type f \( -iname \*.fasta -o -iname \*.fa \))

REF="REF/H37Rv.fasta"
mkdir -p DELTAS PAF PAF_VCF PSEUDOVCF ~/SCRATCH


while read sample_line
do

	sample=$(echo $sample_line | cut -f1 -d ' ')
	sample_file=$(echo $sample_line | cut -f3 -d ' ')
	sample_lineage=$(echo $sample_line | cut -f5 -d ' ')
	sample_source=$(echo $sample_line | cut -f6 -d ' ' | sed -E 's/.$//')

	>&2 echo "minimap2 $sample"
	if [ ! -f "PAF/${sample}.paf" ]
	then
		minimap2 -cx asm10 --cs "REF/H37Rv.fasta" $sample_file > "PAF/${sample}.paf"
	fi

	>&2 echo "PAF $sample"
	if [ ! -f "PAF_VCF/$sample.vcf.gz" ]
	then
		sort -k6,6 -k8,8n "PAF/${sample}.paf" | paftools.js call -f "REF/H37Rv.fasta" -s $sample - | \
		bcftools view -i '(ILEN >= -100 && ILEN <= 100) || TYPE!="INDEL"' > ~/SCRATCH/tmp.vcf && bgzip -f ~/SCRATCH/tmp.vcf
		bcftools norm -a -m- ~/SCRATCH/tmp.vcf.gz -Oz > "PAF_VCF/$sample.vcf.gz"
	fi

	>&2 echo "nucmer $sample"
	if [ ! -f "DELTAS/$sample.delta" ]
	then
		nucmer -p "DELTAS/"$sample "REF/H37Rv.fasta" $sample_file
	fi

	>&2 echo "snps $sample"
	if [ ! -f "DELTAS/$sample.snps" ]
	then
		delta-filter -r -q "DELTAS/$sample.delta" > "DELTAS/$sample.delta2"
		show-snps "DELTAS/$sample.delta2" -T > "DELTAS/$sample.snps"
	fi

	>&2 echo "all2vcf $sample"
	if [ ! -f "PSEUDOVCF/$sample.vcf.gz" ]
	then
		all2vcf mummer --snps "DELTAS/$sample.snps" --input-header --output-header --reference "REF/H37Rv.fasta" > ~/SCRATCH/tmp.vcf && bgzip -f ~/SCRATCH/tmp.vcf
		bcftools norm -a -m- ~/SCRATCH/tmp.vcf.gz -Oz > "PSEUDOVCF/$sample.vcf.gz"
	fi


done < <( grep -v '#' "natural_strains_info.txt")

