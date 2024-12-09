#!/usr/bin/env bash

all_strains=$(find ./ -type f \( -iname \*.fasta -o -iname \*.fa \) | grep -v maketube_strains)
all_strains_maketube=$(find ./ -type f \( -iname \*.fasta -o -iname \*.fa \) | grep pop1 | grep -e "H1\." -e "H10\." )

REF="REF/H37Rv.fasta"

#rm -rf DELTAS PAF PAF_VCF PSEUDOVCF
mkdir -p DELTAS PAF PAF_VCF PSEUDOVCF

for one_strain in $all_strains $all_strains_maketube
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

	#>&2 echo "doing $one_strain"
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




