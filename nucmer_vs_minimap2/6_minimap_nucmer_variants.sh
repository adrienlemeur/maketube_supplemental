#!/usr/bin/env bash

#Align the sample genome to the reference using nucmer and minimap2, calling the variant and converting the result to a vcf

#IMPORTANT : ALL2VCF to convert .snps to .vcf
#https://github.com/MatteoSchiavinato/all2vcf

all_strains=$(find ./ -type f \( -iname \*.fasta -o -iname \*.fa \))

REF="REF/H37Rv.fasta"
mkdir -p DELTAS PAF PAF_VCF PSEUDOVCF ~/SCRATCH

while read line
do
	sample=$(echo $line | cut -f1 -d ' ')
	fasta=$(echo $line | cut -f9 -d ' ')

	unzipped_fasta=$(echo $fasta | sed 's/\.gz//g')

	>&2 echo "gunzip $sample"
	if [ ! -f $unzipped_fasta ]
	then
		zcat $fasta > $unzipped_fasta
	fi

	>&2 echo "minimap2 $sample"
	if [ ! -f "PAF/${sample}.paf" ]
	then
		minimap2 -cx asm10 --cs "REF/H37Rv.fasta" $unzipped_fasta > "PAF/${sample}.paf"
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
		nucmer -p "DELTAS/"$sample $REF $unzipped_fasta
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



done < <(grep -v 'strain' "maketube_run/H37Rv_run/my_run_arborescence.tsv")

exit



#rm -rf DELTAS PAF PAF_VCF PSEUDOVCF
mkdir -p DELTAS PAF PAF_VCF PSEUDOVCF ~/SCRATCH

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




