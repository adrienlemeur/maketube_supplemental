#!/usr/bin/env bash

all_strains=$(find ./ -type f \( -iname \*.fasta -o -iname \*.fa \))
REF="REF/H37Rv.fasta"

mkdir -p ~/SCRATCH

while read LINE
do
	sample=$(echo $LINE | cut -f1 -d ' ')
	FASTA=$(echo $LINE | cut -f2 -d ' ')
	
	source=$(echo $LINE | cut -f3 -d ' ')
	reference=$(echo $LINE | cut -f4 -d ' ')

	EQUIVALENCE=$(echo $LINE | cut -f5 -d ' ')
	annotation=$(echo $LINE | cut -f6 -d ' ')

	RAW_FREEBAYES=$(echo $LINE | cut -f7 -d ' ')
	GENOTUBE=$(echo $LINE | cut -f8 -d ' ')

	RAW_GATK=$(echo $LINE | cut -f9 -d ' ')
	MTBseq=$(echo $LINE | cut -f10 -d ' ')

	TBprofiler=$(echo -n $LINE | cut -f11 -d ' ')

	cat $annotation | grep "dupli" > ~/SCRATCH/dupli_region.bed


	bcftools norm -a -m- $reference -Oz > ~/SCRATCH/reference.vcf.gz #2>> /dev/null

	bcftools norm -a -m- $RAW_FREEBAYES -Oz > ~/SCRATCH/raw_freebayes.vcf.gz #2>> /dev/null


	vcf2metrics.py -i ~/SCRATCH/raw_freebayes.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--backtrack $EQUIVALENCE --sample $sample --add_col $source "freebayes" "RAW" "dupli"

	bcftools norm -a -m- $GENOTUBE -Oz > ~/SCRATCH/genotube_freebayes.vcf.gz #2>> /dev/null

	vcf2metrics.py -i ~/SCRATCH/genotube_freebayes.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--backtrack $EQUIVALENCE --sample $sample --add_col $source "freebayes" "genotube" "dupli"


	bcftools norm -a -m- $TBprofiler -Oz > ~/SCRATCH/tbprofiler.vcf.gz #2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/tbprofiler.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--backtrack $EQUIVALENCE --sample $sample --add_col $source "freebayes" "TBprofiler" "dupli"


	bcftools norm -a -m- $RAW_GATK -Oz > ~/SCRATCH/raw_gatk.vcf.gz #2>> /dev/null	
	vcf2metrics.py -i ~/SCRATCH/raw_gatk.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--backtrack $EQUIVALENCE --sample $sample --add_col $source "mpileup" "RAW" "dupli"

	bcftools norm -a -m- $MTBseq -Oz > ~/SCRATCH/mtbseq.vcf.gz #2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/mtbseq.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--backtrack $EQUIVALENCE --sample $sample --add_col $source "mpileup" "MTBseq" "dupli"

	cat $annotation | grep "duplicated" > ~/SCRATCH/dupli_region.bed

	#subtract
	vcf2metrics.py -i ~/SCRATCH/raw_freebayes.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--backtrack $EQUIVALENCE --subtract ~/SCRATCH/dupli_region.bed --sample $sample --add_col $source "freebayes" "RAW" "wo_dupli"

	vcf2metrics.py -i ~/SCRATCH/genotube_freebayes.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--backtrack $EQUIVALENCE --sample $sample --subtract ~/SCRATCH/dupli_region.bed --add_col $source "freebayes" "genotube" "wo_dupli"

	vcf2metrics.py -i ~/SCRATCH/tbprofiler.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--backtrack $EQUIVALENCE --sample $sample --subtract ~/SCRATCH/dupli_region.bed --add_col $source "freebayes" "TBprofiler" "wo_dupli"

	vcf2metrics.py -i ~/SCRATCH/raw_gatk.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--backtrack $EQUIVALENCE --sample $sample --subtract ~/SCRATCH/dupli_region.bed --add_col $source "mpileup" "RAW" "wo_dupli"

	vcf2metrics.py -i ~/SCRATCH/mtbseq.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--backtrack $EQUIVALENCE --sample $sample --subtract ~/SCRATCH/dupli_region.bed --add_col $source "mpileup" "MTBseq" "wo_dupli"


done < <(grep -v '#' a_database_to_rule_them_all.tsv | grep "maketube" | grep -e "H1\." -e "H10\.")

while read LINE
do
	sample=$(echo $LINE | cut -f1 -d ' ')
	FASTA=$(echo $LINE | cut -f2 -d ' ')
	
	source=$(echo $LINE | cut -f3 -d ' ')
	reference=$(echo $LINE | cut -f4 -d ' ')

	EQUIVALENCE=$(echo $LINE | cut -f5 -d ' ')
	annotation=$(echo $LINE | cut -f6 -d ' ')

	RAW_FREEBAYES=$(echo $LINE | cut -f7 -d ' ')
	GENOTUBE=$(echo $LINE | cut -f8 -d ' ')

	RAW_GATK=$(echo $LINE | cut -f9 -d ' ')
	MTBseq=$(echo $LINE | cut -f10 -d ' ')

	TBprofiler=$(echo -n $LINE | cut -f11 -d ' ')

	bcftools norm -a -m- $reference -Oz -o ~/SCRATCH/reference.vcf.gz #2>> /dev/null

	bcftools norm -a -m- $RAW_FREEBAYES -Oz > ~/SCRATCH/raw_freebayes.vcf.gz #2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/raw_freebayes.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--sample $sample --add_col $source "freebayes" "RAW" "dupli"

	bcftools norm -a -m- $GENOTUBE -Oz > ~/SCRATCH/genotube_freebayes.vcf.gz #2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/genotube_freebayes.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--sample $sample --add_col $source "freebayes" "genotube" "dupli"

	bcftools norm -a -m- $TBprofiler | bcftools view -v snps,mnps,indels -Oz > ~/SCRATCH/tbprofiler.vcf.gz #2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/tbprofiler.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--sample $sample --add_col $source "freebayes" "TBprofiler" "dupli"

	bcftools norm -a -m- $RAW_GATK -Oz > ~/SCRATCH/raw_gatk.vcf.gz #2>> /dev/null	
	vcf2metrics.py -i ~/SCRATCH/raw_gatk.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--sample $sample --add_col $source "mpileup" "RAW" "dupli"

	bcftools norm -a -m- $MTBseq -Oz > ~/SCRATCH/mtbseq.vcf.gz #2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/mtbseq.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--sample $sample --add_col $source "mpileup" "MTBseq" "dupli"

done < <(grep -v '#' a_database_to_rule_them_all.tsv | grep "snpmutator" | grep -e "H1\." -e "H10\.")

