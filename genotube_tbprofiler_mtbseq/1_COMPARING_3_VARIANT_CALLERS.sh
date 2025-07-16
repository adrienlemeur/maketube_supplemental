#!/usr/bin/env bash

REF="REF/H37Rv.fasta"

mkdir -p ~/SCRATCH

while read LINE
do

	sample=$(echo $LINE | cut -f1 -d ' ')
	strain=$(echo $LINE | cut -f2 -d ' ')

	fasta=$(echo $LINE | cut -f12 -d ' ')
	source=$(echo $LINE | cut -f3 -d ' ')

	reference_vcf=$(echo $LINE | cut -f11 -d ' ')

	maketube_equivalence=$(echo $LINE | cut -f4 -d ' ')
	maketube_structural_variants=$(echo $LINE | cut -f5 -d ' ')

	freebayes_raw=$(echo $LINE | cut -f6 -d ' ')
	genotube=$(echo $LINE | cut -f7 -d ' ')

	samtools_raw=$(echo $LINE | cut -f8 -d ' ')
	MTBseq=$(echo $LINE | cut -f9 -d ' ')

	TBprofiler=$(echo -n $LINE | cut -f10 -d ' ')
	
	fasta_length=$(cat $fasta | tr -d '\n' | grep -Pi "A|T|C|G|N" | wc -c)

	cat $maketube_structural_variants | grep "dupli" > ~/SCRATCH/dupli_region.bed

	bcftools norm -a -m- $reference_vcf -Ob | bcftools view -v snps -s $strain | grep -v 'AC=0' | bcftools view -Oz > ~/SCRATCH/reference.vcf.gz #2>> /dev/null

	bcftools norm -a -m- $samtools_raw -Oz > ~/SCRATCH/samtools_raw.vcf.gz #2>> /dev/null	
	bcftools norm -a -m- $freebayes_raw -Oz > ~/SCRATCH/raw_freebayes.vcf.gz #2>> /dev/null

	bcftools norm -a -m- $genotube -Oz > ~/SCRATCH/genotube.vcf.gz #2>> /dev/null
	bcftools norm -a -m- $MTBseq -Oz > ~/SCRATCH/mtbseq.vcf.gz #2>> /dev/null
	bcftools norm -a -m- $TBprofiler -Oz > ~/SCRATCH/tbprofiler.vcf.gz #2>> /dev/null

	python3 vcf2metrics.py -i ~/SCRATCH/raw_freebayes.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --backtrack $maketube_equivalence --sample $sample --add_col $source "freebayes" "freebayes_raw" "dupli" | grep -P "\tsnp"
	python3 vcf2metrics.py -i ~/SCRATCH/samtools_raw.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --backtrack $maketube_equivalence --sample $sample --add_col $source "samtools" "samtools_raw" "dupli" | grep -P "\tsnp"

	python3 vcf2metrics.py -i ~/SCRATCH/genotube.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --backtrack $maketube_equivalence --sample $sample --add_col $source "freebayes" "genotube" "dupli" | grep -P "\tsnp"
	python3 vcf2metrics.py -i ~/SCRATCH/mtbseq.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --backtrack $maketube_equivalence --sample $sample --add_col $source "samtools" "MTBseq" "dupli" | grep -P "\tsnp"
	python3 vcf2metrics.py -i ~/SCRATCH/tbprofiler.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --backtrack $maketube_equivalence --sample $sample --add_col $source "freebayes" "TBprofiler" "dupli" | grep -P "\tsnp"

	python3 vcf2metrics.py -i ~/SCRATCH/raw_freebayes.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --subtract ~/SCRATCH/dupli_region.bed --backtrack $maketube_equivalence --sample $sample --add_col $source "freebayes" "freebayes_raw" "wo_dupli" | grep -P "\tsnp"
	python3 vcf2metrics.py -i ~/SCRATCH/samtools_raw.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --subtract ~/SCRATCH/dupli_region.bed --backtrack $maketube_equivalence --sample $sample --add_col $source "samtools" "samtools_raw" "wo_dupli" | grep -P "\tsnp"

	python3 vcf2metrics.py -i ~/SCRATCH/genotube.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --subtract ~/SCRATCH/dupli_region.bed --backtrack $maketube_equivalence --sample $sample --add_col $source "freebayes" "genotube" "wo_dupli" | grep -P "\tsnp"
	python3 vcf2metrics.py -i ~/SCRATCH/mtbseq.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --subtract ~/SCRATCH/dupli_region.bed --backtrack $maketube_equivalence --sample $sample --add_col $source "samtools" "MTBseq" "wo_dupli" | grep -P "\tsnp"
	python3 vcf2metrics.py -i ~/SCRATCH/tbprofiler.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --subtract ~/SCRATCH/dupli_region.bed --backtrack $maketube_equivalence --sample $sample --add_col $source "freebayes" "TBprofiler" "wo_dupli" | grep -P "\tsnp"

done < <(grep -v '#' a_database_to_rule_them_all.tsv | grep "maketube")


while read LINE
do

	sample=$(echo $LINE | cut -f1 -d ' ')
	strain=$(echo $LINE | cut -f2 -d ' ')

	fasta=$(echo $LINE | cut -f12 -d ' ')
	source=$(echo $LINE | cut -f3 -d ' ')

	reference_vcf=$(echo $LINE | cut -f11 -d ' ')

	maketube_equivalence=$(echo $LINE | cut -f4 -d ' ')
	maketube_structural_variants=$(echo $LINE | cut -f5 -d ' ')

	freebayes_raw=$(echo $LINE | cut -f6 -d ' ')
	genotube=$(echo $LINE | cut -f7 -d ' ')

	samtools_raw=$(echo $LINE | cut -f8 -d ' ')
	MTBseq=$(echo $LINE | cut -f9 -d ' ')

	TBprofiler=$(echo -n $LINE | cut -f10 -d ' ')

	bcftools norm -a -m- $reference_vcf -Oz > ~/SCRATCH/reference.vcf.gz #2>> /dev/null

	bcftools norm -a -m- $samtools_raw -Oz > ~/SCRATCH/samtools_raw.vcf.gz #2>> /dev/null	
	bcftools norm -a -m- $freebayes_raw -Oz > ~/SCRATCH/raw_freebayes.vcf.gz #2>> /dev/null

	bcftools norm -a -m- $genotube -Oz > ~/SCRATCH/genotube.vcf.gz #2>> /dev/null
	bcftools norm -a -m- $MTBseq -Oz > ~/SCRATCH/mtbseq.vcf.gz #2>> /dev/null
	bcftools norm -a -m- $TBprofiler -Oz > ~/SCRATCH/tbprofiler.vcf.gz #2>> /dev/null

	python3 vcf2metrics.py -i ~/SCRATCH/raw_freebayes.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --sample $sample --add_col $source "freebayes" "freebayes_raw" "dupli" | grep -P "\tsnp"
	python3 vcf2metrics.py -i ~/SCRATCH/samtools_raw.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --sample $sample --add_col $source "samtools" "samtools_raw" "dupli" | grep -P "\tsnp"

	python3 vcf2metrics.py -i ~/SCRATCH/genotube.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --sample $sample --add_col $source "freebayes" "genotube" "dupli" | grep -P "\tsnp"
	python3 vcf2metrics.py -i ~/SCRATCH/mtbseq.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --sample $sample --add_col $source "samtools" "MTBseq" "dupli" | grep -P "\tsnp"
	python3 vcf2metrics.py -i ~/SCRATCH/tbprofiler.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --sample $sample --add_col $source "freebayes" "TBprofiler" "dupli" | grep -P "\tsnp"
done < <(grep -v '#' a_database_to_rule_them_all.tsv | grep "snpmutator")

