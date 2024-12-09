#!/usr/bin/env bash

mkdir -p BAM FILTERED_VCF MPILEUP_FILTERED_VCF

if [[ "$1" == "--filter" ]]
then
	for vcf in VCF/*.vcf.gz
	do
		sample=$(basename $vcf .vcf.gz)
		bcftools index -f $vcf
		bcftools norm -a -m- $vcf -Oz -o ~/SCRATCH/tmp.vcf.gz

		bcftools view ~/SCRATCH/tmp.vcf.gz -i "MAX(INFO/AF)>0.80 && MAX(FORMAT/AD)>=5 && QUAL>20" -Oz -o FILTERED_VCF/${sample}.vcf.gz
	done

	for vcf in VCF_MPILEUP/*.vcf.gz
	do
		sample=$(basename $vcf .vcf.gz)
		bcftools index -f $vcf
		bcftools norm -a -m- $vcf -Oz -o ~/SCRATCH/tmp.vcf.gz

		bcftools view ~/SCRATCH/tmp.vcf.gz -v snps,indels,mnps -Oz -o MPILEUP_FILTERED_VCF/${sample}.vcf.gz
	done

fi

#rm -rf region_sizes.tsv
while read LINE
do
	sample=$(echo $LINE | cut -f1 -d ' ')
	#echo "Doing... $sample"

	fasta_gz=$(echo $LINE | cut -f9 -d ' ')

	fasta=$(echo $fasta_gz | sed 's/\.gz//g')

	R1=$(echo $LINE | cut -f10 -d ' ')
	R2=$(echo $LINE | cut -f11 -d ' ')
	fastq=$(echo $R1 | sed 's/1.fq.gz//g')

	reference_vcf=$(echo $LINE | cut -f8 -d ' ')
	test_vcf=$(echo $LINE | cut -f12 -d ' ')

	equivalence=$(echo $LINE | cut -f6 -d ' ')

	strain_regions=$(echo $LINE | cut -f7 -d ' ')

	H37Rv_regions="REF/all_regions.bed"
	antibiotic_resistance="REF/antibiotic_resistance.bed"

	#"REF/H37Rv_non_independant_region.bed"
	cat $strain_regions $antibiotic_resistance $H37Rv_regions | \
		grep -v "NEUTRAL" | sed 's/NC_000962.3/H37Rv/g' | \
		grep -v "duplicata_region" | \
		awk '{if($4 == "insertion"){print($1, $2-250, $3+250, $4)} else {print($0)} }' | sed 's/ /\t/g' | sort -k1,1 -k2,2n > ~/SCRATCH/${sample}.tmp

	bedtools window \
		-a <(sort -k1,1 -k2,2n $strain_regions) \
		-b <(grep -v NEUTRAL REF/all_regions.bed | sort -k1,1 -k2,2) \
		-w 50 -u > ~/SCRATCH/${sample}.SV_close_region.bed
	bedtools window \
		-a <(sort -k1,1 -k2,2n $strain_regions) \
		-b <(grep -v NEUTRAL REF/all_regions.bed | sort -k1,1 -k2,2) \
		-w 50 -v > ~/SCRATCH/${sample}.SV_not_close_region.bed


#	vcf2metrics.py --sample $sample -i FILTERED_VCF/${sample}.vcf.gz \
#		--pipeline "genotube" \
#		--SV "NEAR" \
#		--reference ~/SCRATCH/${sample}.vcf.gz \
#		--backtrack $equivalence \
#		--bed <(sort -k1,1 -k2,2n ~/SCRATCH/${sample}.SV_close_region.bed)

	vcf2metrics.py --sample $sample -i FILTERED_VCF/${sample}.vcf.gz \
		--pipeline "TBprofiler" \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence \
		--bed ~/SCRATCH/${sample}.tmp

done < <(grep -v 'strain' maketube_run/my_run_arborescence.tsv | grep -v "NA" )




