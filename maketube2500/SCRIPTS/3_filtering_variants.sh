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

	cat $strain_regions $antibiotic_resistance $H37Rv_regions | \
		grep -v "NEUTRAL" | sed 's/NC_000962.3/H37Rv/g' | \
		grep -v "duplicata_region" | \
		awk '{if($4 == "insertion"){print($1, $2-250, $3+250, $4)} else {print($0)} }' | sed 's/ /\t/g' | sort -k1,1 -k2,2n > ~/SCRATCH/${sample}.tmp

	bedtools complement -i ~/SCRATCH/${sample}.tmp -g "REF/H37Rv.genome" | awk '{print $0, "NEUTRAL"}' | sed 's/ /\t/g' >> ~/SCRATCH/${sample}.tmp2
	cat ~/SCRATCH/${sample}.tmp ~/SCRATCH/${sample}.tmp2 | sort -k1,1 -k2,2n > ~/SCRATCH/${sample}.bed

	haplotype=$(echo $LINE | cut -f2 -d ' ')

	#python3 get_interval_size.py ~/SCRATCH/${sample}.bed $sample >> region_sizes.tsv

	#get the reference
	bcftools view $reference_vcf --samples $haplotype | grep -v "0:" > ~/SCRATCH/ref1.vcf
	bcftools norm -a -m- ~/SCRATCH/ref1.vcf -Oz -o ~/SCRATCH/${sample}.vcf.gz | grep -v "0:" 2> /dev/null

	bcftools view ~/SCRATCH/${sample}.vcf.gz -H -v indels | awk '{print($1, $2 - 25, $2 + length($5) + 25, "INDELS")}' | sed 's/ /\t/g' > ~/SCRATCH/indels.bed
	bedtools complement -i ~/SCRATCH/indels.bed -g "REF/H37Rv.genome" | awk '{print $0, "NOT_INDELS"}' | sed 's/ /\t/g' >> ~/SCRATCH/indels.bed

	vcf2metrics.py --sample $sample -i MPILEUP_FILTERED_VCF/${sample}.vcf.gz \
		--pipeline "mpileup" \
		--SV "INDELS" \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence \
		--bed ~/SCRATCH/indels.bed

	vcf2metrics.py --sample $sample -i FILTERED_VCF/${sample}.vcf.gz \
		--pipeline "genotube" \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence \
		--bed <(sort -k1,1 -k2,2n ~/SCRATCH/${sample}.bed)

	vcf2metrics.py --sample $sample -i MPILEUP_FILTERED_VCF/${sample}.vcf.gz \
		--pipeline "genotube" \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence \
		--bed <(sort -k1,1 -k2,2n ~/SCRATCH/${sample}.bed)

done < <(grep -v 'strain' maketube_run/my_run_arborescence.tsv | grep -v "NA" )
