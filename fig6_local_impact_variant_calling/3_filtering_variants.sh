 #!/usr/bin/env bash

mkdir -p BAM FILTERED_VCF

if [[ "$1" == "--filter" ]]
then
	for vcf in VCF/*.vcf.gz
	do
		sample=$(basename $vcf .vcf.gz)
		bcftools index -f $vcf
		bcftools norm -a -m- $vcf -Oz -o ~/SCRATCH/tmp.vcf.gz

		bcftools view ~/SCRATCH/tmp.vcf.gz -i "MAX(FORMAT/AO)/FORMAT/DP>0.8 && MAX(FORMAT/AD)>5 && QUAL>30" -Oz > FILTERED_VCF/${sample}.vcf.gz
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

	cat $strain_regions | grep 300 | cut -f1-4 > ~/SCRATCH/$sample.regions.300.bed
	cat $strain_regions | grep duplicated >> ~/SCRATCH/$sample.regions.300.bed

	bedtools makewindows -b ~/SCRATCH/$sample.regions.300.bed -w 1000 -i src > ~/SCRATCH/${sample}.subsample.bed

	for region_class in $(cut -f4 ~/SCRATCH/${sample}.subsample.bed | sort -u)
	do
		grep $region_class ~/SCRATCH/${sample}.subsample.bed | shuf | head -n10
	done > ~/SCRATCH/${sample}.all_regions_are_10k.bed

	bedtools complement -i ~/SCRATCH/${sample}.300.tmp -g "REF/H37Rv.genome" > ~/SCRATCH/$sample.neutral.bed

	bedtools makewindows -b ~/SCRATCH/$sample.neutral.bed -w 1000 | shuf | head -n10 | \
		awk '{print $0, "NEUTRAL"}' | sed 's/ /\t/g' >> ~/SCRATCH/${sample}.all_regions_are_10k.bed

	haplotype=$(echo $LINE | cut -f2 -d ' ')

	cat $strain_regions | grep "duplicated" > ~/SCRATCH/dupli_region.bed
	fasta_length=$(cat $fasta | grep -E -o "G|C|T|A|N" | wc -l)
	echo | awk -v fasta_lgth=$fasta_length '{print "H37Rv\t"fasta_lgth-150000"\t"fasta_lgth"\tduplicata_region"}' >> ~/SCRATCH/dupli_region.bed

	#get the reference
	bcftools view $reference_vcf --samples $haplotype -v snps | grep -v "0:" > ~/SCRATCH/ref1.vcf
	bcftools norm -a -m- ~/SCRATCH/ref1.vcf -o ~/SCRATCH/${sample}.vcf.gz 2> /dev/null

	bcftools norm -a -m- FILTERED_VCF/${sample}.vcf.gz -Oz | bcftools view -v snps -Oz > ~/SCRATCH/genotube.vcf.gz

	vcf2metrics.py --sample $sample -i ~/SCRATCH/genotube.vcf.gz \
		--add_col "genotube" $sample "with_dupli" \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence \
		--bed ~/SCRATCH/${sample}.all_regions_are_10k.bed

	vcf2metrics.py --sample $sample -i ~/SCRATCH/genotube.vcf.gz \
		--add_col "genotube" $sample "without_dupli" \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence --subtract ~/SCRATCH/dupli_region.bed \
		--bed ~/SCRATCH/${sample}.all_regions_are_10k.bed

done < <(grep -v 'strain' maketube_run/my_run_arborescence.tsv | grep -v "NA" )

exit


<<FANCY_STUFF

	cat $strain_regions | grep "50" | sort -k1,1 -k2,2n > ~/SCRATCH/${sample}.50.tmp
	cat $strain_regions | grep "100" | sort -k1,1 -k2,2n > ~/SCRATCH/${sample}.100.tmp
	cat $strain_regions | grep "150" | sort -k1,1 -k2,2n > ~/SCRATCH/${sample}.150.tmp
	cat $strain_regions | grep "200" | sort -k1,1 -k2,2n > ~/SCRATCH/${sample}.200.tmp
	cat $strain_regions | grep "250" | sort -k1,1 -k2,2n > ~/SCRATCH/${sample}.250.tmp
	vcf2metrics.py --sample $sample -i FILTERED_VCF/${sample}.vcf.gz \
		--add_col "genotube" $sample 50 \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence \
		--bed ~/SCRATCH/${sample}.50.tmp

	vcf2metrics.py --sample $sample -i FILTERED_VCF/${sample}.vcf.gz \
		--add_col "genotube" $sample 100 \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence \
		--bed ~/SCRATCH/${sample}.100.tmp

	vcf2metrics.py --sample $sample -i FILTERED_VCF/${sample}.vcf.gz \
		--add_col "genotube" $sample 150 \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence \
		--bed ~/SCRATCH/${sample}.150.tmp

	vcf2metrics.py --sample $sample -i FILTERED_VCF/${sample}.vcf.gz \
		--add_col "genotube" $sample 200 \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence \
		--bed ~/SCRATCH/${sample}.200.tmp

	vcf2metrics.py --sample $sample -i FILTERED_VCF/${sample}.vcf.gz \
		--add_col "genotube" $sample 250 \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence \
		--bed ~/SCRATCH/${sample}.250.tmp

	vcf2metrics.py --sample $sample -i FILTERED_VCF/${sample}.vcf.gz \
		--add_col "genotube" $sample 300 \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence \
		--bed ~/SCRATCH/${sample}.300.tmp


	vcf2metrics.py --sample $sample -i MPILEUP_FILTERED_VCF/${sample}.vcf.gz \
		--add_col "TBprofiler" $sample 50 \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence \
		--bed ~/SCRATCH/${sample}.50.tmp

	vcf2metrics.py --sample $sample -i MPILEUP_FILTERED_VCF/${sample}.vcf.gz \
		--add_col "TBprofiler" $sample 100 \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence \
		--bed ~/SCRATCH/${sample}.100.tmp

	vcf2metrics.py --sample $sample -i MPILEUP_FILTERED_VCF/${sample}.vcf.gz \
		--add_col "TBprofiler" $sample 150 \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence \
		--bed ~/SCRATCH/${sample}.150.tmp

	vcf2metrics.py --sample $sample -i MPILEUP_FILTERED_VCF/${sample}.vcf.gz \
		--add_col "TBprofiler" $sample 200 \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence \
		--bed ~/SCRATCH/${sample}.200.tmp

	vcf2metrics.py --sample $sample -i MPILEUP_FILTERED_VCF/${sample}.vcf.gz \
		--add_col "TBprofiler" $sample 250 \
		--reference ~/SCRATCH/${sample}.vcf.gz \
		--backtrack $equivalence \
		--bed ~/SCRATCH/${sample}.250.tmp
FANCY_STUFF

#bcftools view ~/SCRATCH/${sample}.vcf.gz -H -v indels | awk '{print($1, $2 - 25, $2 + length($5) + 25, "INDELS")}' | sed 's/ /\t/g' > ~/SCRATCH/indels.bed
#bedtools complement -i ~/SCRATCH/indels.bed -g "REF/H37Rv.genome" | awk '{print $0, "NOT_INDELS"}' | sed 's/ /\t/g' >> ~/SCRATCH/indels.bed

#vcf2metrics.py --sample $sample -i MPILEUP_FILTERED_VCF/${sample}.vcf.gz \
#	--pipeline "mpileup" \
#	--SV "INDELS" \
#	--reference ~/SCRATCH/${sample}.vcf.gz \
#	--backtrack $equivalence \
#	--bed ~/SCRATCH/indels.bed
