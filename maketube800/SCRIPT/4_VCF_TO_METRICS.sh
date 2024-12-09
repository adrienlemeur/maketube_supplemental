#!/usr/bin/env bash

all_strains=$(find ./ -type f \( -iname \*.fasta -o -iname \*.fa \))
REF="REF/H37Rv.fasta"

echo -e "pipeline\tsample\tfilter\tvariant_caller\tregion\tregion2\tvariant_class\tTP\tFP\tFN\tRECALL\tPRECISION\tF1"

while read LINE
do
	sample=$(echo $LINE | cut -f1 -d ' ')
	source=$(echo $LINE | cut -f2 -d ' ')
	reference=$(echo $LINE | cut -f3 -d ' ')
	EQUIVALENCE=$(echo $LINE | cut -f4 -d ' ')
	BED=$(echo $LINE | cut -f5 -d ' ')
	PSEUDOVCF=$(echo $LINE | cut -f6 -d ' ')
	PAF_VCF=$(echo $LINE | cut -f7 -d ' ')

	bcftools norm -a -m- $reference -Oz -o ~/SCRATCH/reference.vcf.gz 2>> /dev/null

	bcftools norm -a -m- $PAF_VCF -Oz -o ~/SCRATCH/PAF_VCF.vcf.gz 2>> /dev/null
	bcftools norm -a -m- $PSEUDOVCF -Oz -o ~/SCRATCH/PSEUDOVCF.vcf.gz 2>> /dev/null

	vcf2metrics.py --sample $sample -i ~/SCRATCH/PSEUDOVCF.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --backtrack $EQUIVALENCE --pipeline $source --SV "nucmer"
	vcf2metrics.py --sample $sample -i ~/SCRATCH/PAF_VCF.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --backtrack $EQUIVALENCE --pipeline $source --SV "PAF"

done < <(grep -v '#' a_database_to_rule_them_all.tsv | grep "maketube_strains" | grep -e "H1\." -e "H10\." )

while read LINE
do
	sample=$(echo $LINE | cut -f1 -d ' ')
	source=$(echo $LINE | cut -f2 -d ' ')
	reference=$(echo $LINE | cut -f3 -d ' ')
	EQUIVALENCE=$(echo $LINE | cut -f4 -d ' ')
	BED=$(echo $LINE | cut -f5 -d ' ')
	PSEUDOVCF=$(echo $LINE | cut -f6 -d ' ')
	PAF_VCF=$(echo $LINE | cut -f7 -d ' ')

	bcftools norm -a -m- $reference -Oz -o ~/SCRATCH/reference.vcf.gz 2>> /dev/null

	bcftools norm -a -m- $PAF_VCF -Oz -o ~/SCRATCH/PAF_VCF.vcf.gz 2>> /dev/null
	bcftools norm -a -m- $PSEUDOVCF -Oz -o ~/SCRATCH/PSEUDOVCF.vcf.gz 2>> /dev/null

	vcf2metrics.py --sample $sample -i ~/SCRATCH/PSEUDOVCF.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --pipeline $source --SV "nucmer"
	vcf2metrics.py --sample $sample -i ~/SCRATCH/PAF_VCF.vcf.gz --reference ~/SCRATCH/reference.vcf.gz --pipeline $source --SV "PAF"

done < <(grep -v '#' a_database_to_rule_them_all.tsv | grep "snpmutator" | grep -e "H1\." -e "H10\." )


while read LINE
do
	sample=$(echo $LINE | cut -f1 -d ' ')
	source=$(echo $LINE | cut -f2 -d ' ')
	reference=$(echo $LINE | cut -f3 -d ' ')

	MUMMER_VCF=$(echo $LINE | cut -f6 -d ' ')
	PAF_VCF=$(echo $LINE | cut -f7 -d ' ')

	RAW_FREEBAYES=$(echo $LINE | cut -f8 -d ' ')
	GENOTUBE=$(echo $LINE | cut -f11 -d ' ')

	RAW_GATK=$(echo $LINE | cut -f12 -d ' ')
	MTBseq=$(echo $LINE | cut -f13 -d ' ')

	MPILEUP=$(echo -n $LINE | cut -f14 -d ' ' | sed -E 's/.$//')

	bcftools norm -a $reference -Oz > ~/SCRATCH/reference.vcf.gz 2>> /dev/null

	bcftools norm -a $MUMMER_VCF -Oz > ~/SCRATCH/mummer.vcf.gz 2>> /dev/null
	bcftools norm -a $PAF_VCF -Oz > ~/SCRATCH/paf.vcf.gz 2>> /dev/null

	bcftools norm -a $RAW_FREEBAYES -Oz > ~/SCRATCH/raw_freebayes.vcf.gz 2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/raw_freebayes.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--depth "freebayes" --sample $sample --pipeline $source --SV "RAW"

	bcftools norm -a $GENOTUBE -Oz > ~/SCRATCH/genotube_freebayes.vcf.gz 2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/genotube_freebayes.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--depth "freebayes" --sample $sample --pipeline $source --SV "genotube"

	bcftools norm -a $MPILEUP | bcftools view -v snps,mnps,indels -Oz > ~/SCRATCH/tbprofiler.vcf.gz 2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/tbprofiler.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--depth "mpileup" --sample $sample --pipeline $source --SV "TBprofiler"

	bcftools norm -a $RAW_GATK -Oz > ~/SCRATCH/raw_gatk.vcf.gz 2>> /dev/null	
	vcf2metrics.py -i ~/SCRATCH/raw_gatk.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--depth "GATK" --sample $sample --pipeline $source --SV "RAW"

	bcftools norm -a $MTBseq -Oz > ~/SCRATCH/mtbseq.vcf.gz 2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/mtbseq.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--depth "GATK" --sample $sample --pipeline $source --SV "MTBseq"

done < <(grep -v '#' a_database_to_rule_them_all.tsv | grep "snpmutator" | grep -e "H1\." -e "H10\.")

while read LINE
do
	sample=$(echo $LINE | cut -f1 -d ' ')
	source=$(echo $LINE | cut -f2 -d ' ')
	reference=$(echo $LINE | cut -f3 -d ' ')

	EQUIVALENCE=$(echo $LINE | cut -f4 -d ' ')
	REGIONS_OF_INTEREST=$(echo $LINE | cut -f5 -d ' ')

	MUMMER_VCF=$(echo $LINE | cut -f6 -d ' ')
	PAF_VCF=$(echo $LINE | cut -f7 -d ' ')

	RAW_FREEBAYES=$(echo $LINE | cut -f8 -d ' ')
	GENOTUBE=$(echo $LINE | cut -f11 -d ' ')

	RAW_GATK=$(echo $LINE | cut -f12 -d ' ')
	MTBseq=$(echo $LINE | cut -f13 -d ' ')

	MPILEUP=$(echo $LINE | cut -f14 -d ' ' | sed -E 's/.$//')

	bcftools norm -a $reference -Oz > ~/SCRATCH/reference.vcf.gz 2>> /dev/null

	bcftools norm -a $MUMMER_VCF -Oz > ~/SCRATCH/mummer.vcf.gz 2>> /dev/null
	bcftools norm -a $PAF_VCF -Oz > ~/SCRATCH/paf.vcf.gz 2>> /dev/null

	bcftools norm -a $RAW_FREEBAYES -Oz > ~/SCRATCH/raw_freebayes.vcf.gz 2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/raw_freebayes.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--backtrack $EQUIVALENCE --depth "freebayes" \
		--sample $sample --pipeline $source --SV "RAW"

	bcftools norm -a $GENOTUBE -Oz > ~/SCRATCH/genotube_freebayes.vcf.gz 2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/genotube_freebayes.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--backtrack $EQUIVALENCE --depth "freebayes" \
		--sample $sample --pipeline $source --SV "genotube"

	bcftools norm -a $MPILEUP | bcftools view -v snps,mnps,indels -Oz > ~/SCRATCH/tbprofiler.vcf.gz 2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/tbprofiler.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--backtrack $EQUIVALENCE --depth "mpileup" \
		--sample $sample --pipeline $source --SV "TBprofiler"

	bcftools norm -a $RAW_GATK -Oz > ~/SCRATCH/raw_gatk.vcf.gz 2>> /dev/null	
	vcf2metrics.py -i ~/SCRATCH/raw_gatk.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--backtrack $EQUIVALENCE --depth "GATK" \
		--sample $sample --pipeline $source --SV "RAW"

	bcftools norm -a $MTBseq -Oz > ~/SCRATCH/mtbseq.vcf.gz 2>> /dev/null
	vcf2metrics.py -i ~/SCRATCH/mtbseq.vcf.gz --reference ~/SCRATCH/reference.vcf.gz \
		--backtrack $EQUIVALENCE --depth "GATK" \
		--sample $sample --pipeline $source --SV "MTBseq"
done < <(grep -v '#' a_database_to_rule_them_all.tsv | grep "maketube_strains" | grep -e "H1\." -e "H10\.")

