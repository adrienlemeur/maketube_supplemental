#!/usr/bin/env bash

all_strains=$(find ./ -type f \( -iname \*.fasta -o -iname \*.fa \))

mkdir -p VCF
mkdir -p F/F1 F/F2 F/F3 PSEUDOVCF
mkdir -p /home/alemeur/SCRATCH/


rm -rf BAM_30/*.nopg.bam GATK_VCF
mkdir -p GATK_VCF

ls BAM_30/*.bam | \
	parallel -j 5 '''export sample=$(basename {} .bam) ; echo $sample ; \
	samtools addreplacerg BAM_30/${sample}.bam \
		-r "@RG\tID:$sample\tLB:LIB\tPL:Illumina\tPU:UNIT\tSM:$sample" --output-fmt "BAM" > BAM_30/${sample}.nopg.bam;
	samtools index BAM_30/${sample}.nopg.bam
	'''

ls BAM_30/*.nopg.bam | \
	parallel -j 5 '''export sample=$(basename {} .nopg.bam) ; echo $sample ; \
	singularity exec ~/CONTAINERS/gatk.img gatk HaplotypeCaller \
		--reference REF/H37Rv.fasta --input BAM_30/${sample}.nopg.bam \
		--sample-ploidy 2 --native-pair-hmm-threads 2 \
		--annotation StrandBiasBySample \
		--output GATK_VCF/$sample.vcf'''

for i in GATK_VCF/*.vcf ; do bgzip $i ; done

ls GATK_VCF/*.vcf.gz | grep -v MTBseq | \
	parallel -j 5 '''export sample=$(basename {} .vcf.gz) ; echo $sample ; \
	bcftools view GATK_VCF/$sample.vcf.gz -i "FORMAT/AD[0:1]/FORMAT/DP >= 0.75 && SB[0:2]>4 && SB[0:3]>4" \
	-Oz > GATK_VCF/$sample.MTBseq.vcf.gz
	'''
