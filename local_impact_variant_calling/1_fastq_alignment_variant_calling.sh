#!/usr/bin/env bash

mkdir -p BAM

#creating fastq archives and aligning them on the reference genome
while read LINE
do

	sample=$(echo $LINE | cut -f1 -d ' ')
	echo Doing... $sample

	fasta_gz=$(echo $LINE | cut -f9 -d ' ')

	fasta=$(echo $fasta_gz | sed 's/\.gz//g')

	R1=$(echo $LINE | cut -f10 -d ' ')
	R2=$(echo $LINE | cut -f11 -d ' ')
	fastq=$(echo $R1 | sed 's/1.fq.gz//g')

	rm -f $R1 $R2

	gunzip --keep $fasta_gz 2> /dev/null

	art_illumina -i $fasta -na \
		-o $fastq -ss "MSv3" -l '250' \
		-f '30' -p -m '500' -s '30' -rs '167631' \
	&& bgzip "${fastq}1.fq" && bgzip "${fastq}2.fq"

	bwa-mem2 mem "REF/H37Rv.fasta" $R1 $R2 -t 5 2>> /dev/null | \
		samtools view -@ 5 -bSh --no-PG 2>> /dev/null | \
		samtools sort -@ 5 --no-PG -o BAM/${sample}.bam 2> /dev/null
	#-q 30
done < <(grep -v 'strain' maketube_run/my_run_arborescence.tsv)
