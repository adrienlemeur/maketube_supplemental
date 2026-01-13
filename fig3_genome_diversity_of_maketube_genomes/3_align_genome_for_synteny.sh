#!/usr/bin/env bash

mkdir -p PAF

#compare all genome to the 3 reference genome annotated reference in natural_strains_info
while read reference_line
do
	reference=$(echo $reference_line | cut -f1 -d ' ')
	reference_file=$(echo $reference_line | cut -f3 -d ' ')
	reference_lineage=$(echo $reference_line | cut -f5 -d ' ')
	reference_source=$(echo $reference_line | cut -f6 -d ' ' | sed -E 's/.$//')

	while read sample_line
	do

		sample=$(echo $sample_line | cut -f1 -d ' ')
		sample_file=$(echo $sample_line | cut -f3 -d ' ')
		sample_lineage=$(echo $sample_line | cut -f5 -d ' ')
		sample_source=$(echo $sample_line | cut -f6 -d ' ' | sed -E 's/.$//')

		unique_sample=${sample}_${sample_lineage}

		#echo "Doing $unique_sample..."
		echo "reference $reference sample $sample"
		minimap2 $reference_file -x asm5 $sample_file -r 150,1000 --secondary=no > PAF/${unique_sample}_${reference}.paf

	done < <(cat natural_strains_info.txt | grep -v '#' | grep "sample")
done < <(cat natural_strains_info.txt | grep -v '#' | grep "reference_genome")

