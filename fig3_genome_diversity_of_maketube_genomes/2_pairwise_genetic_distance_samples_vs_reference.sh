#!/usr/bin/env bash

mkdir -p DELTAS BAM

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

		#all done accordingly to nucmer documentation
		if [ ! -f "DELTAS/${unique_sample}_${reference}.delta" ]
		then
			nucmer -p "DELTAS/${unique_sample}_${reference}" $reference_file $sample_file 2> /dev/null
			delta-filter -1 "DELTAS/${unique_sample}_${reference}.delta" > "DELTAS/${unique_sample}_${reference}.filter"
		fi

		if [ ! -f "DELTAS/${unique_sample}_${reference}.report" ]
		then
			dnadiff -d "DELTAS/${unique_sample}_${reference}.delta" -p "DELTAS/${unique_sample}_${reference}" 2> /dev/null
		fi

		#grep aligned bases and total bases to avoid the nucmer rounding
		ref_length=$(sed -E 's/\s+/\t/g' "DELTAS/${unique_sample}_${reference}.report" | grep "TotalBases" | cut -f2)
		seq_length=$(sed -E 's/\s+/\t/g' "DELTAS/${unique_sample}_${reference}.report" | grep "TotalBases" | cut -f3)

		#clean up the report line to extract only the value
		aligned_ref=$(head -n20 "DELTAS/${unique_sample}_${reference}.report" | sed -E 's/\s+/\t/g' | grep "AlignedBases" | cut -f2 | cut -f1 -d '(')
		aligned_sample=$(head -n20 "DELTAS/${unique_sample}_${reference}.report" | sed -E 's/\s+/\t/g' | grep "AlignedBases" | cut -f3 | cut -f1 -d '(')

		echo -e $reference $reference_lineage $reference_source $unique_sample $sample_lineage  $sample_source $ref_length $aligned_ref $seq_length $aligned_sample | sed 's/ /\t/g'

	done < <(cat natural_strains_info.txt | grep -v '#' | grep "sample")
done < <(cat natural_strains_info.txt | grep -v '#' | grep "reference_genome")

