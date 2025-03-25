		if [ ! -f "DELTAS/${sample}_${reference}.delta" ]
		then
			nucmer -p "DELTAS/${sample}_${reference}" $reference_file $sample_file 2> /dev/null
			delta-filter -1 "DELTAS/${sample}_${reference}.delta" > "DELTAS/${sample}_${reference}.filter"

		if [ ! -f "DELTAS/${sample}_${reference}.delta" ]
		then
			dnadiff -d "DELTAS/${sample}_${reference}.delta" -p "DELTAS/${sample}_${reference}" 2> /dev/null
		fi

		#ref_length=$(sed -E 's/\s+/\t/g' "DELTAS/${sample}_${reference}.report" | grep "TotalBases" | cut -f2)
		#seq_length=$(sed -E 's/\s+/\t/g' "DELTAS/${sample}_${reference}.report" | grep "TotalBases" | cut -f3)

		#aligned_ref=$(head -n20 "DELTAS/${sample}_${reference}.report" | sed -E 's/\s+/\t/g' | grep "TotalLength" | cut -f2)
		#aligned_sample=$(head -n20 "DELTAS/${sample}_${reference}.report" | sed -E 's/\s+/\t/g' | grep "TotalLength" | cut -f3)

