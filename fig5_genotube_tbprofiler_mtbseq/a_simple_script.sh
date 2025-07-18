#!/usr/bin/env bash

while read LINE
do
	sample=$(echo $LINE | cut -f1 -d ' ')
	strain=$(echo $LINE | cut -f2 -d ' ')

	reference_vcf=$(echo $LINE | cut -f11 -d ' ')
	maketube_equivalence=$(echo $LINE | cut -f4 -d ' ')

	cat $maketube_equivalence | grep deletion | awk '{print $3-$2}' #| awk 'BEGIN {FS = "|"} ; {sum+=$1} END {print sum}'

done < <(grep -v '#' a_database_to_rule_them_all.tsv | grep "maketube")
