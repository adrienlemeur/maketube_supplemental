#!/usr/bin/env bash


all_strains=$(find ./ -type f \( -iname \*.fasta -o -iname \*.fa \) | grep -v maketube_strains)
all_strains_maketube=$(find ./ -type f \( -iname \*.fasta -o -iname \*.fa \) | grep pop1 | grep -e "H1\." -e "H10\." )
rm -rf DELTAS
mkdir -p DELTAS
mkdir -p ~/SCRATCH

REF="REF/H37Rv.fasta"
echo -e "sample\treference\taverage_identity\tsample_to_ref\tref_to_sample\tnucmer_snp_count\tnucmer_indel_count" > all_strains_dnadiff.txt

for one_strain in $all_strains $all_strains_maketube
do
	sample=$(basename $one_strain ".fasta" | sed 's/\.fa//g')
	#>&2 echo "doing $one_strain"

	if [ ! -f "DELTAS/$sample.delta" ]
	then
		nucmer -p "DELTAS/"$sample $REF $one_strain
	fi
 
	dnadiff -d "DELTAS/$sample.delta" -p "/home/alemeur/SCRATCH/tmp"

	average_identity=$(cat ~/SCRATCH/tmp.report | grep "AvgIdentity" | head -n1 | sed -E 's/\s+/\t/g' | cut -f2)
	reference_to_sample=$(cat ~/SCRATCH/tmp.report | grep "AlignedBases" | sed -E 's/[(|)]/\t/g' | cut -f2,4 | tr -d '%\n' | cut -f1)
	sample_to_reference=$(cat ~/SCRATCH/tmp.report | grep "AlignedBases" | sed -E 's/[(|)]/\t/g' | cut -f2,4 | tr -d '%\n' | cut -f2)
	nucmer_snp_count=$(cat ~/SCRATCH/tmp.report | grep "TotalSNPs" | head -n1 | sed -E 's/\s+/\t/g' | cut -f2)
	nucmer_indel_count=$(cat ~/SCRATCH/tmp.report | grep "TotalIndels" | head -n1 | sed -E 's/\s+/\t/g' | cut -f2)
	echo -e "$sample\tH37Rv\t$average_identity\t$reference_to_sample\t$sample_to_reference\t$nucmer_snp_count\t$nucmer_indel_count"

done | sort -h -k1 >> all_strains_dnadiff.txt
