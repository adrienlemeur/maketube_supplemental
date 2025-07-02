#!/usr/bin/env bash


all_strains=$(find ./ -type f \( -iname \*.fasta -o -iname \*.fa \) | grep -v maketube_strains)
all_strains_maketube=$(find ./ -type f \( -iname \*.fasta -o -iname \*.fa \) | grep pop1 | grep -e "H1\." -e "H10\." )

rm -f coverage_distribution.txt

rm -rf BAM_30 FASTQ_30
mkdir -p BAM_30 FASTQ_30
REF="REF/H37Rv.fasta"

for one_strain in $all_strains $all_strains_maketube
do
	sample=$(basename $one_strain ".fasta" | sed 's/\.fa//g')

	>&2 echo "doing $one_strain"

	if [ ! -f "FASTQ_30/${sample}_1.fq.gz" ]
	then
		>&2 echo -e "\tgenesis"
		art_illumina -i $one_strain -na \
			-o FASTQ_30/"${sample}_" -ss "MSv3" -l '250' \
			-ir 0 -ir2 0 -dr 0 -dr2 0 -k 0 \
			-f '30' -p -m '500' -s '30' -rs '167631' > /dev/null 2> /dev/null

		>&2 echo -e "\tcompression"
		bgzip -f -l 0 "FASTQ_30/${sample}_1.fq"
		bgzip -f -l 0 "FASTQ_30/${sample}_2.fq"
	fi
done
