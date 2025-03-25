#!/usr/bin/env bash

<<ONE_TIME

cat REAL_STRAINS/*.fasta > BLAST_DB/genome_concatenation.fasta
makeblastdb -in BLAST_DB/genome_concatenation.fasta -title all_tuberculosis_genomes -dbtype nucl -out BLAST_DB/genome_concatenation
blastn -query BLAST_DB/IS6110.fasta -db BLAST_DB/genome_concatenation -outfmt "6 qseqid sseqid pident sstart send" > BLAST_DB/IS6110_query_results.tsv

while read reference_line
do
	reference=$(echo $reference_line | cut -f1 -d ' ')
	reference_file=$(echo $reference_line | cut -f3 -d ' ')
	reference_transposon=$(echo $reference_line | cut -f4 -d ' ')
	reference_lineage=$(echo $reference_line | cut -f5 -d ' ')

	Rscript maketube.R \
		--reference $reference_file --transposon $reference_transposon --nonhomoseq "REF/nonH37Rv_pool_sequence_cp.fasta" \
		--output "maketube_run/${reference}_run" \
		--structural_variants 5 --haplotype_count 10 --pop_count 1 \
		--unmuted \
		--mutation_rate 1.23e-7 --scaling_factor 0.125 \
		--effective_pop 700 \
		--TCAGI "0.172,0.329,0.172,0.329,0" \
		--ABCDEF "0.65,0.05,0.21,0.27,0.02,0.65" \
		--slope "50,100,150,200,250,300"

	snpmutator "$reference_file" \
		--num-simulations 15 \
		--num-substitutions 800 \
		--num-insertions 50 \
		--num-deletions 50 \
		--mono \
		--seqid "${reference}" \
		--fasta-dir "snpmutator_run/${reference}"

done < <(grep "reference_genome" natural_strains_info.txt | grep -v "#")
ONE_TIME


mkdir -p DELTAS BAM
while read reference_line
do
	reference=$(echo $reference_line | cut -f1 -d ' ')
	reference_file=$(echo $reference_line | cut -f3 -d ' ')
	reference_lineage=$(echo $reference_line | cut -f5 -d ' ' | sed -E 's/.$//')

	while read sample_line
	do
		sample=$(echo $sample_line | cut -f1 -d ' ')
		sample_file=$(echo $sample_line | cut -f3 -d ' ')
		sample_lineage=$(echo $sample_line | cut -f5 -d ' ' | sed -E 's/.$//')

		if [ ! -f "BAM/${sample}_${reference}.bam" ]
		then
			minimap2 -ax asm5 $sample_file $reference_file | samtools sort -@ 5 --no-PG -o "BAM/${sample}_${reference}.bam" 2> /dev/null
		fi

		if [ ! -f "BAM/${reference}_${sample}.bam" ]
		then
			minimap2 -ax asm5 $reference_file $sample_file | samtools sort -@ 5 --no-PG -o "BAM/${reference}_${sample}.bam" 2> /dev/null
		fi

		sample_to_reference=$(samtools coverage "BAM/${reference}_${sample}.bam" | tail -n1 | cut -f3,5,6)
		reference_to_sample=$(samtools coverage "BAM/${sample}_${reference}.bam" | tail -n1 | cut -f3,5,6)
	
		echo -e $reference $reference_lineage $sample $sample_lineage $reference_to_sample $sample_to_reference | sed 's/ /\t/g'

	done < <(cat natural_strains_info.txt | grep "sample" | grep -v "#")
done < <(cat natural_strains_info.txt | grep "reference_genome" | grep -v "#")

exit



