#!/usr/bin/env bash

#can be skipped altogether, all necessary files are precomputed
if [ "$1" == "-r" ]; then
	echo "rerun the complete analysis"

	#should be run once, the blast db file was too big
	cat REAL_STRAINS/*.fasta > BLAST_DB/genome_concatenation.fasta
	makeblastdb -in BLAST_DB/genome_concatenation.fasta -title all_tuberculosis_genomes -dbtype nucl -out BLAST_DB/genome_concatenation
	blastn -query BLAST_DB/IS6110.fasta -db BLAST_DB/genome_concatenation -outfmt "6 qseqid sseqid pident sstart send" > BLAST_DB/IS6110_query_results.tsv

	#parse reference line by line from natural_strains_infos
	#for each genome generate a run of 
	while read reference_line
	do
		reference=$(echo $reference_line | cut -f1 -d ' ')
		reference_file=$(echo $reference_line | cut -f3 -d ' ')
		reference_transposon=$(echo $reference_line | cut -f4 -d ' ')
		reference_lineage=$(echo $reference_line | cut -f5 -d ' ')
		
		#maketube genome are manually copy pasted
		Rscript maketube.R \
			--reference $reference_file --transposon $reference_transposon --nonhomoseq "REF/nonH37Rv_pool_sequence_cp.fasta" \
			--output "maketube_run/${reference}_run" \
			--structural_variants 3 --haplotype_count 10 --pop_count 1 \
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
fi

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


		#all done accordingly to nucmer documentation
		if [ ! -f "DELTAS/${sample}_${reference}.delta" ]
		then
			nucmer -p "DELTAS/${sample}_${reference}" $reference_file $sample_file 2> /dev/null
			delta-filter -1 "DELTAS/${sample}_${reference}.delta" > "DELTAS/${sample}_${reference}.filter"
		fi

		if [ ! -f "DELTAS/${sample}_${reference}.report" ]
		then
			dnadiff -d "DELTAS/${sample}_${reference}.delta" -p "DELTAS/${sample}_${reference}" 2> /dev/null
		fi

		#grep aligned bases and total bases to avoid the nucmer rounding
		ref_length=$(sed -E 's/\s+/\t/g' "DELTAS/${sample}_${reference}.report" | grep "TotalBases" | cut -f2)
		seq_length=$(sed -E 's/\s+/\t/g' "DELTAS/${sample}_${reference}.report" | grep "TotalBases" | cut -f3)

		#clean up the report line to extract only the value
		aligned_ref=$(head -n20 "DELTAS/${sample}_${reference}.report" | sed -E 's/\s+/\t/g' | grep "AlignedBases" | cut -f2 | cut -f1 -d '(')
		aligned_sample=$(head -n20 "DELTAS/${sample}_${reference}.report" | sed -E 's/\s+/\t/g' | grep "AlignedBases" | cut -f3 | cut -f1 -d '(')

		echo -e $reference $reference_lineage $reference_source $sample $sample_lineage  $sample_source $ref_length $aligned_ref $seq_length $aligned_sample | sed 's/ /\t/g'

	done < <(cat natural_strains_info.txt | grep -v '#' | grep "sample")
done < <(cat natural_strains_info.txt | grep -v '#' | grep "reference_genome")

