#!/usr/bin/env bash

#you may skip this altogether, all necessary files are precomputed
echo "rerun the complete analysis"

#should be run once, the blast db file was too big
#cat REAL_STRAINS/*.fasta > BLAST_DB/genome_concatenation.fasta
#makeblastdb -in BLAST_DB/genome_concatenation.fasta -title all_tuberculosis_genomes -dbtype nucl -out BLAST_DB/genome_concatenation
#blastn -query BLAST_DB/IS6110.fasta -db BLAST_DB/genome_concatenation -outfmt "6 qseqid sseqid pident sstart send" > BLAST_DB/IS6110_query_results.tsv

#parse reference line by line from natural_strains_infos
#for each genome generate a run of maketube and snp mutator
while read reference_line
do
	reference=$(echo $reference_line | cut -f1 -d ' ')
	reference_file=$(echo $reference_line | cut -f3 -d ' ')
	reference_transposon=$(echo $reference_line | cut -f4 -d ' ')
	reference_lineage=$(echo $reference_line | cut -f5 -d ' ')

	#maketube genome are THEN manually copy pasted into artificial strains
	sudo Rscript ./maketube.R \
		--reference $reference_file --transposon $reference_transposon --nonhomoseq "REF/nonH37Rv_pool_sequence_cp.fasta" \
		--output "maketube_run/${reference}_run" \
		--structural_variants 3 --haplotype_count 10 --pop_count 1 \
		--unmuted \
		--mutation_rate 1.23e-7 --scaling_factor 0.125 \
		--effective_pop 700 \
		--TCAGI "0.172,0.329,0.172,0.329,0" \
		--ABCDEF "0.65,0.05,0.21,0.27,0.02,0.65" \
		--slope "50,100,150,200,250,300"

done < <(grep "reference_genome" natural_strains_info.txt | grep -v "#")

while read reference_line
do
	reference=$(echo $reference_line | cut -f1 -d ' ')
	reference_file=$(echo $reference_line | cut -f3 -d ' ')
	reference_transposon=$(echo $reference_line | cut -f4 -d ' ')
	reference_lineage=$(echo $reference_line | cut -f5 -d ' ')

	mkdir -p "snpmutator_run/$reference/tmp"

	snpmutator "$reference_file" \
		--num-simulations 30 \
		--num-substitutions 800 \
		--num-insertions 50 \
		--num-deletions 50 \
		--vcf ${reference} \
		--mono \
		--seqid "${reference}" \
		--fasta-dir "snpmutator_run/${reference}" \
		--vcf "snpmutator_run/${reference}/tmp/tmp.vcf"

	bgzip -f "snpmutator_run/${reference}/tmp/tmp.vcf"
	bcftools index "snpmutator_run/${reference}/tmp/tmp.vcf.gz"

	for ((i=1;i<=30;i++))
	do
		bcftools view -s "${reference}_mutated_$i" "snpmutator_run/${reference}/tmp/tmp.vcf.gz" | grep -v "GT	0" | bcftools view -Oz -o "snpmutator_run/${reference}//H37Rv_mutated_$i.vcf.gz"
	done

done < <(grep "reference_genome" natural_strains_info.txt | grep -v "#")
