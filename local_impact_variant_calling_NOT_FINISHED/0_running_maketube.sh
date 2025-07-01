#!/usr/bin/env bash

#generating a maketube run from reference sequence
Rscript maketube.R \
	--reference "REF/H37Rv.fasta" --transposon "REF/H37Rv_transposon.bed" --nonhomoseq "REF/nonH37Rv_pool_sequence_cp.fasta" \
	--output "maketube_run" \
	--structural_variants 20 --haplotype_count 10 --pop_count 1 \
	--unmuted \
	--mutation_rate 1.23e-7 --scaling_factor 0.125 \
	--effective_pop 2500 \
	--TCAGI "0.172,0.329,0.172,0.329,0" \
	--ABCDEF "0.65,0.05,0.21,0.27,0.02,0.65" \
	--slope "50,100,150,200,250,300"
