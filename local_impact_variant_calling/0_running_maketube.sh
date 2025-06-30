#!/usr/bin/env bash

#generating a maketube run from reference sequence

#bcftools norm -a -m- maketube_run_BLRg/SV1/pop1/pop1.vcf.gz | bcftools view -s H1 -v snps,indels | grep -v "	0:441453" > reference.vcf ; bgzip reference.vcf.gz
rm -rf maketube_run_with_no_homo_seqs

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
