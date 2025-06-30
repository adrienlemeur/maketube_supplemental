rm -rf 'output'

mkdir -p 'output'
mkdir -p 'output/info'
mkdir -p 'output/tmp'

END=20
for ((i=1;i<=END;i++))
do
snpmutator "../REF/H37Rv.fasta" \
	--num-simulations 2 \
	--num-substitutions 800 \
	--num-insertions 50 \
	--num-deletions 50 \
	--mono \
	--seqid "snpmutator_SV${i}" \
	--fasta-dir "output/classic_benchmark_fasta_$i" \
	--vcf "output/tmp/tmp_$i.vcf" \
	--metrics "output/info/classic_metrics_$i.txt"

	mv "output/classic_benchmark_fasta_$i/H37Rv_mutated_1.fasta" "output/snpmutator_SV${i}_H1.fasta"
	mv "output/classic_benchmark_fasta_$i/H37Rv_mutated_2.fasta" "output/snpmutator_SV${i}_H2.fasta"

	bgzip "output/tmp/tmp_$i.vcf"
	bcftools index "output/tmp/tmp_$i.vcf.gz"

	bcftools view -s H37Rv_mutated_1 "output/tmp/tmp_$i.vcf.gz" | grep -v "GT	0" | bcftools view -Oz -o "output/snpmutator_SV${i}_H1.vcf.gz"
	bcftools view -s H37Rv_mutated_2 "output/tmp/tmp_$i.vcf.gz" | grep -v "GT	0" | bcftools view -Oz -o "output/snpmutator_SV${i}_H2.vcf.gz"
done
