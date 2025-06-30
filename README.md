This is the repository for the maketube paper complete data.
It includes : 
	- the sequences used (with the exception of larger files)
	- the script used for the analysis
	- the figures & the script used to generate the figures.

Exhaustive description :

├── genome_diversity_of_maketube_genomes			# analysis of sequence diversity between natural and artificial genome
│   ├── REF
│   ├── DELTAS							# nucmer intermediary files .delta, .snps, etc.
│   ├── artificial_strains					# snpmutator and maketube genomes
│   ├── REAL_STRAINS						# natural strains sequences (& IS6110 locations, unused)
│   ├── maketube_run						# maketube files : sequence, evolution partitions, annotations files. 
|   │   ├── AF2122
│   │   ├── 18b
│   │   └── H37Rv
│   └── maketube_with_natural_genomes.sh
└────── natural_strains_info.txt				#ID, lineage and additional information about natural strains

