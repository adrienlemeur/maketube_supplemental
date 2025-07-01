<h1 align="center"> Maketube supplemental </h1>

This is the repository for the maketube article (link).

### Table of contents

<!--ts-->
-  [Figure 3 : ](#genome_diversity_tag)


### <a name="genome_diversity_tag"></a>Diversity of artificial and natural genomes

Dependencies :
Tools : snpmutator, maketube (included)
R libs : 
	maketube : seqinr, jackalope, optparse, Biostrings, dplyr, ape
	visualisation : ggplot2, scales, ggrepel, ggforce

Descriptions :

Maketube and snpmutator are first launched to generate strains from the 3 reference genomes (natural_strains_info.txt).
A subsample of 10 strains for each reference, was manually copied into the "artificial_strains" folder.
These genomes, and natural genomes are aligned onto the three reference genomes using nucmer, and the resulting delta files are filtered. A report is generated with dnadiff.
The number of base of the considered genome aligned on the reference, as well as the number of base of the reference aligned on the considered genome, are grepped and conserved in distance_to_reference.tsv.
The visualisation is made from distance_to_reference.tsv.

File :
```
├── genome_diversity_of_maketube_genomes			# source folder
│   ├── REF
│   ├── DELTAS							# nucmer intermediary files .delta, .snps, etc.
│   ├── artificial_strains					# snpmutator and maketube genomes
│   ├── REAL_STRAINS						# natural strains sequences (& IS6110 locations, unused)
│   ├── maketube_run						# maketube files : sequence, evolution partitions, annotations files. 
|   │   ├── AF2122
│   │   ├── 18b
│   │   └── H37Rv
│   ├────── distance_to_reference.tsv			# bases aligned on the reference, reference aligned on the sample genome
│   └────── maketube_with_natural_genomes.sh		# script with all the analysis
├────────── genome_pairwise_comparison.R			# visualisation
├────────── maketube.R
└────── natural_strains_info.txt				#ID, lineage and additional information about natural strains
```


