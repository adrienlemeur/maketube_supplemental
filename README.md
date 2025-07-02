<h1 align="center"> Maketube supplemental </h1>

This is the repository for the maketube article (link).

#### Table of contents

<!--ts-->
-  [Figure 3 : Diversity of artificial and natural genomes](#genome_diversity_tag)


## <a name="genome_diversity_tag"></a> Diversity of artificial and natural genomes

### Dependencies :

#### Tools :
	* snpmutator
	* maketube (included)

#### R libs :
	* seqinr, jackalope, optparse, Biostrings, dplyr, ape (maketube)
	* ggplot2, scales, ggrepel, ggforce

### Descriptions :

<p align="justify">
Maketube and snpmutator are first launched to generate strains from the 3 reference genomes (`natural_strains_info.txt`).
A subsample of 10 strains for each reference, was manually copied into the `artificial_strains` folder.
These genomes, and natural genomes are aligned onto the three reference genomes using nucmer, and the resulting delta files are filtered. A report is generated with dnadiff.
The number of base of the considered genome aligned on the reference, as well as the number of base of the reference aligned on the considered genome, are grepped and conserved in distance_to_reference.tsv.
The visualisation is made from `distance_to_reference.tsv` with `genome_pairwise_comparison.R`
</p>

### File :
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
## <a name="nucmer_minimap2"></a> Evaluation of TBprofiler, MTBseq & genotube
### Dependencies :

#### Tools : 
	* all2vcf, nucmer, minimap2, bcftools

#### R libs :
	* ggplot2, gridExtra, scales, patchwork, gghalves, ggpp, ggpubr

### Descriptions :

<p align="justify">
Artificial genomes were 
</p>




## <a name="nucmer_minimap2"></a> Evaluating the performance of nucmer & minimap2
### Dependencies :

#### Tools : 
	- all2vcf
	- nucmer
	- minimap2
	- bcftools

#### R libs :
	- ggplot2
	- gridExtra
	- scales
	- patchwork
	- gghalves
	- ggpp
	- ggpubr

#### Descriptions :

<p align="justify">

Artificial genomes were 
</p>

```
nucmer_vs_minimap2/
├── DELTAS								# nucmer deltas & temporary files
├── HEUPINK_STRAINS								# snpmutator genomes
├── maketube_strains								# maketube genomes
├── PAF								# PAF alignments
├── PAF_VCF								# PAFtools.js produced VCF
├── PSEUDOVCF								# nucmer derived VCF
├── 6_minimap_nucmer_variants.sh								# generate VCF from nucmer & minimap2 files
├── 7_COMPARING_3_VARIANT_CALLERS.sh								# compute the performance metrics from the VCF
├── a_database_to_rule_them_all.tsv								# info about the 
├── nucmer_minimap2_results.tsv								# results from script 7
└── NUCMER_PAF.R								# visualisation
```
