<h1 align="center"> Maketube supplemental </h1>

This is the repository for the maketube article (link).

#### Table of contents

<!--ts-->
-  [Figure 3 : Diversity of artificial and natural genomes](#genome_diversity_tag)
-  [Figure 5 : Comparing three variant caller](#three_vc)
-  [Figure 4 : Nucmer and minimap2](#nucmer_minimap2)

---
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
The number of base of the considered genome aligned on the reference, as well as the number of base of the reference aligned on the considered genome, are grepped and conserved in `distance_to_reference.tsv`.
The visualisation is made from `distance_to_reference.tsv` with `genome_pairwise_comparison.R`
</p>

### File :
```
├── genome_diversity_of_maketube_genomes		# source folder
│   ├── REF
│   ├── DELTAS						# nucmer intermediary files .delta, .snps, etc.
│   ├── artificial_strains				# snpmutator and maketube genomes
│   ├── REAL_STRAINS					# natural strains sequences (& IS6110 locations, unused)
│   ├── maketube_run					# maketube files : sequence, evolution partitions, annotations files. 
|   │   ├── AF2122
│   │   ├── 18b
│   │   └── H37Rv
│   ├────── distance_to_reference.tsv			# bases aligned on the reference, reference aligned on the sample genome
│   └────── maketube_with_natural_genomes.sh		# script with all the analysis
├────────── genome_pairwise_comparison.R		# visualisation
├────────── maketube.R
└────── natural_strains_info.txt			#ID, lineage and additional information about natural strains
```
---
## <a name="three_vc"></a> Evaluation of TBprofiler, MTBseq & genotube
### Dependencies :

#### Tools : 
	* all2vcf, nucmer, minimap2, bcftools

#### R libs :
	* ggplot2, gridExtra, scales, patchwork, gghalves, ggpp, ggpubr

### Descriptions :
Variants were called 


```
└─ genotube_tbprofiler_mtbseq
   ├── 0_alignement_30X_run.sh				 		# generate the FASTQ from the fasta using art
   ├── 1_COMPARING_3_VARIANT_CALLERS.sh					# compare the variants called by the pipelines to the reference VCF
   ├── a_database_to_rule_them_all.tsv					# information about strains
   ├── all_strains_info.tsv						# reduced strain information
   ├── fig5_precision_recall_3VC_only_violin.svg			# pipeline performances with duplication region
   ├── fig5_precision_recall_3VC_only_violin_wo_dupli.svg		# pipeline performances without duplication region
   ├── FREEBAYES_RAW_CALLS						# freebayes raw call without filtering
   ├── GATK_MTBseq							# MTBseq VCF
   ├── GATK_RAW_CALLS							# GATK raw call without additional filtering
   ├── genotube								# genotube VCF
   ├── TBprofiler							# TBprofiler VCF
   ├── HEUPINK_STRAINS							# snpmutator fasta
   ├── maketube_strains							# maketube fasta
   ├── three_variant_caller_results.tsv					# output from 1_COMPARING_3_VARIANT_CALLERS.sh
   └── THREE_VC.R							# visualisation
```
<p align="justify">
Artificial genomes were 
</p>


---
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
├── DELTAS						# nucmer deltas & temporary files
├── HEUPINK_STRAINS					# snpmutator genomes
├── maketube_strains					# maketube genomes
├── PAF							# PAF alignments
├── PAF_VCF						# PAFtools.js produced VCF
├── PSEUDOVCF						# nucmer derived VCF
├── 6_minimap_nucmer_variants.sh			# generate VCF from nucmer & minimap2 files
├── 7_COMPARING_3_VARIANT_CALLERS.sh			# compute the performance metrics from the VCF
├── a_database_to_rule_them_all.tsv			# info about the 
├── nucmer_minimap2_results.tsv				# results from script 7
└── NUCMER_PAF.R					# visualisation
```
