<h1 align="center"> Maketube supplemental </h1>

This is the repository for the maketube article (link).

#### Table of contents

<!--ts-->
-  [Figure 3 : Diversity of artificial and natural genomes](#genome_diversity_tag)
-  [Figure 4 : Nucmer and minimap2](#nucmer_minimap2)
-  [Figure 5 : Comparing three variant caller](#three_vc)
-  [Figure 6 : Impact of structural variants](#structural_variants)
---
## <a name="forewords"></a> Forewords

This repository contains all the files and script used in the maketube genome article.
Heavy files such as 30X FASTQ and bam from `fig6_local_impact_variant_calling` were omitted.
They can be generated straight away using the included script.

Due to the lack of seed system for maketube, re-doing a run will not produce the same exact genomes.
The strains generated for the article are generated in [Figure 3](#genome_diversity_tag).
The strains generated from H37Rv are then reused for [Figure 4](#nucmer_minimap2) and [Figure 5](#three_vc).
A new set of strains, with different maketube parameters are generated in [Figure 6](#structural_variants).

If you encounter an issue when trying to reproduce the analysis, or if you have a question, please do send me an email at alemeur at biophylo.com so I can look into it. 

---
## <a name="genome_diversity_tag"></a> Diversity of artificial and natural genomes

### Dependencies :

#### Tools :
	* snpmutator (v1.2.0, installed via mamba)
	* maketube (included)

#### R libs :
	* seqinr, jackalope, optparse, Biostrings, dplyr, ape (maketube)
	* ggplot2, scales, ggrepel, ggforce (INSTALL WITH DEVTOOLS v0.5)

### Descriptions :

<p align="justify">

Artificial strains from maketube and snpmutator are first generated.
Maketube and snpmutator are first launched to generate strains from the 3 reference genomes (`natural_strains_info.txt`).
A subsample of 10 strains for each reference, was manually copied into the `artificial_strains` folder.
These genomes, and natural genomes are aligned onto the three reference genomes using nucmer, and the resulting delta files are filtered. A report is generated with dnadiff.
The number of base of the considered genome aligned on the reference, as well as the number of bases of the reference aligned on the considered genome, are grepped and conserved in `distance_to_reference.tsv`.
The visualisation is made from `distance_to_reference.tsv` with `genome_pairwise_comparison.R`.
</p>

### File :
```
└── genome_diversity_of_maketube_genomes		# source folder
    ├── REF
    ├── DELTAS						# nucmer intermediary files .delta, .snps, etc.
    ├── artificial_strains				# snpmutator and maketube genomes
    ├── REAL_STRAINS					# natural strains sequences (& IS6110 locations, unused)
    ├── maketube_run					# maketube files : sequence, evolution partitions, annotations files. 
    │   ├── AF2122
    │   ├── 18b
    │   └── H37Rv
    ├────── distance_to_reference.tsv			# bases aligned on the reference, reference aligned on the sample genome
    ├────── 1_creating_artificial_genomes.sh		# creating the artificial genomes
    ├────── 2_pairwise_genetic_distance_samples_vs_reference.sh		# compare the distance of artificial and natural strains to the references
    ├────── genome_pairwise_comparison.R		# visualisation
    ├────── maketube.R
    └── natural_strains_info.txt			#ID, lineage and additional information about natural strains
```

---
## <a name="nucmer_minimap2"></a> Evaluating the performance of nucmer & minimap2
### Dependencies :

#### Tools : 
	- all2vcf (v0.7.8, installed manually), nucmer (v3.1, installed manually), minimap2 (v2.28-r1209, compiled manually), bcftools (v1.19, htslib : v1.19, installed manually)

#### R libs :
	- ggplot2, gridExtra, scales, patchwork, gghalves, ggpp, ggpubr, flipr
#### Python libs :
	* cyvcf2 (vcf2metrics.py)

#### Descriptions :

<p align="justify">
	Artificial genomes [previously generated](#genome_diversity_tag) were also used to evaluate the performance of nucmer and minimap2.
	First the genomes are aligned unto H37Rv using minimap2 and nucmer. Variants were called from minimap2 paf with paftools into `PAF_VCF`. 
	nucmer variants were called using show-snps and converted into VCF in `PSEUDOVCF` with all2vcf (`6_minimap_nucmer_variants.sh`).
	Variants called were compared to reference VCF in `maketube_strains` and `HEUPINK_STRAINS` using vcf2metrics.py. 
</p>

```
nucmer_vs_minimap2/
├── DELTAS						# nucmer deltas alignment files & other nucmer files
├── HEUPINK_STRAINS					# snpmutator genomes (fasta, reference VCF)
├── maketube_strains					# maketube genomes (fasta, backtrack, annotations and reference VCF)
├── PAF							# PAF alignments
├── PAF_VCF						# PAFtools.js produced VCF from script 6
├── PSEUDOVCF						# nucmer derived VCF using all2vcf
├── 1_call_variants_nucmer_minimap2.sh			# call variants from the artificial sequences FASTA
├── 2_compute_nucmer_minimap2_performance_metrics.sh			# compare the VCF created by script n°1 to the reference VCF
├── a_database_to_rule_them_all.tsv			# information about the strains (source, artificial or natural, maketube or snpmutator, etc.)
├── nucmer_minimap2_results.tsv				# piped results from script 7
└── NUCMER_PAF.R					# visualisation & statistical analysis
```

---
## <a name="three_vc"></a> Evaluation of TBprofiler, MTBseq & genotube
### Dependencies :

#### Tools : 
	* art_illumina (v2.5.8, installed via apt) , bcftools (v1.19, htslib : v1.19, installed manually)

#### R libs :
	* ggplot2, gridExtra, tidyr, flipr, patchwork
	* seqinr, jackalope, optparse, Biostrings, dplyr, ape (maketube)

#### Python libs :
	* cyvcf2 (vcf2metrics.py)

### Descriptions :

<p align="justify">
	Artificial genomes were created using maketube and snpmutator and respectively put in `maketube_strains` and `HEUPINK_STRAINS`.
	Their variants were called using TBprofiler (v6.6.3, installed via mamba), MTBseq (v1.1.0, installed via mamba) and genotube (v???, installed manually).
	VCF from TBprofiler were put in `TBprofiler`, VCF from MTBseq were put in `GATK_MTBseq` and VCF from genotube in `genotube`. 
	Raw calls, ie. unfiltered call from GATK and freebayes are respectively in `GATK_RAW_CALLS` and `FREEBAYES_RAW_CALLS`.
	Variants are then compared to the reference variants using `vcf2metrics.py` in script `1_COMPARING_3_VARIANTS_CALLERS.sh`.

</p>


```
└─ genotube_tbprofiler_mtbseq
   ├── 1_three_pipeline_variant_metrics.sh					# compare the variants called by the pipelines to the reference VCF
   ├── a_database_to_rule_them_all.tsv					# strain information and file location
   ├── all_strains_info.tsv						# reduced strain information
   ├── fig5_precision_recall_3VC_only_violin.svg			# pipeline performances with duplication region
   ├── fig5_precision_recall_3VC_only_violin_wo_dupli.svg		# pipeline performances without duplication region
   ├── FREEBAYES_RAW_CALLS						# freebayes raw call
   ├── GATK_MTBseq							# MTBseq VCF
   ├── GATK_RAW_CALLS							# GATK raw call
   ├── genotube								# genotube VCF
   ├── TBprofiler							# TBprofiler VCF
   ├── HEUPINK_STRAINS							# snpmutator genomes
   ├── maketube_strains							# maketube genomes
   ├── three_variant_caller_results.tsv					# piped output of script n°1. TSV with performance metrics & additional informations
   └── THREE_VC.R							# script for three_variant_caller_results.tsv visualisation
```


---
## <a name="structural_variants"></a> Impact of structural variants
### Dependencies :

#### Tools : 
	* bcftools (v1.19, htslib : v1.19, installed manually), gnu-parallel (v20250422, installed manually)

#### R libs :
	* seqinr, jackalope, optparse, Biostrings, dplyr, ape (maketube)
	* ggplot2, gridExtra, scales, grDevices, flipr, patchwork, gghalves, ggpp

#### Python libs :
	* cyvcf2 (vcf2metrics.py)

<p align="justify">
	Maketube genomes are generated by script 0.
	FASTQ are generated from the fasta and the reads are aligned on the reference genome by script 1.
	Variant are called using freebayes by script 2.
	Variant are compared to the test VCF using `vcf2metrics.py`. 10kbp regions composed of 1kbp windows are taken across the different regions.
	vcf2metrics.py computes the precision and recall specificaly in these 10kbp regions windows.
</p>

```
local_impact_variant_calling/
├── 0_running_maketube.sh					# creation of maketube genome
├── 1_fastq_alignment_variant_calling.sh			# fastq creation and alignment on H37Rv
├── 2_actually_calling_variants.sh				# variant calling
├── 3_filtering_variants.sh					# comparison of the called variant to the true variants
├── BAM								# bwa-mem2 alignment files on H37Rv from script n°1
├── FILTERED_VCF							# filtered vcf from script n°3
├── maketube.R
├── maketube_run							# maketube genomes (fasta, backtrack, annotations and reference VCF)
├── REF								# reference sequence and source annotation
├── results_by_regions_with_without_duplicated.every_at_10k.bed	# 
├── structural_variants_metrics.R				# visualisation using artificial precision and recall
├── structural_variants_metrics_vs_real_regions.R	# visualisation using the 10k region across the different regions
├── VCF								# unfiltered VCF from script n°2
└── vcf2metrics.py
```
