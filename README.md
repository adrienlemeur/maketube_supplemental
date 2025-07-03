<h1 align="center"> Maketube supplemental </h1>

This is the repository for the maketube article (link).

#### Table of contents

<!--ts-->
-  [Figure 3 : Diversity of artificial and natural genomes](#genome_diversity_tag)
-  [Figure 5 : Comparing three variant caller](#three_vc)
-  [Figure 4 : Nucmer and minimap2](#nucmer_minimap2)
-  [Figure 6 : Impact of structural variants](#structural_variants)

---
## <a name="genome_diversity_tag"></a> Diversity of artificial and natural genomes

### Dependencies :

#### Tools :
	* snpmutator (v1.2.0, installed via mamba)
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
The visualisation is made from `distance_to_reference.tsv` with `genome_pairwise_comparison.R`.
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
   ├── 0_alignement_30X_run.sh				 		# create illumina FASTQ from fasta art
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
   ├── three_variant_caller_results.tsv					# piped output from 1_COMPARING_3_VARIANT_CALLERS.sh
   └── THREE_VC.R							# visualisation
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
	Artificial genomes [used for comparing genotube, TBprofiler and MTBseq](#three_vc) were also used to evaluate the performance of nucmer and minimap2.
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
├── 6_minimap_nucmer_variants.sh			# from fasta sequence to alignment to VCF format
├── 7_COMPARING_3_VARIANT_CALLERS.sh			# compute the performance metrics from the VCF computed in script 6
├── a_database_to_rule_them_all.tsv			# information about the strains (source, artificial or natural, maketube or snpmutator, etc.)
├── nucmer_minimap2_results.tsv				# piped results from script 7
└── NUCMER_PAF.R					# visualisation & statistical analysis
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
	Variant are compared to the test VCF using vcf2metrics.py

</p>
```
local_impact_variant_calling/
├── 0_running_maketube.sh						#
├── 1_fastq_alignment_variant_calling.sh
├── 2_actually_calling_variants.sh
├── 3_filtering_variants.sh
├── BAM
├── FILTERED_VCF
├── maketube.R
├── maketube_run
├── REF
├── results_by_regions_with_without_duplicated.every_at_10k.bed
├── structural_variants_metrics.R
├── structural_variants_metrics_vs_real_regions.R
├── VCF
└── vcf2metrics.py
```
