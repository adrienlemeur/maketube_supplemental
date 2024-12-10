
You will find in this repository most of the files used in Le Meur - 2024

This repository contains :

	- All input data (reference genomes) as well as most of the intermediary files (to the exeption of FASTQ and BAM for file size limitations).
	- All command lines and scripts used for generating those files, with their exact parameters for this study
	- Final table used for building the different figures presented in this study, as well as the R script used

If you have any question / observation, feel free to send me an email at :
> alemeur@biophylo.com

### Complete folder and file description
```
├── maketube800								# First batch of strains with ~ 800 SNP. Fig
│   ├── ARTIFICIALY_EVOLVED_STRAINS					# Maketube & SNPmutator strains 
│   │   ├── MAKETUBE
│   │   │   ├── FASTA							# Mutated sequence with structural variants 
│   │   │   ├── PARTITION						# Partition file for backtracking the position
│   │   │   ├── SV							# Structural variants on the source genome
│   │   │   └── VCF							# True variants
│   │   └── SNPMUTATOR
│   │       ├── FASTA							# Mutated sequence
│   │       ├── info							# Informations on the strain construction
│   │       └── VCF							# True variants
│   ├── DATA_AND_VISUALISATION						# Figure and scripts for creating the figures
│   │   ├── a_database_to_rule_them_all.tsv				# Strain info
│   │   ├── all_strains_info.tsv					# Strain info
│   │   ├── dnadiff_metrics.tsv						# Grepped lines of dnadiff output from the delta files (/maketube800/DELTAS)
│   │   ├── Fig3_dotplot_dist_to_H37Rv.R				# Code for figure 3 & stats
│   │   ├── Fig4_MINIMAP2_PAF_plot.R					# Code for figure 4 & stats
│   │   ├── Fig_5_THREE_VARIANT_CALLER.R				# Code for Fig 5 & stats
│   │   └── maketube_metrics.tsv					# Data (True Positive, False Positive, False Negative) for Genotube, TBprofiler and MTBseq for snpmutator and maketube genomes
│   ├── DELTAS								# DNA diff output for both natural and artificial strains
│   ├── NATURAL_STRAINS							# Genome sequence of natural strains
│   ├── SCRIPT								# Bash script for generating FASTQ, align reads on the ref. genome, variant calling, and generating various metrisc
│   │   ├── 0_FASTA_TO_METRICS.sh					# Generating delta files and metrics from natural and artificial strains
│   │   ├── 1_FASTA_TO_BAM.sh						# FASTQ construction and alignement on the reference genome
│   │   ├── 2_BAM_TO_VCF_FREEBAYES.sh					# Calling variants with freebayes and filtering them using the Genotube filters
│   │   ├── 2_BAM_TO_VCF_GATK.sh					# Calling variants with GATK and filtering them using the MTBseq filters
│   │   ├── 2_BAM_TO_VCF_MPILEUP.sh					# Calling variants with mpileup and filtering them using the TBprofilers filters
│   │   ├── 3_FASTA_TO_VCF.sh						# Calling variants with nucmer and minimap2
│   │   ├── 4_VCF_TO_METRICS.sh						# Computation of Precision and Recall for Maketube and SNPmutator strains
│   │   ├── all2vcf							# py script for converting nucmer files to vcf
│   │   ├── snpmutator_genomes_genesis.sh				# snpmutator parameters
│   │   └── vcf2metrics.py						# maketube script for computing precision and recall of both Maketube and SNPmutator strains
│   └── VARIANT_CALLING_RESULTS						# Result of variant calling. For Maketube, only H1 and H10 variants were considered
│       ├── FREEBAYES_RAW_VCF						# Freebayes raw calls before filtering
│       ├── GATK_VCF							# GATK raw calls and filtered calls
│       ├── GENOTUBE_VCF						# Freebayes filtered calls
│       ├── MINIMAP2_VCF						# Variants called by Minimap2
│       ├── NUCMER_VCF							# Variants called by nucmer (with all2vcf)
│       └── VCF_MPILEUP							# Mpileup raw calls
│           └── F3							# Mpileup filtered calls to keep only variants
│
│
│
├── maketube2500							# Second batch of maketube strains with more variants to study the impact of structural variants
│   ├── ARTIFICIAL_EVOLVED_STRAINS					# Maketube strains
│   │   └── MAKETUBE
│   │       ├── FASTA							# Mutated sequence with structural variants
│   │       ├── PARTITIONS						# Partition file for backtracking the position
│   │       ├── SV							# Structural variants on the source genome
│   │       └── VCF							# True variants
│   ├── DATA_AND_VISUALISATION						# Figure and scripts for creating the figures
│   │   ├── A_metrics_by_regions_by_SV.R				# Creating the graph of precision and recall in the vicinity of structural variants
│   │   ├── B_statistical_exploitation.R				# Creating the graph of precision and recall within regions of interest
│   │   ├── maketube_genomes_genotube_performance_metrics_by_regions.tsv	# Performance metrics for maketube genomes called with freebayes
│   │   ├── maketube_genomes_TBprofiler_performance_metrics_by_regions.tsv	# Performance metrics for maketube genomes called with mpileup
│   │   ├── precision_and_recall_region.svg				# Figure 7 : Graph of precision and recall within regions of interest
│   │   └── precision_and_recall_SV.svg					# Figure 6 : Graph of precision and recall in the vicinity of structural variants
│   ├── SCRIPTS
│   │   ├── 0_running_maketube.sh					# generating maketube strains
│   │   ├── 1_fastq_alignment_variant_calling.sh			# generating fastq and alignment on the reference genome
│   │   ├── 2_actually_calling_variants.sh				# variant calling
│   │   ├── 3_filtering_variants.sh					# filtering variants
│   │   └── 4_SV_region.sh						# Precision and recall within SV and regions of interest
│   └── VARIANT_CALLING_RESULTS						# TBprofiler and Maketube variants calling
│       ├── GENOTUBE_FREEBAYES_FILTERED
│       ├── GENOTUBE_FREEBAYES_RAW_CALLS
│       ├── TBPROFILER_MPILEUP_FILTERED
│       └── TBPROFILER_MPILEUP_RAW_CALL
└── README.txt								# You are here !
```
