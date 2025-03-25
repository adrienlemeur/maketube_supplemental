#!/usr/bin/env python3

from cyvcf2 import VCF
import numpy as np
import argparse
import re, sys, gc, os

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(
			prog="TB-detective",
			description="""TB-detective is a script to identify the lineage, sublineage and antibiotic resistance of a Mycobacterium tuberculosis sample from a VCF annotated with snpEff. It was written in python with cyvcf2.""",
			epilog="Written by Adrien Le Meur, v.1.2")

parser.add_argument('-i', type=str, nargs='+', required=True, help='a single sample VCF aligned on H37Rv genome')
parser.add_argument('-lin', type=str, nargs=1, required=False, help='a tab separated table with 1-based SNP position, the ALT nucleotide and the associated lineage')
parser.add_argument('-ab', type=str, nargs=1, required=False, help='a tab separated table with genes, mutations (snpEff prot. or nuc. mutation annotation), associated antibiotic resistance and confidence threshold')
parser.add_argument('-cf', type=int, default=2, help='threshold level of confidence (only AMR associated variants with the same or higher confidence will be considered)')

args = parser.parse_args()

lineage = args.lin
input_vcf = args.i
antibio = args.ab

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

def flatten(xs):
    result = []
    if isinstance(xs, (list, tuple)):
        for x in xs:
            result.extend(flatten(x))
    else:
        result.append(xs)
    return result

#antibioresistance & lineage information are stored in two dictionaries for fast accession
amr_dict = {}
sample_antibio_dict = {}

lineage_dict = {}
for line in open(lineage[0]):
	line = line.rstrip().split("\t")
	lineage_dict[line[0]] = {}
	lineage_dict[line[0]]["VAR"] = line[1]
	lineage_dict[line[0]]["LIN"] = line[2]

for line in open(antibio[0]):
	line = line.rstrip().split("\t")
	if(int(line[3]) <= args.cf):
		key = "_".join( [ line[1], line[0] ] )
		amr_dict[key] = {}
		amr_dict[key]["gene"] = line[0]
		amr_dict[key]["position"] = line[1]
		amr_dict[key]["resistance"] = line[2].lower()
		amr_dict[key]["confidence"] = int(line[3])
		sample_antibio_dict[line[2].lower()] = "FALSE" #additional dict with antibiotic resistance class as key and TRUE (resistant) or FALSE (sensitive)


#print the header
print("sample", "signal_type", "lineages", "sublineages", "drug_res_type", "resistance_count", "\t".join(sample_antibio_dict.keys()), sep = "\t")

if os.path.exists("all_samples_full_barcode.txt"):
	os.remove("all_samples_full_barcode.txt")

if os.path.exists("all_samples_AB.txt"):
	os.remove("all_samples_AB.txt")

#print(amr_dict.keys())
#sys.exit()
#for every VCF input in -i
for vcf_file in input_vcf:
	vcf = VCF(vcf_file)

	if vcf.samples[0] == "unknown":
		sample = os.path.basename(vcf_file).split(".")[0]
	else:
		sample = vcf.samples[0]

	eprint("Doing "+sample+" ...")

	sample_lineage_dict = {}

	#printing all lineage SNP found
	LBC = open("all_samples_full_barcode.txt", "a")

	antibio_count = 0

	#printing all antibiotic SNP found
	AB = open("all_samples_AB.txt", "a")

	for i in sample_antibio_dict:
		sample_antibio_dict[i] = {}
		sample_antibio_dict[i]["CLASS"] = 'FALSE'
		sample_antibio_dict[i]["confidence"] = -1

	for v in vcf:

		if(v.FILTER is None):
			FILTER = 'PASS'
		else:
			FILTER = v.FILTERS[0]

		if v.INFO.get("ANN"):
			ANN = v.INFO["ANN"].split("|")

			if FILTER == 'PASS':
				protein_key = ["_".join( [ ANN[i], ANN[i-7] ]) for i, v in enumerate(ANN) if re.match("^p\.", v)]
				nucleotide_key = ["_".join( [ ANN[i], ANN[i-6] ]) for i, v in enumerate(ANN) if re.match("^n\.", v)]
				upstream_key = ["_".join( [ ANN[i], ANN[i-6] ]) for i, v in enumerate(ANN) if re.match("^c\.", v)]

				for one_key in protein_key + nucleotide_key + upstream_key:
					if(one_key in amr_dict):
						AB.write(sample+"\t"+amr_dict[one_key]["gene"]+"\t"+amr_dict[one_key]["position"]+"\t"+FILTER+"\t"+amr_dict[one_key]["resistance"]+"\t"+str(amr_dict[one_key]["confidence"])+'\n')
						if(sample_antibio_dict[amr_dict[one_key]["resistance"]]["CLASS"] == "FALSE"):
							antibio_count+=1
							sample_antibio_dict[amr_dict[one_key]["resistance"]]["CLASS"] = "TRUE"
						if(amr_dict[one_key]["confidence"] < sample_antibio_dict[amr_dict[one_key]["resistance"]]["confidence"] or sample_antibio_dict[amr_dict[one_key]["resistance"]]["confidence"] == -1):
							sample_antibio_dict[amr_dict[one_key]["resistance"]]["confidence"] = amr_dict[one_key]["confidence"]

		if v.is_snp & (FILTER == 'PASS'):
			variant_pos = str(v.POS)

			if(variant_pos in lineage_dict):
				if lineage_dict[variant_pos]["LIN"] not in sample_lineage_dict:
					sample_lineage_dict[lineage_dict[variant_pos]["LIN"]] = 1
				else:
					sample_lineage_dict[lineage_dict[variant_pos]["LIN"]] = sample_lineage_dict[lineage_dict[variant_pos]["LIN"]] + 1
				LBC.write(sample+"\t"+variant_pos+"\t"+v.REF+"\t"+lineage_dict[variant_pos]["VAR"]+"\t"+lineage_dict[variant_pos]["LIN"]+"\t"+FILTER+"\n")
	vcf.close()
	AB.close()
	LBC.close()

	if len(sample_lineage_dict) == 0:
		signal = "NO_SNP"
		main_sublineages = ["NONE"]
		main_lineages = ["NONE"]
	else:
		sample_lineage_list = list(sample_lineage_dict)
		sample_lineage_list.sort(reverse=True, key=lambda x: (len(x), x))
		main_sublineages = [ sample_lineage_list.pop(0) ]

		while len(sample_lineage_list) > 0:
			new = True
			first_key = sample_lineage_list[0]
			for i in main_sublineages:
				if re.search(first_key, i) is not None:
					new = False
			if new == True:
				main_sublineages.append(first_key)
			sample_lineage_list.pop(0)

			if len(main_sublineages) == 1:
				signal = "CLEAR"
			elif len(main_sublineages) == 2:
				signal = "COINFECTION"
			else:
				signal = "CONTAMINATION"

		main_lineages = set([ i.split('.')[0] for i in main_sublineages ])

	DR_type = "SENSITIVE"
	antibiolist = [x for x in sample_antibio_dict]

	if("isoniazid" in antibiolist and "rifampicin" in antibiolist):
		if(sample_antibio_dict["isoniazid"]["CLASS"] == "TRUE" and sample_antibio_dict["rifampicin"]["CLASS"] == "TRUE"):
			DR_type = "MDR"
			if("fluoroquinolones" in antibiolist):
				if(sample_antibio_dict["fluoroquinolones"]["CLASS"] == "TRUE"):
					if("capreomycin" in antibiolist):
						if(sample_antibio_dict["capreomycin"]["CLASS"] == "TRUE"):
							DR_type = "XDR"
					if("kanamycin" in antibiolist):
						if(sample_antibio_dict["kanamycin"]["CLASS"] == "TRUE"):
							DR_type = "XDR"
					if("amikacin" in antibiolist):
						if(sample_antibio_dict["amikacin"]["CLASS"] == "TRUE"):
							DR_type = "XDR"
		elif "rifampicin" in antibiolist:
			if sample_antibio_dict["rifampicin"]["CLASS"] == "TRUE":
				DR_type = "RR"

	output = "\t".join([sample, signal, ";".join(list(main_lineages)), ";".join(list(main_sublineages)), DR_type, str(antibio_count), '\t'.join(str(sample_antibio_dict[x]["CLASS"]+" ("+str(sample_antibio_dict[x]["confidence"])+") ") for x in sample_antibio_dict)])
	print(re.sub("FALSE \(-1\)", "FALSE", output))

