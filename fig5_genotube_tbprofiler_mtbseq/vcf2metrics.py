#!/usr/bin/env python3

from cyvcf2 import VCF
import numpy as np
import argparse
import re, sys, gc, os
import math

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(
			prog="vcf2metrics",
			description=""" vcf2metrics compares a sample vcf to the reference VCF built by maketube""",
			epilog="Written by Adrien Le Meur, v??")

parser.add_argument('-i', type=str, nargs='+', required=True, help='folder containing VCF to test against the reference')
parser.add_argument('--reference', type=str, nargs='+', required=True, help='reference vcf')
parser.add_argument('--bed', type=str, nargs='+', required=False, default=None, help='bed file with the region annotations')
parser.add_argument('--backtrack', type=str, nargs='+', required=False, help='bed backtrack file produced by maketube')
parser.add_argument('--subtract', type=str, nargs='+', required=False, help='bed file to remove regions')

parser.add_argument('--trace', action='store_true')

parser.add_argument('--add_col', nargs='+', required=False, help='Space separated list of columns specifying different parameters of the experiment if needed')


#parser.add_argument('--pipeline', type=str, nargs=1, required=False, default="NA", help='single string stating the pipeline name')
#parser.add_argument('--depth', type=str, nargs=1, required=False, default="NA", help='single string stating the coverage depth')
#parser.add_argument('--SV', type=str, nargs=1, required=False, default="NA", help='single string stating the structural variant')
parser.add_argument('--sample', type=str, nargs=1, required=False, default="NA", help='single string stating the sample name')
parser.add_argument('-o', type=str, nargs='+', required=False, default = "./", help='output directory')

def sort_dict(my_dict):
	return(dict(sorted(my_dict.items(), key=len, reverse=True)))

def bed2dict(files):
	bed_dict = {}

	for file in files:
		with open(file) as bed_lines:
			for line in bed_lines:
				line = line.rstrip('\n').split("\t")

				if line[3] not in bed_dict:
					bed_dict[line[3]] = []
				bed_dict[line[3]].append([ min(line[1],line[2]), max(line[1],line[2]) ])
	return(bed_dict)

def vcf2dict(vcf, correction:int):
	vcf_dict = {}

	#already bed dict key, so de facto unique
	for v in vcf:
		#should output the variant anyway
		if v.FILTERS == [] or v.FILTERS[0] == 'PASS':
			variant_start = int(v.POS) - correction
			ALT = v.ALT[0]

			variant_key = "_".join( [ str(variant_start), v.REF, ALT])


			vcf_dict[variant_key] = {}
			vcf_dict[variant_key]['START'] = variant_start
			vcf_dict[variant_key]['CLASS'] = v.var_type
			vcf_dict[variant_key]['END'] = variant_start+len(ALT)
			vcf_dict[variant_key]['REF'] = v.REF
			vcf_dict[variant_key]['ALT'] = ALT
	return(vcf_dict)


def vcfbed2dict4table(vcf_dict_list:list, ref_dict_list:list, backtrack_list:list, variant_class:list, additional_columns):
	my_merge_dictionnary = {}

	for vcf_dict in ref_dict_list:
		for variant in vcf_dict:
			if vcf_dict[variant]['CLASS'] in variant_class:
				variant_new_key = variant+"_ref"
				my_merge_dictionnary[variant_new_key] = {}
				my_merge_dictionnary[variant_new_key]['START'] = vcf_dict[variant]['START']
				my_merge_dictionnary[variant_new_key]['END'] = vcf_dict[variant]['END']
				my_merge_dictionnary[variant_new_key]['PREVIOUS_START'] = "REF"
				my_merge_dictionnary[variant_new_key]['LENGTH'] = vcf_dict[variant]['END'] - vcf_dict[variant]['START']
				my_merge_dictionnary[variant_new_key]['CLASS'] = vcf_dict[variant]['CLASS']
				my_merge_dictionnary[variant_new_key]['REF'] = vcf_dict[variant]['REF']
				my_merge_dictionnary[variant_new_key]['ALT'] = vcf_dict[variant]['ALT']
				my_merge_dictionnary[variant_new_key]['STATUS'] = "NA"


	for vcf_dict in vcf_dict_list:
		for variant in vcf_dict:
			if vcf_dict[variant]['CLASS'] in variant_class:
				my_merge_dictionnary[variant] = {}
				my_merge_dictionnary[variant]['START'] = vcf_dict[variant]['START']
				my_merge_dictionnary[variant]['END'] = vcf_dict[variant]['END']
				if 'PREVIOUS_START' in vcf_dict[variant]:
					my_merge_dictionnary[variant]['PREVIOUS_START'] = vcf_dict[variant]['PREVIOUS_START']
				else:
					my_merge_dictionnary[variant]['PREVIOUS_START'] = vcf_dict[variant]['PREVIOUS_START']
				my_merge_dictionnary[variant]['LENGTH'] = vcf_dict[variant]['END'] - vcf_dict[variant]['START']
				my_merge_dictionnary[variant]['CLASS'] = vcf_dict[variant]['CLASS']
				my_merge_dictionnary[variant]['REF'] = vcf_dict[variant]['REF']
				my_merge_dictionnary[variant]['ALT'] = vcf_dict[variant]['ALT']
				my_merge_dictionnary[variant]['STATUS'] = "NA"

	for backtrack_dict in backtrack_list:
		for region in backtrack_dict:
			if backtrack_dict[region]['class'] == "insertion":
				my_merge_dictionnary[region] = {}
				my_merge_dictionnary[region]['PREVIOUS_START'] = "NA"
				my_merge_dictionnary[region]['REF'] = "NA"
				my_merge_dictionnary[region]['ALT'] = "NA"

				my_merge_dictionnary[region]['START'] = backtrack_dict[region]['ins_start']
				my_merge_dictionnary[region]['END'] = backtrack_dict[region]['ins_stop']
				my_merge_dictionnary[region]['LENGTH'] = my_merge_dictionnary[region]['END'] - my_merge_dictionnary[region]['START']
				my_merge_dictionnary[region]['CLASS'] = "insertion"
				my_merge_dictionnary[region]['STATUS'] = "STRUCT_VAR"

			if backtrack_dict[region]['class'] == "deletion":
				my_merge_dictionnary[region] = {}
				my_merge_dictionnary[region]['PREVIOUS_START'] = "NA"
				my_merge_dictionnary[region]['REF'] = "NA"
				my_merge_dictionnary[region]['ALT'] = "NA"

				my_merge_dictionnary[region]['START'] = backtrack_dict[region]['start']
				my_merge_dictionnary[region]['END'] = backtrack_dict[region]['stop']
				my_merge_dictionnary[region]['LENGTH'] = my_merge_dictionnary[region]['START'] - my_merge_dictionnary[region]['END']
				my_merge_dictionnary[region]['CLASS'] = "deletion"
				my_merge_dictionnary[region]['STATUS'] = "STRUCT_VAR"

			if backtrack_dict[region]['class'] == "transposon":

				initial_pos_key = region+"_initial_position"
				my_merge_dictionnary[initial_pos_key] = {}
				my_merge_dictionnary[initial_pos_key]['PREVIOUS_START'] = "NA"
				my_merge_dictionnary[initial_pos_key]['REF'] = "NA"
				my_merge_dictionnary[initial_pos_key]['ALT'] = "NA"

				my_merge_dictionnary[initial_pos_key]['START'] = backtrack_dict[region]['start']
				my_merge_dictionnary[initial_pos_key]['END'] = backtrack_dict[region]['stop']
				my_merge_dictionnary[initial_pos_key]['LENGTH'] = my_merge_dictionnary[initial_pos_key]['START'] - my_merge_dictionnary[initial_pos_key]['END']
				my_merge_dictionnary[initial_pos_key]['CLASS'] = "transposon_initial_position"
				my_merge_dictionnary[initial_pos_key]['STATUS'] = "STRUCT_VAR"

				final_pos_key = region+"_post_jump_position"
				my_merge_dictionnary[final_pos_key] = {}
				my_merge_dictionnary[final_pos_key]['PREVIOUS_START'] = "NA"
				my_merge_dictionnary[final_pos_key]['REF'] = "NA"
				my_merge_dictionnary[final_pos_key]['ALT'] = "NA"
				
				my_merge_dictionnary[final_pos_key]['START'] = backtrack_dict[region]['ins_start']
				my_merge_dictionnary[final_pos_key]['END'] = backtrack_dict[region]['ins_stop']
				my_merge_dictionnary[final_pos_key]['LENGTH'] = my_merge_dictionnary[initial_pos_key]['END'] - my_merge_dictionnary[initial_pos_key]['START']
				my_merge_dictionnary[final_pos_key]['CLASS'] = "transposon_post_jump"
				my_merge_dictionnary[final_pos_key]['STATUS'] = "STRUCT_VAR"


	for variant in my_merge_dictionnary:
		variant_key = "_".join( [ str(my_merge_dictionnary[variant]['START']), my_merge_dictionnary[variant]['REF'], my_merge_dictionnary[variant]['ALT']])

		if my_merge_dictionnary[variant]['PREVIOUS_START'] == "REF":
			if variant_key in my_merge_dictionnary:
				my_merge_dictionnary[variant]["STATUS"] = "TP"
			else:
				my_merge_dictionnary[variant]["STATUS"] = "FN"
		elif my_merge_dictionnary[variant]['STATUS'] != "STRUCT_VAR":
			if variant_key+"_ref" in my_merge_dictionnary:
				my_merge_dictionnary[variant]["STATUS"] = "TP"
			else:
				my_merge_dictionnary[variant]["STATUS"] = "FP"

	sort_my_merge_dictionnary = dict(sorted(my_merge_dictionnary.items(),key=lambda x: int(x[1]['START'])))

	for key in sort_my_merge_dictionnary:
		if additional_columns is not None :
			print("\t".join(additional_columns), sort_my_merge_dictionnary[key]['START'], sort_my_merge_dictionnary[key]['END'], sort_my_merge_dictionnary[key]['PREVIOUS_START'], sort_my_merge_dictionnary[key]['CLASS'], sort_my_merge_dictionnary[key]['LENGTH'], sort_my_merge_dictionnary[key]['REF'], sort_my_merge_dictionnary[key]['ALT'], sort_my_merge_dictionnary[key]['STATUS'], sep = "\t")
		else:
			print(sort_my_merge_dictionnary[key]['START'], sort_my_merge_dictionnary[key]['END'], sort_my_merge_dictionnary[key]['PREVIOUS_START'], sort_my_merge_dictionnary[key]['CLASS'], sort_my_merge_dictionnary[key]['LENGTH'], sort_my_merge_dictionnary[key]['REF'], sort_my_merge_dictionnary[key]['ALT'], sort_my_merge_dictionnary[key]['STATUS'], sep = "\t")
def compareVCF(reference_vcf_dict, test_vcf_dict, bed):

	comparison_dict = {}

	if bed is None:
		all_bed_region_and_total = [ "TOTAL" ]
	else:
		all_bed_region_and_total = list(bed.keys())
		all_bed_region_and_total.append("TOTAL")

	for region in all_bed_region_and_total:
		comparison_dict[region] = {}
		comparison_dict[region]["snp"] = {}
		comparison_dict[region]["mnp"] = {}
		comparison_dict[region]["indel"] = {}
		comparison_dict[region]["other"] = {}
		comparison_dict[region]["unknown"] = {}

		for variant_type in comparison_dict[region]:
			comparison_dict[region][variant_type] = {}
			comparison_dict[region][variant_type]["TP"] = 0
			comparison_dict[region][variant_type]["FN"] = 0
			comparison_dict[region][variant_type]["FP"] = 0

	for variant in test_vcf_dict:
		if variant in reference_vcf_dict:
			start = int(test_vcf_dict[variant]['START'])
			comparison_dict["TOTAL"][test_vcf_dict[variant]['CLASS']]["TP"] += 1
			if bed:
				for regions in bed:
					for region in bed[regions]:
						if( (start > int(region[0]) and start < int(region[1])) or (start < int(region[0]) and (start > int(region[1]))) ):
							comparison_dict[regions][test_vcf_dict[variant]['CLASS']]["TP"] += 1
							break

		else:
			start = int(test_vcf_dict[variant]['START'])
			comparison_dict["TOTAL"][test_vcf_dict[variant]['CLASS']]["FP"] += 1

			#print(test_vcf_dict[variant])
			if bed:
				for regions in bed:
					for region in bed[regions]:
						if( (start > int(region[0]) and start < int(region[1])) or (start < int(region[0]) and (start > int(region[1]))) ):
							comparison_dict[regions][test_vcf_dict[variant]['CLASS']]["FP"] += 1
							break

	for variant in reference_vcf_dict:
		if variant not in test_vcf_dict:
			start = int(reference_vcf_dict[variant]['START'])
			comparison_dict["TOTAL"][reference_vcf_dict[variant]['CLASS']]["FN"] += 1

			if bed:
				for regions in bed:
					for region in bed[regions]:
						if( (start > int(region[0]) and start < int(region[1])) or (start < int(region[0]) and start > int(region[1])) ):
							comparison_dict[regions][reference_vcf_dict[variant]['CLASS']]["FN"] += 1
							break
	return(comparison_dict)




def read_backtrack(file):
	bed_dict = {}

	with open(file) as bed_lines:
		for line in bed_lines:
			if line.startswith('#'):
				continue
			line = line.rstrip().split()
			key = "_".join([ line[3], line[1], line[2] ])

			if(line[3] == "transposon"):
				bed_dict[key] = {"class":line[3], "start":int(line[1]), "stop":int(line[2]), "ins_start":int(line[4]), "ins_stop":int(line[5])}
			elif(line[3] == "insertion"):
				bed_dict[key] = {"class":line[3], "start":int(line[1]), "stop":int(line[2]), "ins_start":int(line[4]), "ins_stop":int(line[5])}
			elif(line[3] == "deletion"):
				bed_dict[key] = {"class":line[3], "start":int(line[1]), "stop":int(line[2])}
	return(bed_dict)

def subtract_bed_from_vcf(vcf_dict, subtrack):
	new_dict = {}
	for variant_key in list(vcf_dict.keys()):
		vcf_dict[variant_key]["REMOVE"] = 0

		for line in subtrack:
			for interval in subtrack[line]:
				if( (vcf_dict[variant_key]["START"] <= max(int(interval[0]), int(interval[1])) ) and (vcf_dict[variant_key]["START"] >= min(int(interval[0]), int(interval[1])) ) ):
					vcf_dict[variant_key]["REMOVE"] = 1

	for variant_key in vcf_dict:
		if(vcf_dict[variant_key]["REMOVE"] == 0):
			new_dict[variant_key] = vcf_dict[variant_key]
	return(new_dict)

def backtrack(vcf_dict, backtrack_dict):
	new_dict = {}

	for variant_key in vcf_dict:
		ovPOS = vcf_dict[variant_key]["START"]

		for interval in backtrack_dict:
			start = backtrack_dict[interval]["start"]
			stop = backtrack_dict[interval]["stop"]

			if( (vcf_dict[variant_key]["START"] >= stop) and (backtrack_dict[interval]["class"] == "insertion") ):
				#print(ovPOS, sep = "\t")
				ovPOS = ovPOS - backtrack_dict[interval]["ins_start"] + backtrack_dict[interval]["ins_stop"] - 1
				#print(backtrack_dict[interval], vcf_dict[variant_key], ovPOS)

			if( (vcf_dict[variant_key]["START"] >= stop) and (backtrack_dict[interval]["class"] == "deletion") ):
				ovPOS = ovPOS + start - backtrack_dict[interval]["stop"] - 1
			
			if(backtrack_dict[interval]["class"] == "transposon"):
				start = backtrack_dict[interval]["start"]
				stop = backtrack_dict[interval]["stop"]
				insertion_start = backtrack_dict[interval]["ins_start"]
				insertion_stop = backtrack_dict[interval]["ins_stop"]

				if( vcf_dict[variant_key]["START"] <= insertion_stop & vcf_dict[variant_key]["START"] >= insertion_start ):
					#print("a transposon was muted after its jump"
					continue
				#I
				elif((vcf_dict[variant_key]["START"] > start) and (vcf_dict[variant_key]["START"] < insertion_start) and (insertion_start > start) ):
					ovPOS = ovPOS - (stop - start) - 1
				#II
				elif((vcf_dict[variant_key]["START"] < start) and (vcf_dict[variant_key]["START"] > insertion_start) and (insertion_start < start)):
					ovPOS = ovPOS + (insertion_stop - insertion_start)
		new_key = "_".join([ str(ovPOS), vcf_dict[variant_key]["REF"], vcf_dict[variant_key]["ALT"] ])
		new_dict[new_key] = {}
		new_dict[new_key]["START"] = ovPOS
		new_dict[new_key]["END"] = ovPOS + len(vcf_dict[variant_key]["ALT"])
		new_dict[new_key]["CLASS"] = vcf_dict[variant_key]["CLASS"]
		new_dict[new_key]["REF"] = vcf_dict[variant_key]["REF"]
		new_dict[new_key]["ALT"] = vcf_dict[variant_key]["ALT"]
		new_dict[new_key]["PREVIOUS_START"] = vcf_dict[variant_key]["START"]

	return(new_dict)


#THE MAIN STARTS 

args = parser.parse_args()

output_directory = str(args.o[0])


if(args.o):
	output_directory = str(args.o[0])
	if not os.path.exists(output_directory):
		os.mkdir(output_directory)


if args.subtract:
	regions_to_remove = bed2dict(args.subtract)

if args.bed is not None:
	regions_dictionnary = bed2dict(args.bed)
else:
	regions_dictionnary = None

if(args.backtrack):
	backtrack_dict = read_backtrack(args.backtrack[0])

if args.subtract is not None:
	reference_dictionnary = subtract_bed_from_vcf(sort_dict(vcf2dict(VCF(args.reference[0]), 0)), regions_to_remove)
else:
	reference_dictionnary = sort_dict(vcf2dict(VCF(args.reference[0]), 0))

comparison_dict = {}

for sample_vcf in args.i:

	sample_dictionnary = sort_dict(vcf2dict(VCF(sample_vcf), 0))

	if args.subtract is not None and args.backtrack is not None:
		backtracked_sample_vcf_dictionnary = subtract_bed_from_vcf(sort_dict(backtrack(sample_dictionnary, backtrack_dict)), regions_to_remove)
		comparison_dict = compareVCF(reference_dictionnary, backtracked_sample_vcf_dictionnary, regions_dictionnary)
	elif args.backtrack is not None :
		backtracked_sample_vcf_dictionnary = sort_dict(backtrack(sample_dictionnary, backtrack_dict))
		comparison_dict = compareVCF(reference_dictionnary, backtracked_sample_vcf_dictionnary, regions_dictionnary)
	else:
		comparison_dict = compareVCF(reference_dictionnary, sample_dictionnary, regions_dictionnary)

	if(args.trace):
		if args.add_col:
			vcfbed2dict4table(vcf_dict_list = [ backtracked_sample_vcf_dictionnary ], ref_dict_list = [ reference_dictionnary ], backtrack_list = [ backtrack_dict ], variant_class = 'snp', additional_columns = args.add_col )
		else:
			vcfbed2dict4table(vcf_dict_list = [ backtracked_sample_vcf_dictionnary ], ref_dict_list = [ reference_dictionnary ], backtrack_list = [ backtrack_dict ], variant_class = 'snp', additional_columns = None)
	for region_type in comparison_dict:
		for variant_type in comparison_dict[region_type]:
			TP = comparison_dict[region_type][variant_type]["TP"]
			FP = comparison_dict[region_type][variant_type]["FP"]
			FN = comparison_dict[region_type][variant_type]["FN"]
			if (TP+FN) > 0:
				RECALL = TP/(TP+FN)
			else:
				RECALL = 'NA'
			if (TP+FP) > 0:
				PPV=TP/(TP+FP)
			else:
				PPV = 'NA'
			if PPV != 'NA' and RECALL != 'NA':
				F1 = (PPV + RECALL)/2
			else:
				F1 = 'NA'

			if args.add_col:
				print("\t".join(args.add_col), "OUTCOME", region_type.replace(" ", "_"), variant_type, TP, FP, FN, RECALL, PPV, F1, sep = "\t")
			else:
				print("OUTCOME", region_type.replace(" ", "_"), variant_type, TP, FP, FN, RECALL, PPV, F1, sep = "\t")

