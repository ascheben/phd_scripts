#!/usr/bin/env python
from __future__ import division
import argparse
from cyvcf2 import VCF
import csv
from scipy.stats import chisquare

# Outputs VCF as filtered geno table for linkage mapping
# Genotypes are converted to parental genotype (A|B|H|-)
parser = argparse.ArgumentParser('My program')
# Provide list of parental individuals
parser.add_argument('-p1', '--mother')
parser.add_argument('-p2', '--father')
# Provide VCF
parser.add_argument('-g', '--genofile')
 

args = vars(parser.parse_args())

p1 = args['mother'].split(',')
p2 = args['father'].split(',')
geno = args['genofile']

vcf = VCF(geno)
DAD = []
MUM = []

for sample in p1:
	DAD.append(vcf.samples.index(sample))
for sample in p2:
	MUM.append(vcf.samples.index(sample))
PAR = p1 + p2
progeny = [e for e in vcf.samples if e not in PAR]
pindex = [] #progeny indices
for pro in progeny:
	pindex.append(vcf.samples.index(pro))
cols = ["ID"] #first column of header

header = cols + progeny #append ID and parent names to start of header

if geno.endswith('.vcf'):
	gtfile = geno[:-4] + ".gt"
else:
	gtfile = geno + ".gt"
# open output file
with open(gtfile, 'w') as gtout:
	writer = csv.writer(gtout)
	writer.writerow(header)
	for var in vcf: 
		var_id = var.CHROM + "__" + str(var.end)
		var_id.encode('utf8', 'replace')
		dad_gt = []
		mum_gt = []
		for sindex in DAD:
			dad_gt.append(var.gt_types[sindex])
		for sindex in MUM:
			mum_gt.append(var.gt_types[sindex])
		gt = []
		
		# if MUM is HOM_REF and DAD is HOM_ALT or vice versa
                # gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
		# Complicated bool checks to ensure that genotypes can be missing
		if (any(v == 0 for v in mum_gt) and all(v == 0 or v == 2 for v in mum_gt) and any(v == 3 for v in dad_gt) and all(v == 3 or v == 2 for v in dad_gt)) or (any(v == 3 for v in mum_gt) and all(v == 3 or v == 2 for v in mum_gt) and any(v == 0 for v in dad_gt) and all(v == 0 or v == 2 for v in dad_gt)): 
		#if (mum_gt == 0 and dad_gt == 3) or (mum_gt == 3 and dad_gt == 0):
			# Append all mum and then dad genotypes
                        # These should be excluded for mapping
			#gt.extend(dad_gt + mum_gt)
	
			#Get consensus mum allele
			if (any(v == 0 for v in mum_gt)):
				con_mum = 0
			elif (any(v == 3 for v in mum_gt)): 
				con_mum = 3
			else:
				print("Error! Mum is not REF or ALT")
				print(mum_gt)
			if (any(v == 0 for v in dad_gt)):
				con_dad = 0
			elif (any(v == 3 for v in dad_gt)): 
				con_dad = 3
			else:
				print("Error! Dad is not REF or ALT")
				print(dad_gt)
			print(var)
			print("Mum (A): " + str(mum_gt))
			print("Dad (B): " + str(dad_gt))
			for idx in pindex:
				gt.append(var.gt_types[idx])
			for n, i in enumerate(gt): #change to AB format
				if i == con_mum:
					gt[n] = "A"
				if i == con_dad:
					gt[n] = "B"
				if i == 2:
					gt[n] = "-"
				if i == 1:
					gt[n] = "H"
				print(str(i) + " assigned to " + str(gt[n]))

			missing_prop = gt.count('-') / len(gt)
			total_progeny = len(gt) - gt.count('-')
			exp_a = total_progeny / 4
			exp_b = exp_a
			exp_ab = exp_a + exp_b
			obs_a = gt.count('A')
			obs_b = gt.count('B')
			obs_ab = gt.count('H')
			chisq_p = chisquare([obs_a, obs_ab, obs_b], [exp_a, exp_ab, exp_b]).pvalue
			if missing_prop < 0.2 and chisq_p > 0.01:
				genotypes = gt.insert(0,var_id)
				writer.writerow(gt)

