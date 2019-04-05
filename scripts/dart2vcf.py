#!/usr/bin/env python
import sys
import csv
from collections import Iterable
import datetime

dart_in = sys.argv[1]
now = datetime.datetime.now()

def flatten(mylist):
     for item in mylist:
         if isinstance(item, Iterable) and not isinstance(item, basestring):
             for x in flatten(item):
                 yield x
         else:        
             yield item

def apply_conditions(i):
	if i  == '0': return '0/0'
	elif i == '1': return '1/1'
	elif i == '2': return '0/1'
	elif i == '-': return './.'

	return i

count = 0
with open(dart_in) as csvf:
	header_line = next(csvf)
	filedate = "##fileDate=" + now.strftime("%Y-%m-%d")
	fileformat = "##fileformat=VCFv4.0"
	source = "##source=dart2vcf.py"
	vcf_head = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', header_line.split(',')[22:]]
	vcf_head_flat = list(flatten(vcf_head))
	vcfout = open("mydart.vcf","w") 
	vcfout.write(filedate + "\n" + fileformat + "\n" + source + "\n" + '\t'.join(vcf_head_flat).replace("\r", "")) #remove carriage return with replace

	for line in csv.reader(csvf):
		chrom = line[4]
		chrom_pos = line[5]
		alleleID = line[1]
		snp = line[8]
		
		ref = snp.split(":")[1].split(">")[0]
		alt = snp.split(":")[1].split(">")[1]
		empty = "."
		genotypes = line[22:]
		l=[apply_conditions(i) for i in genotypes]
		genofields = [geno + ":.:.,.:.,.,." for geno in l]
		hit = [chrom, chrom_pos, alleleID, ref, alt, empty, 'PASS', empty, 'GT:DP:AD:GL', genofields]
		hit = list(flatten(hit))
		count = count + 1
		vcfout.write('\t'.join(hit)+ "\n")
	print("Done! Wrote " + str(count) + " SNPs to VCF file 'mydart.vcf'.")


