#!/usr/bin/env python
#script to rename scaffolds to chromosomes in a vcf
#inputs are headerless vcf and gff3 files
#output is renamed/renumbered col1 and col2 of the input vcf
import re
import pandas as pd
import argparse
parser = argparse.ArgumentParser('My program')
parser.add_argument('-i', '--vcf')
parser.add_argument('-g', '--gff3')
parser.add_argument('-c', '--chr')
args = vars(parser.parse_args())
vcf = args['vcf']
gff3 = args['gff3']
mychr = args['chr']
#load dictionary file
#dict will contain each scaffold as key and the start coordinate as value
gff3_dict = {}
with open(gff3,'r') as mycvs:
    for line in mycvs:
        line = line.strip()
        start = int(line.split("\t")[3])
        name = line.split("\t")[8]
        scaff = name.split(";")[1][5:]
        gff3_dict[scaff]=start
#remove vcf header
df = pd.read_csv(vcf,sep='\t',header=None)
myfile=open("test.tab",'w')
myerr=open("test.e",'w')
#for snp in df:
for i in df.iterrows():
    snp_scaff = i[1][0]
    if snp_scaff in gff3_dict:
        my_v = gff3_dict[snp_scaff]
        cor_pos = int(my_v) + int(i[1][1])
        print(mychr + "\t" + str(cor_pos))
    else:
        print(i[1][0] + "\t" + str(i[1][1]))
