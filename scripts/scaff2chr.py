#!/usr/bin/env python
import re
import argparse
import pandas as pd

# Covert scaffolds to chromosomes 

parser = argparse.ArgumentParser('My program')
parser.add_argument('-c', '--chrchain')
parser.add_argument('-g', '--genofile')

args = vars(parser.parse_args())
chain = args['chrchain']
geno = args['genofile']

#load dictionary file
scaff2chr_dict = {}
with open(chain,'r') as mycvs:
    for line in mycvs:
        line = line.strip()
        scaff2chr_dict[line.split(",")[0]]=line.split(",")[1]

df = pd.read_csv(geno,sep='\t',header=None)
myfile=open("test.tab",'w')
myerr=open("test.e",'w')
for i in range(len(df)):
    try:        
        my_v = scaff2chr_dict[str(df[0][i])]
        myfile.write("%s\n"%my_v)
    except:
#       print i+1 
	my_e = df[0][i]
	myerr.write("%s\n"%my_e)
