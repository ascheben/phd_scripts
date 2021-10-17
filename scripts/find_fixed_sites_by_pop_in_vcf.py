#!/usr/bin/env python3

import sys
import csv
from collections import defaultdict
from collections import Counter
import itertools
from argparse import ArgumentParser

parser = ArgumentParser(description='''Finds fixed heterozygous alleles within and between populations
                        from a vcf and a popmap file''')
parser.add_argument('-v', '--version', action='version', version='1.0')
parser.add_argument('-i',
                    dest='vcf',
                    help='the input vcf file')
parser.add_argument('-p',
                    dest='popmap',
                    help='Population map file with two tab-separated columns: population\tsample name')
parser.add_argument('-m',
                    dest='minnonmiss',
                    help='''Absolute minimum number of non-missing genotype calls
                    in a population to analyse a site''')
parser.add_argument('-f',
                    dest='minfrequency',
                    help='Fixation frequency cutoff of target allele for a population (0-1)')
args = parser.parse_args()

if None not in [args.vcf, args.popmap]:
    if args.minnonmiss == None:
        mcutoff = 1
    else:
        mcutoff = args.minnonmiss
    if args.minfrequency == None:
        fcutoff = 1
    else:
        fcutoff = args.minfrequency
else:
    print
    parser.print_help()
    print
    sys.exit(1)


# could be switch for 1/1 or 0/0 or ./.
# assume correct VCF format without "1/0" encoding of het
target_geno = "0/1"

vcf = args.vcf #input vcf file
pop = args.popmap #input popmap file

print("Starting analysis of population-based fixed sites...")
print("Using fixation frequency cutoff: ", fcutoff)
print("Minimum non-missing samples per site per population: ", mcutoff)



def combinations(mylist):
    slist = []
    for L in range(0, len(mylist)+1):
        for subset in itertools.combinations(mylist,L):
            subset = list(subset)
            if len(subset)>0:
                subset = sorted(subset)
                subset = "_".join(subset)
                slist.append(subset)
    return slist

## PARSE POPMAP ##
# make a pop dict
pops = defaultdict(str)
plist = []
with open(pop,'r') as popmap:
    for line in popmap:
        line = line.strip()
        line = line.split("\t")
        pops[line[1]] = line[0]
        if line[0] not in plist:
            plist.append(line[0])

print("Site of interest: ", target_geno)
print("Populations being analysed: ", plist)


## MAKE DICTS ##
comb_dict = defaultdict(int)
comb_snps = defaultdict(list)
comb_pops = combinations(plist)
for cp in comb_pops:
    comb_dict[cp] = 0
    comb_snps[cp] = []
names = []
# number of non-genotype fields in standard VCF format
num_info_fields = 9

# genotypes
gdict = defaultdict(lambda: defaultdict(list))
# target allele frequencies
fdict = defaultdict(lambda: defaultdict(float))


## PARSE VCF ##
with open(vcf) as tsv:
    for line in csv.reader(tsv, dialect="excel-tab"):
        # skip comment lines
        if line[0].startswith("#CHROM"):
            names = line[num_info_fields:]
        elif not line[0].startswith("#"):
            # unique SNP chrom+pos
            chrompos = line[0] + "__" + line[1]
            # get number of samples
            num_samples = len(line) - num_info_fields
            # loop through alleles for each SNP and count genotypes
            for index,allele in enumerate(line[num_info_fields:]):
                # sample
                name = names[index]
                if name in pops:
                    population = pops[name]
                else:
                    population = "Undefined"
                #genotype
                g = allele.split(":")[0]
                gdict[chrompos][population].append(g)
            for p in gdict[chrompos]:
                p_alleles = gdict[chrompos][p]
                tot_alleles = len(p_alleles)
                p_counts = Counter(p_alleles)
                missing = p_counts["./."]
                tot_nonmissing = tot_alleles - missing
                #print(p,tot_alleles,p_counts,tot_nonmissing)
                if tot_nonmissing >= mcutoff:
                    prop_target = p_counts[target_geno] / tot_nonmissing
                else:
                    prop_target = None
                fdict[chrompos][p] = prop_target
            fixed_pops = []
            for key,value in fdict[chrompos].items():
                if value and value >= fcutoff:
                    fixed_pops.append(key)
                    #print("fixed_pops", fixed_pops)
            if len(fixed_pops) > 0:
                fpop_comb = combinations(fixed_pops)
                for fpc in fpop_comb:
                    if fpc in comb_dict:
                        comb_dict[fpc] += 1
                        comb_snps[fpc].append(chrompos)

## OUTPUT ##

print("Fixed target sites per population: ")
for k,v in comb_dict.items():
    print(k,v)
print("\n")
print("Writing output files...")
with open("populations_shared_counts.txt", "w") as counts:
    header = "pop_group" + "\t" + "shared_fixed_sites" + "_" + target_geno + "\n"
    counts.write(header)
with open("populations_shared_counts.txt", "a") as counts:
    for k,v in comb_dict.items():
        out_comb = k + "\t" + str(v) + "\n"
        counts.write(out_comb)
with open("populations_fixed_sites.txt", "w") as sites:
    header = "pop" + "\t" + "chromosome" + "\t" + "position" + "\n"
    sites.write(header)
with open("populations_fixed_sites.txt", "a") as sites:
    for k,v in comb_snps.items():
        for pos in v:
            chromosome = pos.split("__")[0]
            position = pos.split("__")[1]
            outsnp = k + "\t" + chromosome + "\t" + position + "\n"
            sites.write(outsnp)
# make output file with col1 pop, col2 chrom, col3 pos
with open("populations_frequencies.txt", "w") as freqs:
    header = "chrom" + "\t" + "position" + "\t" + "pop" + "\t" + "frequency" + "_" + target_geno + "\n"
    freqs.write(header)
with open("populations_frequencies.txt", "a") as freqs:
    for k,v in fdict.items():
        for pop,f in v.items():
            if f is not None:
                fround = round(f,5)
            else:
                fround = "NA"
            chromosome = k.split("__")[0]
            position = k.split("__")[1]
            outfreq = chromosome + "\t" + position + "\t" + pop + "\t" + str(fround) + "\n"
            freqs.write(outfreq)
