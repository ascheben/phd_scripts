#!/usr/bin/env python
"""
.. module:: run
   :synopsis: Convert VCF to flapjack.

This is the VCF2FJ converter.

To get help, ``python vcf2fj.py -h``.

"""

import sys
import os
import argparse
from argparse import RawTextHelpFormatter
import subprocess
import shlex

def parse_user_args():
    """
    This function parses the arguments provided by the user
    :return: a dictionary having a key for each arguments
    :rtype: dict
    """
    parser = argparse.ArgumentParser(
        description='Welcome to VCF2FJ.\n'
                    'Before executing this script, Java/1.8.0_45 or higher must be loaded using module load ',
        usage='python vcf2fj.py -i [INPUT_VCF] -r [REFERENCE] -g [GATK_PATH]',
        formatter_class=RawTextHelpFormatter, add_help=False)

    parser.add_argument(
        '-i', '--input', dest='input', required=True, metavar='VCF FILE', help='Input variant file in VCF format. ')
    parser.add_argument(
        '-r', '--reference', dest='reference', required=True, metavar='FASTA FILE', help='Input variant file in FASTA format, muste be faidx indexed and have a corresponding picard dictionary file. ')
    parser.add_argument(
        '-g', '--gatk_path', dest='gatk_path', required=True, metavar='PATH', help='Path to GATK jar file. ')
    parser.add_argument('-h', '--help', action="help", help="Show this help message and exit")


    return vars(parser.parse_args())

def var2tab(): 
    arg = parse_user_args()
    infile = arg['input']
    ref = arg['reference']
    gatk_path = arg['gatk_path']
    outfile = os.path.splitext(infile)[0] + ".tab"

    gatk_cmd = shlex.split("java -jar %s -R %s -T VariantsToTable -V %s -F CHROM -F POS -F REF -F ALT -GF GT -o %s" % (gatk_path, ref, infile, outfile))
    run_gatk = subprocess.Popen(gatk_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #stdout, stderr =  _map.communicate()

def tab2geno(): 
    arg = parse_user_args()
    tab = os.path.splitext(arg['input'])[0] + ".tab"
    outgeno = tab[:-3] + "fj.geno" 

    for ourline in open(tab):
        ourline = ourline.split()
        
        if ourline[0] == "CHROM":
            name = ""
        else:
            name = ourline[0] + '_' + ourline[1]
      
        new_ll = [name]
        for element in ourline[4:]:
            if element == 'A/A':
                element = 'A'
            elif element == 'C/C':
                element = 'C'
            elif element == 'G/G':
                element = 'G'
            elif element == 'T/T':
                element = 'T'
            elif element == './.':
                element = '-'
            else:
                pass
            new_ll.append(element)
      
        #print '\t'.join(new_ll)
        with open(outgeno, "a") as f:
            print >> f, '\t'.join(new_ll) 

def tab2map():
    arg = parse_user_args()
    tab = os.path.splitext(arg['input'])[0] + ".tab"
    outmap = tab[:-3] + "fj.map"

    for ourline in open(tab):
        ourline = ourline.split()
        if ourline[0] != "CHROM":
            name = ourline[0] + '_' + ourline[1] + '\t' + ourline[0] + '\t' + ourline[1]
            new_ll = [name]
          
            with open(outmap, "a") as f:
                print >> f, '\t'.join(new_ll)

    os.remove(tab)

def main():
    var2tab()
    tab2geno()
    tab2map()

if __name__ == "__main__":
    main()

