import sys
import os
import argparse
from argparse import RawTextHelpFormatter
import subprocess
import shlex
import csv
from itertools import izip


def parse_user_args():
    """
    This function parses the arguments provided by the user
    :return: a dictionary having a key for each arguments
    :rtype: dict
    """
    parser = argparse.ArgumentParser(
        description='Welcome to VCF2MST.\n'
                    'Before executing this script, Java/1.8.0_45 or higher must be loaded using module load ',
        usage='python vcf2mst.py -i [INPUT_VCF] -r [REFERENCE] -g [GATK_PATH]',
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
    run_gatk.communicate()


var2tab()
arg = parse_user_args()
tab = os.path.splitext(arg['input'])[0] + ".tab"
mstout = os.path.splitext(arg['input'])[0] + ".mst"
rqtltemp = os.path.splitext(arg['input'])[0] + ".tmp"
rqtlout = os.path.splitext(arg['input'])[0] + ".rqtl.csv"

with open(tab, 'r') as tab:
  header_line = next(tab)
  header_line = header_line.split()
  #print("locus_name" + "\t" + '\t'.join(header_line[4:]))
  with open(mstout, "w") as f:
    print >> f, "locus_name" + "\t" + '\t'.join(header_line[4:])
  with open(rqtltemp, "w") as f:
    print >> f, "Sample" + ",," + ','.join(header_line[4:])

  for ourline in tab:
    ourline = ourline.split()
    ageno = set(ourline[2].split('/'))
    bgeno = set(ourline[3].split('/'))
      
    if ageno == set(['.']) or bgeno == set(['.']):
      pass #continue
    
    if len(ageno) > 1 or len(bgeno) > 1:
      pass # continue
    
    a = ageno.pop()
    b = bgeno.pop()
    name = ourline[0] + '_' + ourline[1]
    qtlname = ourline[1] + ',' + ourline[0].split('_')[0] 

    qtl_ll = [qtlname]
    new_ll = [name]
    for element in ourline[4:]:
      transl_element = element.replace(b, 'B').replace(a, 'A')
      qtl_element = element.replace(b, 'B').replace(a, 'A')
      if transl_element == 'A/A':
        transl_element = 'A'
      elif transl_element == 'B/B':
        transl_element = 'B'
      elif transl_element == 'B/A':
        transl_element = 'X'

      elif transl_element == 'A/B':
        transl_element = 'X'
      else:
        transl_element = 'U'
      new_ll.append(transl_element)

      if qtl_element == 'A/A':
        qtl_element = 'A'
      elif qtl_element == 'B/B':
        qtl_element = 'B'
      elif qtl_element == 'B/A':
        qtl_element = 'H'
      elif qtl_element == 'A/B':
        qtl_element = 'H'
      else:
        qtl_element = '-'
      qtl_ll.append(qtl_element)

    #print '\t'.join(new_ll)
    with open(mstout, "a") as f:
      print >> f, '\t'.join(new_ll)

    with open(rqtltemp, "a") as f:
      print >> f, ','.join(qtl_ll)
a = izip(*csv.reader(open(rqtltemp, "rb")))
csv.writer(open(rqtlout, "wb")).writerows(a)
os.remove(rqtltemp)


