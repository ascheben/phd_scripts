#!/usr/bin/python3

# Make SNP calling pipeline config scripts
# Scripts are submitted to SLURM using runJobs.sh

import argparse
import shutil
import os

parser = argparse.ArgumentParser('My program')
parser.add_argument('-d', '--directory')
parser.add_argument('-r', '--reference')
parser.add_argument('-n', '--names')
# Optional argument if SGS AutoSNP will be used
parser.add_argument('-p', '--SGS')

args = vars(parser.parse_args())

rootDir = args['directory']
refDir = args['reference']
rawDir = rootDir+'/raw'
bamDir = rootDir+'/bam'
seqNames = args['names']
SGS = args ['SGS']
lines = [line.rstrip('\n') for line in open(seqNames)]

sabre = """
CPU=1 
MEM=8G 
HOUR=24 
OUT_DIR='{}'
MODULES=() 
PROG='cd {}; /home/ascheben/sabre/sabre pe -m 1 -f R1.fastq -r R2.fastq -b indices.list -u unknown_barcode_R1.fastq -w unknown_barcode_R2.fastq'
OPTION_ARG=()
""".format(rawDir,rawDir)
f = open('sabre.conf', 'w')
f.write(sabre)
f.close()

fastqc = """
CPU=1 
MEM=2G 
HOUR=12 
OUT_DIR='{}'
MODULES=(
fastqc) 
PROG='cd {}; fastqc -o {} -f fastq'
OPTION_ARG=(
""".format(rawDir,rawDir,rawDir)
f = open('fastqc.conf', 'w')
f.write(fastqc)
for name in lines:
  f.write ("'{}_1.fastq'\n".format(name))
  f.write ("'{}_2.fastq'\n".format(name))
f.write(")")
f.close()

trimmomatic = """
CPU=1 
MEM=8G 
HOUR=24 
OUT_DIR='{}'
MODULES=(
java
trimmomatic) 
PROG='cd {}; java -jar /home/ascheben/Trimmomatic-0.36/trimmomatic-0.36.jar'
OPTION_ARG=(
""".format(rawDir,rawDir)
f = open('trimmomatic.conf', 'w')
f.write(trimmomatic)
for name in lines:
  f.write ("'PE {}_1.fastq {}_2.fastq {}_1_pe.fastq {}_1_se.fastq {}_2_pe.fastq {}_2_se.fastq ILLUMINACLIP:/home/ascheben/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'\n".format(name, name, name, name, name, name))
f.write(")")
f.close()

soap = """
CPU=1 
MEM=10G 
HOUR=24 
OUT_DIR='{}'
MODULES=() 
PROG='cd {}; /group/pawsey0149/groupEnv/ivec/bin/soap -D {}.index -p 4 -r 0 -v 2 -m 0 -x 1000'
OPTION_ARG=(
""".format(bamDir,rawDir,refDir)
f = open('soap.conf', 'w')
f.write(soap)
for name in lines:
  f.write ("'-a {}_1_pe.fastq -b {}_2_pe.fastq -o {}/{}.soap -2 {}.fastq.unpaired'\n".format(name, name, bamDir, name, name))
f.write(")")
f.close()

if SGS == "y":
  rename_soap = """
  CPU=1 
  MEM=1G 
  HOUR=5 
  OUT_DIR='{}'
  MODULES=() 
  PROG='cd {}; sed -i'
  OPTION_ARG=(
  """.format(bamDir,bamDir,refDir)
  f = open('rename_soap0.conf', 'w')
  f.write(rename_soap)
  for name in lines:
   f.write ("'s/^/{}P1/g {}.soap'\n".format(name,name))
  f.write(")")
  f.close()
#Replace leading whitespaces with nothing
  o = open("rename_soap.conf","a")
  for line in open("rename_soap0.conf"):
   line = line.replace("  ","")
   o.write(line + "\n") 
  o.close()
#Remove blank lines
  clean_lines = []
  with open("rename_soap.conf", "r") as f:
   dirty_lines = f.readlines()
  clean_lines = [l.strip() for l in dirty_lines if l.strip()]
  with open("rename_soap.conf", "w") as f:
   f.writelines('\n'.join(clean_lines))
#Delete old conf file
  os.popen('rm rename_soap0.conf')


soap2bam = """
CPU=1 
MEM=6G 
HOUR=18 
OUT_DIR='{}'
MODULES=(samtools) 
PROG='cd {}; /group/pawsey0149/groupEnv/ivec/scripts/soap2bam.sh -c -f {}' 
OPTION_ARG=(
""".format(bamDir,bamDir,refDir)
f = open('soap2bam.conf', 'w')
f.write(soap2bam)
for name in lines:
  f.write ("'-p {}.soap'\n".format(name))
f.write(")")
f.close()

addRG = """
CPU=1 
MEM=4G 
HOUR=18 
OUT_DIR='{}'
MODULES=(java) 
PROG='cd {}; java -jar /pawsey/sles11sp4/bio-apps/binary/picard-tools/1.123/jar/AddOrReplaceReadGroups.jar' 
OPTION_ARG=(
""".format(bamDir,bamDir)
f = open('addRG.conf', 'w')
f.write(addRG)
for name in lines:
  f.write ("'I={}.clrm.bam O={}.rg.clrm.bam RGID={} RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM={}'\n".format(name, name, name, name))
f.write(")")
f.close()

addRG = """
CPU=1
MEM=10G
HOUR=10
OUT_DIR='{}'
MODULES=(samtools)
PROG='cd {}; samtools index'
OPTION_ARG=(
""".format(bamDir,bamDir)
f = open('index.conf', 'w')
f.write(addRG)
for name in lines:
  f.write ("'{}.rg.clrm.bam'\n".format(name))
f.write(")")
f.close()

markDuplicates = """
CPU=1 
MEM=60G 
HOUR=24 
OUT_DIR='{}'
MODULES=(java) 
PROG='cd {}; java -Xmx9g -jar /pawsey/sles11sp4/bio-apps/binary/picard-tools/1.123/jar/MarkDuplicates.jar REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT' 
OPTION_ARG=(
""".format(bamDir,bamDir)
f = open('markDuplicates.conf', 'w')
f.write(markDuplicates)
for name in lines:
  f.write ("'INPUT={}_sort.bam OUTPUT={}.clrm.bam METRICS_FILE={}.stats'\n".format(name, name, name))
f.write(")")
f.close()

merge_in = """
CPU=1
MEM=8G
HOUR=18
OUT_DIR='{}'
MODULES=(samtools)
PROG='cd {}; samtools merge merged.bam *rg.clrm.bam; samtools index merged.bam'
OPTION_ARG=(
)""".format(bamDir, bamDir)
f = open('merge_in.conf', 'w')
f.write(merge_in)
f.close()

call = """
CPU=4
MEM=60G
HOUR=24
OUT_DIR='{}'
MODULES=(java/8u66)
PROG='cd {}; java -jar /home/ascheben/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R /scratch/pawsey0149/ascheben/GBS_Thesis_Project/ref/Darmor_v81_assembly_RI.fasta --emitRefConfidence GVCF -ploidy 2 -variant_index_type LINEAR -variant_index_parameter 128000'
OPTION_ARG=(
""".format(bamDir,bamDir)
f = open('call.conf', 'w')
f.write(call)
for name in lines:
  f.write ("'-I {}.rg.clrm.bam -o {}.g.vcf'\n".format(name, name))
f.write(")")
f.close()

if SGS == "y":
  seqNames1 = ','.join(lines)
  SGSautoSNP = """
  CPU=12
  MEM=12G
  HOUR=48
  OUT_DIR='{}'
  MODULES=() 
  python /group/pawsey0149/groupEnv/ivec/scripts/SGSautoSNP.py --bam {}/merged.bam --cultivars {} --fasta {} --chr_offset {}/contig.gff3 --contig {}/contig_output --chr_output {}/chr_output --cpu 12 --pl myStudy --snp_id_prefix  mySNP
  OPTION_ARG=()
  """.format(bamDir, bamDir, seqNames1, refDir, refDir, bamDir, bamDir)
  f = open('SGSautoSNP0.conf', 'w')
  f.write(SGSautoSNP)
  f.close()
#Replace leading whitespaces with nothing
  o = open("SGSautoSNP.conf","a") #open for append
  for line in open("SGSautoSNP0.conf"):
   line = line.replace("  ","")
   o.write(line + "\n")
  o.close()
#Remove blank lines
  clean_lines = []
  with open("SGSautoSNP.conf", "r") as f:
   dirty_lines2 = f.readlines()
  clean_lines = [l.strip() for l in dirty_lines2 if l.strip()]
  with open("SGSautoSNP.conf", "w") as f:
   f.writelines('\n'.join(clean_lines))
#Delete old conf file
  os.popen('rm SGSautoSNP0.conf')

folders = ['analysis', 'raw', 'bam', 'vcf']

for folder in folders:

    os.mkdir(os.path.join(rootDir, folder))

raw = '{}/raw'.format(rootDir)
bam = '{}/bam'.format(rootDir)

files = os.listdir(rootDir)

for f in files:
    if (f.startswith("fastqc") or f.startswith("sabre") or f.startswith("trimmomatic")):
        shutil.move(f, raw)
    elif (f.startswith("soap") or f.startswith("mark") or f.startswith("call")or f.startswith("merge") or f.startswith("add") or f.startswith("SGS") or f.startswith("rename")):
        shutil.move(f, bam)
