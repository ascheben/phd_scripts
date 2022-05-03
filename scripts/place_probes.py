import sys
import csv
from collections import defaultdict

#note: generate input file with: blastn -query snp_probes.fa -db reference.fa -num_threads 4 -dust no -task blastn-short -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" -word_size 5 -out myseq.blastn
#note: blastn output coordinates are 1-based NOT 0-based
#note: query start coordinate is always < query end coordinate
#note: hit start coordinate can be > hit end coordinate (reverse hit)

blastout = sys.argv[1]
snpdict = defaultdict(list)
def get_positions(probe_name, blast_list):
   """
   Calculate the position of a SNP based on a blast result. This requires qseq and sseq input to account
   for gaps when identifying the precise coordinates. The SNP position on the query is extracted from
   the probe name.
   """
   query_s = int(blast_list[1])
   query_e = int(blast_list[2])
   source_s = int(blast_list[3])
   source_e = int(blast_list[4])
   qseq = blast_list[6]
   sseq = blast_list[7]
   # extract from name
   snp_position = int(probe_name.split("_")[-1])
   # 1) Check whether SNP site is aligned, if not return "NA"
   if snp_position >= query_s and snp_position <= query_e:
       pass
   else:
       print("SNP position out of query bounds ", probe_name, blast_list)
       return("NA","NA")
   # 2) Is source sequence reversed
   rev = False
   if source_s > source_e:
       rev = True
   position = query_s
   for idx,nuc in enumerate(qseq):
       if position == snp_position:
           snp_index = idx
           break   
       if nuc != "-":
           position += 1
   if snp_index:
       seq = sseq[:snp_index].replace("-", "")
       seqlen = len(seq)
       if rev:
           strand = "-"
           source_snp_pos = source_s - seqlen
           return(source_snp_pos,strand)
           #print("Hit was reverse, source SNP position = ", source_snp_pos)
       else:
           strand = "+"
           source_snp_pos = source_s + seqlen
           return(source_snp_pos,strand)
           #print("Hit was on same strand, source SNP position = ", source_snp_pos)
   else:
       print("No SNP index found ", probe_name, blast_list)
       return("NA","NA")

# Input file example line for blastout
#Bn-A01-p19602112__1based_SNP_61	chrA01_contigs_placed_v81	100.000	121	0	0	1	121	20394744	20394864	3.83e-62	240	AACGTTTGTCACATTCTATTTCATATGAAAAAATTGGTATACGAAGTAATATTTCCTATTTCAATACGCAAATCGTTTGTATATATATGAAGTAATATTATTCGTTCTAATGAGGATTTTC	AACGTTTGTCACATTCTATTTCATATGAAAAAATTGGTATACGAAGTAATATTTCCTATTTCAATACGCAAATCGTTTGTATATATATGAAGTAATATTATTCGTTCTAATGAGGATTTTC
with open(blastout) as tsv2:
    for line in csv.reader(tsv2, dialect="excel-tab"):
       #The output was filtered to retain only hits with at least 95% identity
       if float(line[2]) >= 90:
           snpname = line[0]
           chrom = line[1]
           query_s = line[6]
           query_e = line[7]
           hit_s = line[8]
           hit_e = line[9]
           bit = line[11]
           qseq = line[12]
           sseq = line[13]
           hit = [chrom, query_s, query_e, hit_s, hit_e, bit, qseq, sseq]
           snpdict[snpname].append(hit)

snp_positions = []
ambiguous_snps = []
for snp, value in snpdict.items(): #loop over snps and their BLAST hits
       # sort by bit score
       value = sorted(value, key=lambda x: float(x[5]))
       value.reverse()
       if len(value) == 1:
               chrom = value[0][0]
               real_pos = get_positions(snp,value[0])
               snp_positions.append([str(snp),chrom,str(real_pos[0]),real_pos[1]])
       elif len(value) > 1:
               if float(value[0][5]) > float(value[1][5]):
                       chrom =  value[0][0]
                       real_pos = get_positions(snp,value[0])
                       snp_positions.append([str(snp),chrom,str(real_pos[0]),real_pos[1]])
               else:
                       ambiguous_snps.append(snp)
       else:
               print("Warning: " + snp + " has no blast hits")

with open('snp_positions.txt', 'w') as file:
    file.writelines('\t'.join(i) + '\n' for i in snp_positions)
with open('ambiguous_snps.txt', 'w') as file:
       for snp in ambiguous_snps:
               file.writelines(snp + "\n")
