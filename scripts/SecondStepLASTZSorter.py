#!/usr/bin/env python
import sys
sys.path.append("./bx/")
sys.path.append("./bx_extras/")
import bx.align.lav
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
""" Iterates over lav-alignment, prints for each contig highest scoring alignment """

if len(sys.argv) != 3:
    sys.stderr.write("Usage: ./SecondStepLASTZSorter.py lav-alignment.lav contig-sequences.fasta\n")
    sys.exit(1)

def get_score(alignment_list):
    return int(alignment_list[1].replace("score=",""))

def get_chr_start(alignment_list):
    return int(alignment_list[4])


# get all contigs
file_handle = open(sys.argv[1])
alignment_dict = defaultdict(list)
for alignment in bx.align.lav.Reader(file_handle):
    alignment_list = str(alignment).split()
    contig = alignment_list[10]
    try:
        contig = contig[:contig.index("[")]
    except:
        pass
    alignment_dict[contig].append(alignment_list)


# get the best alignment for each contig
best_alignment_dict = defaultdict(list)
for contig in alignment_dict:
    alignments = alignment_dict[contig]
    best_alignment = ""
    best_score = 0
    for alignment in alignments:
        alignment_score = get_score(alignment)
        if alignment_score > best_score:
            best_score = alignment_score
            best_alignment = alignment
    best_alignment_dict[contig] = best_alignment

# sort alignments
sort_dict = defaultdict(list)
for contig in best_alignment_dict:
    best_alignment = best_alignment_dict[contig]
    sort_dict[ get_chr_start(best_alignment)] = best_alignment[13] + contig

output = open( sys.argv[2].replace(".","_second_sorted."), "a")
sequences = SeqIO.index(sys.argv[2], "fasta")
for position in sorted(sort_dict.keys()):
    nice_name = sort_dict[position][1:]
    if sort_dict[position][0] == "+":
        output.write(">%s\n%s\n" %(nice_name,str(sequences[nice_name].seq)))
    else:
        output.write(">%s\n%s\n" %(nice_name,str(Seq(str(sequences[nice_name].seq)).reverse_complement())))
