#!/usr/bin/env python
# Author: Philipp Bayer (UWA)

""" This script parses a tab-delimited (-m 8) BLASTn-output file and assigns the best hit (chromosome) to each query. 
It then checks the best two hits for each query and discards the query if the hits are too similar (+/- 1%) by writing the query
into the file "Strange_seqs.fasta".
For each hit (chromosome), one fasta-file based on the sequence-name is generated which contains the fasta-sequences of the queries
that fit this hit (chromosome) the best. """

import sys
if len(sys.argv) != 3:
    sys.stderr.write("This script parses ONE tab-delimited (m 8) output from BLAST.\nFor each query, it compares the best two hits and writes those that have too similar scores to Strange_seqs.fasta.\nThen, the script writes the query-sequences for the best subject-hit into a fasta-file named after the subject.\n")
    sys.stderr.write("Usage is: ./MakeFASTAForLASTZ.py blast-output.blastn all_queries.fasta\n")
    sys.exit(1)

from Bio import SeqIO
toparse = open(sys.argv[1])
PERCENTAGE = 0.01
# go through blastn-file
# get the first two hits for each files (first 2 rows)
# when the scores are +/- 1% similar, write to strange sequences
# else write to chromosome-file
strange_seq = open("Strange_seqs.fasta","w")
strange_seq_expl = open("Strange_seqs_annotation.csv","w")
all_sequences = SeqIO.index(sys.argv[2], "fasta")

seq_dict = {}
sys.stderr.write("Started parsing the blastn-file.\n")

for line in toparse:
    line_list = line.rstrip().split("\t")
    query = line_list[0]
    chromosome = line_list[1]
    score = float(line_list[11])
    try:
        if len(seq_dict[query]) < 2 and chromosome != seq_dict[query][0][1]:
            seq_dict[query].append(line_list)
    except:
        seq_dict[query] = [line_list]

sys.stderr.write("Finished parsing. Now kicking out the 'weird' ones.\n")
# iterate over teh dictionary, kick out keys where the scores of 
# both alignments are too similar
copy_dict = dict(seq_dict)
for key in seq_dict:
    first_line = seq_dict[key][0]
    if len(seq_dict[key]) == 1:
         continue
    second_line = seq_dict[key][1]
    print first_line
    print second_line
    first_score = float(first_line[11])
    second_score = float(second_line[11])
    lower_bound = first_score - (first_score*PERCENTAGE)
    upper_bound = first_score + (first_score*PERCENTAGE)
    # when the score of the second alignment is inside the boundaries,
    # remove
    if second_score > lower_bound and second_score < upper_bound:
        sys.stderr.write(key + " has too similar alignments. Writing out.\n")
        strange_seq.write(all_sequences[key].format("fasta"))
        strange_seq_expl.write("\t".join(first_line) + "\n" + "\t".join(second_line) + "\n")
        del copy_dict[key]


sys.stderr.write("Kicked out strange alignments. Now writing chromosome-files.\n")

file_handle_dict = {}

for key in copy_dict:
    chromosome = seq_dict[key][0][1]
    try:
        to_write = file_handle_dict[chromosome]
    except:
        to_write = open(chromosome + ".fasta", "a")
        file_handle_dict[chromosome] = to_write

    to_write.write(">" + key + "\n" + str(all_sequences[key].seq)+ "\n")

sys.stderr.write("Done! Have a nice day.\n")
