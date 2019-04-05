#!/usr/bin/env bash
# Obtain cumulative bases covered by BLAST hits
# Can e.g. help identify contaminant scaffolds

cut -f1 myseq.blastn.tsv | sort | uniq > myseq.blastn.contig.list
while read l <&3; do
awk -v contig=$l '$1==contig' myseq.blastn.tsv >> $l.contig.temp
  while read q <&4; do
  qstart=`echo $q | awk '{print $10}'`
  qend=`echo $q | awk '{print $11}'`
  seq $qstart $qend | tr "[[:blank:]]" "\n" >> $l.seq.temp
  done 4<$l.contig.temp
unique_bases=`sort ${l}.seq.temp | uniq | wc -l`
echo "$l : $unique_bases" >> $l.contig.cumulative.bases
done 3<myseq.blastn.contig.list
