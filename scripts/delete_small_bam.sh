#!/usr/bin/env bash
#get average bam file size
avg=$(ls -l bam/*.bam | awk -v N=5 '{ sum += $N } END { if (NR > 0) printf "%.2f", sum / NR }')
#get 10% of average file size
min=$(echo $avg*0.1 | bc)
#convert to integer
num=$(echo "($min+0.5)/1" | bc)
# find (and delete) bam files with size lower than 10% of average
find . -name "*rg.*.bam" -size -${num}c -delete
