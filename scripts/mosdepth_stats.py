#!/usr/bin/env python
from __future__ import division
import sys
import fileinput
# takes a single .per-base.bed mosdepth output file
# calculates mean coverage of non-zero bases
sum_covbases = 0
sum_bases = 0
for l in fileinput.input():
    l_arr=l.rstrip().split("\t")
    scaff = l_arr[0]
    start = int(l_arr[1])
    end = int(l_arr[2])
    cov = int(l_arr[3])
    if cov > 0:
        bases = end - start
        covbases = bases * cov 
        sum_covbases = sum_covbases + covbases
        sum_bases = sum_bases + bases
mean_cov= sum_covbases / sum_bases
print(str(mean_cov))
