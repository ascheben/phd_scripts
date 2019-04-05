#!/usr/bin/env python3
import os
import sys

# Outputs command for plotting with R/VennDiagram

set_a_name = sys.argv[1]
set_b_name = sys.argv[2]
set_c_name = sys.argv[3]
foursets = 0

if len(sys.argv) == 5:
    set_d_name = sys.argv[4]
    foursets = 1

set_a = set(line.strip() for line in open(set_a_name))
set_b = set(line.strip() for line in open(set_b_name))
set_c = set(line.strip() for line in open(set_c_name))

len_a = len(set_a)
len_b = len(set_b)
len_c = len(set_c)

nab = set_a.intersection(set_b) 
nac = set_a.intersection(set_c) 
nbc = set_b.intersection(set_c) 
nabc = set_a.intersection(set_b, set_c) 

if foursets == 1:
    set_d = set(line.strip() for line in open(set_d_name))
    len_d = len(set_d)
    nad  = set_a.intersection(set_d)
    nbd = set_b.intersection(set_d)
    ncd = set_c.intersection(set_d)
    nabd = set_a.intersection(set_b, set_d)
    nacd = set_a.intersection(set_c, set_d)
    nbcd = set_b.intersection(set_c, set_d)
    nabcd = set_a.intersection(set_b, set_c, set_d)

if foursets == 0:
    print('draw.triple.venn(area1 = %s, area2 = %s, area3 = %s, n12 = %s, n23 = %s, n13 = %s, n123 = %s, category = c("%s", "%s", "%s"), lty = "blank", fill = c("skyblue", "pink1", "mediumorchid"))' % (len_a, len_b, len_c, len(nab), len(nbc), len(nac), len(nabc), set_a_name, set_b_name, set_c_name)) 

else:
    print('draw.quad.venn(area1 = %s, area2 = %s, area3 = %s, area4 = %s, n12 = %s, n13 = %s, n14 = %s, n23 = %s, n24 = %s, n34 = %s, n123 = %s, n124 = %s, n134 = %s, n234 = %s, n1234 = %s, category = c("%s", "%s", "%s", "%s"), lty = "blank", fill = c("skyblue", "pink1", "mediumorchid", "red"))' % (len_a, len_b, len_c, len_d, len(nab), len(nac), len(nad), len(nbc), len(nbd), len(ncd), len(nabc), len(nabd), len(nacd), len(nbcd), len(nabcd), set_a_name, set_b_name, set_c_name, set_d_name))

