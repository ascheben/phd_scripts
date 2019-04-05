#!/usr/bin/env bash
# Covert custom geno table to MSTMap input
# Input genotype table
invcf="$1"
# P-value threshold for linkage groups
pval="$2"

ind=`head -1 ${invcf} | cut -d, -f2- | tr ',' '\n' | wc -l`
loci=`tail -n +2 ${invcf} | wc -l`

read -d '' params << EOF 
population_type RIL2
population_name Pop
distance_function kosambi
cut_off_p_value ${pval}
no_map_dist 20.0
no_map_size 2
missing_threshold 0.2
estimation_before_clustering no
detect_bad_data yes
objective_function ML
number_of_loci ${loci}
number_of_individual ${ind}
EOF

echo "${params}" >> ${invcf%%gt}mst_ml

sed 's/H/X/g' ${invcf} | tr ',' '\t' >> ${invcf%%gt}mst_ml

