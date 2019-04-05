#!/usr/bin/env bash

# Removes paired reads shorter than min len
# User args: FASTQ_R1 FASTQ_R2 MINLEN_R1 MINLEN_R2
R1fq=$1
R2fq=$2
minlen1=$3
minlen2=$4

fix_base_count() {
    local counts=($(cat))
    echo "${counts[0]} $((${counts[1]} - ${counts[0]}))"
}

if [ ${R1fq: -3} == ".gz" ]
then
	gzip -d $R1fq
	R1fq=${R1fq%%.gz}
fi

if [ ${R2fq: -3} == ".gz" ]
then
	gzip -d $R2fq
	R2fq=${R2fq%%.gz}
fi

R1before=`awk 'NR % 4 == 2' $R1fq | wc -cl | fix_base_count`
R2before=`awk 'NR % 4 == 2' $R2fq | wc -cl | fix_base_count`
R1before_nreads=`echo $R1before | cut -d' ' -f1`
R2before_nreads=`echo $R2before | cut -d' ' -f1`
R1before_nbases=`echo $R1before | cut -d' ' -f2`
R2before_nbases=`echo $R2before | cut -d' ' -f2`
R1before_mlen=`echo "$R1before_nbases / $R1before_nreads" | bc -l`
R2before_mlen=`echo "$R2before_nbases / $R2before_nreads" | bc -l`


#2. Find all entries with read length less than minimum length and print line numbers, for both R1 and R2
awk -v min=$minlen1 '{if(NR%4==2) if(length($0)<min) print NR"\n"NR-1"\n"NR+1"\n"NR+2}' $R1fq > ${R1fq}temp.lines1
awk -v min=$minlen2 '{if(NR%4==2) if(length($0)<min) print NR"\n"NR-1"\n"NR+1"\n"NR+2}' $R2fq >> ${R1fq}temp.lines1

#3. Combine both line files into one, sort them numerically, and collapse redundant entries
sort -n ${R1fq}temp.lines1 | uniq > ${R1fq}temp.lines
rm ${R1fq}temp.lines1

#4. Remove the line numbers recorded in "lines" from both fastqs
awk 'NR==FNR{l[$0];next;} !(FNR in l)' ${R1fq}temp.lines $R1fq > ${R1fq%%fq}${minlen1}.fq
awk 'NR==FNR{l[$0];next;} !(FNR in l)' ${R1fq}temp.lines $R2fq > ${R2fq%%fq}${minlen2}.fq
rm ${R1fq}temp.lines

R1after=`awk 'NR % 4 == 2' ${R1fq%%fq}${minlen1}.fq | wc -cl | fix_base_count`
R2after=`awk 'NR % 4 == 2' ${R2fq%%fq}${minlen2}.fq | wc -cl | fix_base_count`

R1after_nreads=`echo $R1after | cut -d' ' -f1`
R2after_nreads=`echo $R2after | cut -d' ' -f1`
R1after_nbases=`echo $R1after | cut -d' ' -f2`
R2after_nbases=`echo $R2after | cut -d' ' -f2`
R1after_mlen=`echo "$R1after_nbases / $R1after_nreads" | bc -l`
R2after_mlen=`echo "$R2after_nbases / $R2after_nreads" | bc -l`
nreads_del=`echo "$R1before_nreads - $R1after_nreads" | bc -l`

gzip $R1fq $R2fq ${R1fq%%fq}${minlen1}.fq ${R2fq%%fq}${minlen2}.fq


if [ $R1after_nreads -eq $R2after_nreads ]
then
	printf "Done! %s reads dropped. R1 mean len before: %s, after: %s. R2 mean len before: %s, after: %s. R1 (%s) and R2 (%s) read numbers match after filtering!\n" "$nreads_del" "$R1before_mlen" "$R1after_mlen" "$R2before_mlen" "$R2after_mlen" "$R1after_nreads" "$R2after_nreads" 
else
	printf "Done! %s reads dropped. R1 mean len before: %s, after: %s. R2 mean len before: %s, after: %s. R1 (%s) and R2 (%s) read numbers do not match after filtering!\n" "$nreads_del" "$R1before_mlen" "$R1after_mlen" "$R2before_mlen" "$R2after_mlen" "$R1after_nreads" "$R2after_nreads" 
fi
