#!/usr/bin/env bash
# Assemble pangenome for small genomes
#module load jellyfish masurca bowtie2/2.2.6 samtools java/8u66
bbmap="/path/tp/bbmap"
ref="/path/to/ref/ref.fa"
outmapdir="/path/to/output"
config="/path/to/config/config.txt"
cd $outmapdir
while read base; do
bowtie2-build $ref $ref
bowtie2 -X 1500 --end-to-end --sensitive --threads 16 --mm --rg-id ${base} --rg "SM:LM_${base}" -x $ref -1 ${base}/${base}_R1.fastq -2 ${base}/${base}_R2.fastq | samtools sort -@ 16 -o ${base}/${base}_sorted.bam -O BAM --reference $ref -
samtools fastq -f 4 -1 ${outmapdir}/${base}/${base}_unmapped_R1.fastq -2 ${outmapdir}/${base}/${base}_unmapped_R2.fastq ${outmapdir}/${base}/${base}_sorted.bam
$bbmap/repair.sh in1=${outmapdir}/${base}/${base}_unmapped_R1.fastq in2=${outmapdir}/${base}/${base}_unmapped_R2.fastq \
out1=${outmapdir}/${base}/${base}_unmap_R1.fastq out2=${outmapdir}/${base}/${base}_unmap_R2.fastq \
outsingle=${outmapdir}/${base}/${base}_unmap_SE.fastq
rm -f ${outmapdir}/${base}/${base}_unmapped_R1.fastq ${outmapdir}/${base}/${base}_unmapped_R2.fastq
cp $config ${base}/config.txt
echo "PE= pe 276 127 /scratch/pawsey0149/ascheben/blackleg_genes/raw/assembly/${base}/${base}_unmap_R1.fastq /scratch/pawsey0149/ascheben/blackleg_genes/raw/assembly/${base}/${base}_unmap_R2.fastq" >> ${base}/config.txt
echo "PE= se 276 127 /scratch/pawsey0149/ascheben/blackleg_genes/raw/assembly/${base}/${base}_unmap_SE.fastq" >> ${base}/config.txt
echo "END" >> ${base}/config.txt
masurca ${base}/config.txt
mv assemble.sh ${outmapdir}/${base}
cd ${outmapdir}/${base}
bash assemble.sh
cd ${outmapdir}
$bbmap/rename.sh in=${outmapdir}/${base}/CA/10-gapclose/genome.scf.fasta out=${outmapdir}/${base}/CA/10-gapclose/genome.scf.renamed.fasta prefix=${base} && mv ${outmapdir}/${base}/CA/10-gapclose/genome.scf.renamed.fasta ${outmapdir}/${base}/CA/10-gapclose/genome.scf.fasta
cat ${outmapdir}/${base}/CA/10-gapclose/genome.scf.fasta >> $ref
done<folder.list
