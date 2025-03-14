#!/bin/bash

#this script will align reads to your bowtie2 INDEXED reference then extract the reads and put into seperate files per sample.

reference=all_deep6_virus_seqs_all_progs_1_24_25.okay.mrna
readsdir="/work/alex/florida/RNA/fastp"
sampleR1list=sampleR1.list

echo "$reference ; $readsdir ; $sampleR1list"

for x in $(cat $sampleR1list);
do	echo "$x"
	name="$(echo $x | awk -F ".clean.R1.fq.gz" '{print $1}')"
	R2="$(echo $x | sed 's/.R1./.R2./g')"
	echo "$x and $R2"
	bowtie2 --quiet --very-sensitive-local -x "$reference" -1 "$readsdir"/"$x" -2 "$readsdir"/"$R2" --al-conc "$name"_aligned_pairs --al "$name"_aligned_single_reads.fastq -p 32 
done 



