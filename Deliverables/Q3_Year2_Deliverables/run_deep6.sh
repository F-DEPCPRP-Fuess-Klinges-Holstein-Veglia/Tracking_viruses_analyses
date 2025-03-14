#!/bin/bash

deep6_path="/home/.progs/Deep6/Master/deep6.py";
deep6_models="/home/.progs/Deep6/Models"
length=500
allvall="/home/.progs/Deep6/deep6_compare_all_vs_all.py"

for x in *.fasta; 
do      echo "Starting lenght filtering with $x; keeping $length and up my g."
        name=$(echo "$x" | awk -F ".fasta" '{print $1}')
        bbduk.sh in="$x" out="$name"_transcripts_"$length"up.fasta minlen="$length"     
        echo "--------------------------------------------------------"
	echo "Starting deep6 with "$name"_transcripts_"$length"up.fasta"
        python $deep6_path -i "$name"_transcripts_"$length"up.fasta -l $length -m $deep6_models -o "$name"_d6
	sed 's/\t/,/g' ./"$name"_d6/"$name"_transcripts_"$length"up.fasta_predict_500bp_deep6.txt >> "$name"_transcripts_"$length"up.fasta_predict_500bp_deep6.csv
	echo "--------------------------------------------------------"
	echo "Starting with comparison jawn"
	python $allvall ""$name"_transcripts_"$length"up.fasta_predict_500bp_deep6.csv" "$name"_checkcsv.csv "$name"_d6_results.csv
done
