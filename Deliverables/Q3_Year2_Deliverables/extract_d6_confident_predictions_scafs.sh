

#ok updated script to get all the confident virus predictions and extract the transcripts from the fasta files. make sure to edit the path to the fasta files
#DT25_OFRA_1H_August_2022_scaffold_500.fasta
fastdir="/work/alex/florida/RNA/assemblies/trinity/500"

for csv in *_d6_results.csv;
do	awk -F "," '{print $1}' $csv > sample.list
	name=$(echo "$csv" | awk -F "_d6" '{print $1}')
	#generate high confidence prediction file 
	for x in $(cat sample.list)
	do	#echo "$x"
		echo "Top hit: $(grep -w "$x" "$csv" | awk -F "," '{print $3}')"
		top="$(grep -w "$x" "$csv" | awk -F "," '{print $3}')"
		#echo "Checking if above .7"
		if (( $(echo ""$top" > .7" | bc -l) )); 
		then	#echo "Yes"
			grep -w "$x" "$csv" >> "$name"_deep6_confident_predictions.csv
		fi
	done
	#generate high confidence virus transcript csv
	grep -v "euk" "$name"_deep6_confident_predictions.csv | grep -v "pro" > "$name"_deep6_confident_viral_predictions.csv
	#generate high confidence virus transcript list for extraction
	awk -F "," '{print $1}' "$name"_deep6_confident_viral_predictions.csv > tran.list
	seqtk subseq "$fastdir"/"$name"* tran.list > "$name"_deep6_confident_viral_transcripts.fasta
	rm tran.list
done
