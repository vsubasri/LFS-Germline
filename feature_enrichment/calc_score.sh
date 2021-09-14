#!/bin/bash

marks=(RNAseq H3K4me1 H3K4me3 H3K27me3 H3K9me3 H3K27ac H3K36me3 DNase H2A.Z H3K79me2)

outdir=/hpf/largeprojects/adam/projects/lfs/lfs_germline/analysis/feature_enrichment/output_mean
outfile=final_scores.txt

cd $outdir
printf "%s\t%s\t%s\n" "GROUP" "SAMPLE" "SCORE" > $outfile
for mark in ${marks[@]}
do
	echo "Calculating for $mark..."
	for file in $outdir/$mark/*bed
	do
		group=$(echo $(basename $file .bed) | cut -d '_' -f 1)
		sample=$(echo $(basename $file .bed) | cut -d '_' -f 2)
		score_sum=$(awk '{sum+=$4} END{printf "%10.0f" ,sum}' $file)
		n=$(wc -l $file | cut -f1 -d' ').0
		if [ $n == "0.0" ] ; then
			continue
		fi
		score=$(bc -l <<< "scale=2;$score_sum/$n")
		printf "%s\t%s\t%s\n" "$group" "$sample" "$score" >> $outfile
	done
done

