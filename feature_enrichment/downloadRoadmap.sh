#!/bin/bash
#PBS -S /bin/bash
#PBS -l vmem=4g
#PBS -l mem=4g
#PBS -N Feature_Enrichment
#PBS -e Feature_Enrichment.err
#PBS -o Feature_Enrichment.out

module load kentutils/302.1 bedtools/2.27.1

marks=(RNAseq H3K4me1 H3K4me3 H3K27me3 H3K9me3 H3K27ac H3K36me3 DNase H2A.Z H3K79me2)

function validate_url(){
  if [[ `wget -S --spider $1  2>&1 | grep 'HTTP/1.1 200 OK'` ]]; then
    return 0
  else
    return 1
  fi
}

for mark in ${marks[@]}
do
	if [ ! -d $mark ]; then
		mkdir $mark
	fi
	
	for i in E062 E034 E045 E033 E044 E043 E039 E041 E042 E040 E037 E048 E038 E047 E029 E031 E035 E051 E050 E036 E032 E046 E030;
	do 
		echo "$i $mark"
		file=$mark/$i-$mark.imputed.pval.signal.bigwig

		if [ $mark == "RNAseq" ] ; then
			file=$mark/$i-$mark.imputed.LogRPKM.signal.bigwig
		fi

		if validate_url https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/$file && [ ! -f $file ]; then
	        	echo https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/$file
			wget https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/$file -P $mark
		fi
	done
done



