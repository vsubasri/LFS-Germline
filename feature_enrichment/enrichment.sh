#!/bin/bash
#PBS -S /bin/bash
#PBS -N Feature_Enrichment
#PBS -e Feature_Enrichment.err
#PBS -o Feature_Enrichment.out

module load kentutils/302.1 bedtools/2.27.1

win=10
window=/hpf/largeprojects/adam/projects/lfs/lfs_germline/analysis/feature_enrichment/hg19.windows.${win}kb.bed
predictive=/hpf/largeprojects/adam/projects/lfs/lfs_germline/analysis/feature_enrichment/sigCSCE_2612.bed
sampled=/hpf/largeprojects/adam/projects/lfs/lfs_germline/analysis/feature_enrichment/methProbesSampled100k.bed 
outdir=$wd/output

cd $wd

file=$mark/$sample-$mark.imputed.pval.signal

if [ $mark == "RNAseq" ] ; then
	file=$mark/$sample-$mark.imputed.LogRPKM.signal
fi

bigWigToBedGraph ${file}.bigwig ${file}.bed
bedtools map -a $window -b  ${file}.bed -c 4 -o sum >  ${file}.${win}kb.bed
bedtools intersect -a  ${file}.${win}kb.bed -b $predictive > $outdir/$mark/predictive_$sample-$mark.bed
bedtools intersect -a ${file}.${win}kb.bed -b $sampled > $outdir/$mark/sampled_$sample-$mark.bed

