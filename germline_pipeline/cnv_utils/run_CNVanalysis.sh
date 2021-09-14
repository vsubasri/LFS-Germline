#!/bin/bash
#PBS -S /bin/bash
#PBS -N CNVanalysis
#PBS -l mem=8g
#PBS -l vmem=8g
#PBS -e log/CNVAnalysis.err
#PBS -o log/CNVAnalysis.out
#PBS -l nodes=1:ppn=4

module load bedtools/2.27.1 python/3.5.2
cd $cnvnator_erds

##INPUT: cnvnator_erds
RLCR=/hpf/largeprojects/adam/projects/lfs/resources/RLCRs.bed
blacklisted=/hpf/largeprojects/adam/projects/lfs/resources/SKBlacklistedRegions_AffyCytoScanHD_hg19.txt
annotate_script=/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/cnv_utils/annotate_cnvr.sh
annotated_dir=${cnvnator_erds}/annotated/
mkdir $annotated_dir 

python /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/cnv_utils/processCNV.py -c ${cnvnator_erds}

bedtools sort -i ${cnvnator_erds}/cnvnator.dels.tmp >  ${cnvnator_erds}/cnvnator.dels.sorted.tmp
bedtools sort -i  ${cnvnator_erds}/cnvnator.dups.tmp >  ${cnvnator_erds}/cnvnator.dups.sorted.tmp
bedtools sort -i ${cnvnator_erds}/erds.dels.tmp >  ${cnvnator_erds}/erds.dels.sorted.tmp
bedtools sort -i  ${cnvnator_erds}/erds.dups.tmp >  ${cnvnator_erds}/erds.dups.sorted.tmp

bedtools intersect -a ${cnvnator_erds}/erds.dups.sorted.tmp -b ${cnvnator_erds}/cnvnator.dups.sorted.tmp -loj > ${cnvnator_erds}/dups.combined.tmp
bedtools intersect -a ${cnvnator_erds}/cnvnator.dels.sorted.tmp -b ${cnvnator_erds}/erds.dels.sorted.tmp > ${cnvnator_erds}/dels.combined.tmp

bedtools subtract -a  ${cnvnator_erds}/dups.combined.tmp -b ${RLCR} >  ${cnvnator_erds}/dups_noRLCRs.tmp
bedtools subtract -a  ${cnvnator_erds}/dels.combined.tmp -b ${RLCR} >  ${cnvnator_erds}/dels_noRLCRs.tmp

bedtools subtract -a ${cnvnator_erds}/dups_noRLCRs.tmp -b ${blacklisted} > ${cnvnator_erds}/dups_noRLCRs_noblacklisted.tmp
bedtools subtract -a ${cnvnator_erds}/dels_noRLCRs.tmp -b ${blacklisted} > ${cnvnator_erds}/dels_noRLCRs_noblacklisted.tmp


cut -f1,2,3,4 ${cnvnator_erds}/dups_noRLCRs_noblacklisted.tmp > ${cnvnator_erds}/final_dups.tmp
cut -f1,2,3,4 ${cnvnator_erds}/dels_noRLCRs_noblacklisted.tmp > ${cnvnator_erds}/final_dels.tmp

awk '!seen[$0]++' ${cnvnator_erds}/final_dups.tmp > ${cnvnator_erds}/final_dups
awk '!seen[$0]++' ${cnvnator_erds}/final_dels.tmp > ${cnvnator_erds}/final_dels

sed 's/chr//;s/$/&\t0\t0/' ${cnvnator_erds}/final_dels | cut -f1,2,3,5,6 > ${cnvnator_erds}/dels.annovar
sed 's/chr//;s/$/&\t0\t0/' ${cnvnator_erds}/final_dups | cut -f1,2,3,5,6 > ${cnvnator_erds}/dups.annovar

samp_deps=$(cat ${cnvnator_erds}/deps.tmp)
rm ${cnvnator_erds}/*.tmp

files=( $(find ${cnvnator_erds} -maxdepth 1 -type f -name "*.annovar") )
for file in ${files[@]}
do
        annovar_output=$annotated_dir$(basename $file .annovar).anno
	echo $annovar_output
        if [ ! -f $annovar_output ]; then
                dep=$(qsub -v annovar_input=$file,annovar_output=$annovar_output ${annotate_script})
		samp_deps=${samp_deps}:${dep}
        fi
done

qsub -W depend=afterok:$samp_deps -v cnvs=${cnvnator_erds} /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/cnv_utils/format_final.sh

