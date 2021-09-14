#!/bin/bash
#PBS -S /bin/bash
#PBS -N combinePrefiltered
#PBS -e log/combinePrefiltered.err
#PBS -o log/combinePrefiltered.out
#PBS -l vmem=8g
#PBS -l walltime=23:59:00

module load tcl/8.6.0 bcftools/1.6 bedtools/2.27.1 samtools

outdir=$sv_dir/final_svs

sample=$(echo $(basename $bam_file) | cut -d '.' -f 1)
manta_vcf=${sv_dir}/manta/${sample}/results/variants/diploidSV.vcf

delly_vcf=${sv_dir}/delly/data/${sample}/delly.${sample}.vcf

if [ ! -f $delly_vcf ]; then
	sample_delly=""
	for delly_file in ${sv_dir}/delly/data/${sample}/normal*.tab.vcf
	do
		if [ $(basename $delly_file) == normal_${sample}.tra.tab.vcf ] ; then
			continue
		fi
		sample_delly+="${delly_file}.gz "
		bgzip -c $delly_file > ${delly_file}.gz
		tabix -p vcf ${delly_file}.gz
	done
	bcftools concat -a $sample_delly -o $delly_vcf
fi

if [ ! -d $outdir/$sample ] ; then
	mkdir $outdir/$sample
fi
	
if [ ! -f $outdir/$sample/${sample}.vcf ] ; then
	cp ${manta_vcf}.gz tmp.gz ; gunzip tmp.gz ; mv tmp $manta_vcf
	printf '%s\n' $manta_vcf $delly_vcf > $outdir/$sample/input.list
	/hpf/largeprojects/adam/projects/lfs/code/SURVIVOR/Debug/SURVIVOR merge \
		$outdir/$sample/input.list 1000 2 1 1 1 0 $outdir/$sample/${sample}.vcf	
fi

outFile=$outdir/$sample/${sample}.annotated.tsv

if [ ! -f $outFile ] ; then
	echo [ $sample ]
	$ANNOTSV/bin/AnnotSV -SVinputFile $outdir/$sample/${sample}.vcf -SVinputInfo 1 -outputFile ${outFile}
fi

