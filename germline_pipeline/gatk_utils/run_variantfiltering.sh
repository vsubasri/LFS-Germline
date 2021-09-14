#!/bin/bash
#PBS -S /bin/bash
#PBS -N RunVariantFiltering
#PBS -o log/RunVariantFiltering.out
#PBS -e log/RunVariantFiltering.err

bcftools_filter=/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/bcftools_filter.sh
split_vcf=/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/split_vcf.sh
combine_variants=/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/combine_variants.sh
combine_snvs_indels=/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/combine_snvs_indels.sh
fpfilter=/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/submit_fpfilter.sh
intervar=/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/run_intervar.sh
postfiltering=/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/run_postfiltering.sh

cd $GATKDIR

first=$(qsub -v gatkdir=$GATKDIR,multisample=$MULTISAMPLE ${split_vcf})
echo $first


if [[ -d $BAMDIR ]]; then
        bam_files=( $(find $BAMDIR -type f -name "*.bam") )
else
        mapfile -t bam_files < $BAMDIR
fi

output_dir=${GATKDIR}/filtered_variants/
if [ ! -d ${output_dir} ]; then
        mkdir ${output_dir}
fi
indel_dir=${output_dir}sample_indels/
if [ ! -d ${indel_dir} ]; then
        mkdir ${indel_dir}
fi
snv_dir=${output_dir}sample_snvs/
if [ ! -d ${snv_dir} ]; then
        mkdir ${snv_dir}
fi
filtered_dir=${output_dir}fpfilter_snvs/
if [ ! -d ${filtered_dir} ]; then
        mkdir ${filtered_dir}
fi
final_dir=${output_dir}final_variants/
if [ ! -d ${final_dir} ]; then
        mkdir ${final_dir}
fi
intervar_dir=${output_dir}intervar/
if [ ! -d ${intervar_dir} ]; then
        mkdir ${intervar_dir}
fi
if [ ! -d ${output_dir}log/ ]; then
	mkdir ${output_dir}log/
fi

cd $output_dir

COUNTER=0
for bam_file in ${bam_files[*]}
do

	COUNTER=$[$COUNTER +1]
	##sample name
	sample=$(echo $(basename $bam_file) | cut -d '.' -f 1)
	##input file
        sample_path="${GATKDIR}/${sample}/${sample}.vcf"
        ##bcftool filter output file
	bcftools_snvs="${snv_dir}${sample}_SNV.vcf"
	bcftools_indels="${indel_dir}${sample}_indel.vcf"
	##fpfilter output file 
	fpfilter_snvs="${filtered_dir}${sample}_fpfilter.vcf"
	final_variants="${final_dir}${sample}.vcf"

	if [ -f $final_variants ] ; then
		continue
	fi

	echo "Processing $sample ..."
	second=$(qsub -W depend=afterok:$first -v indel_dir=$indel_dir,snv_dir=$snv_dir,sample_path=$sample_path,sample=$sample ${bcftools_filter})
        echo $second
	third=$(qsub -W depend=afterok:$second -v sample=$sample,bam=$bam_file,vcf=$bcftools_snvs,output=$fpfilter_snvs ${fpfilter})
        echo $third

	fourth=$(qsub -W depend=afterok:$third -v snvs=$fpfilter_snvs,indels=$bcftools_indels,output=$final_variants ${combine_snvs_indels})
	echo $fourth

        logoutfile=${intervar_dir}/log/${sample}.out
        logerrfile=${intervar_dir}/log/${sample}.err
	qsub -W depend=afterok:$fourth -o $logoutfile -e $logerrfile -v input=$final_variants,output=${intervar_dir}${sample}.anno ${intervar}

	if [ ${#bam_files[@]} = $COUNTER ]; then
		dependencies=$dependencies$fourth
	else 
		dependencies=$dependencies$fourth:
	fi
done

fifth=$(qsub -v intervar_dir=${intervar_dir} -W depend=afterok:${dependencies} ${postfiltering}) 
echo $fifth
sixth=$(qsub -v var_dir=${output_dir} -W depend=afterok:${dependencies} ${combine_variants})
echo $sixth

