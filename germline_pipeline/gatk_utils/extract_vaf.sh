#!/bin/bash
#PBS -S /bin/bash
#PBS -N ExtractVAF
#PBS -e extract_vaf.err
#PBS -o extract_vaf.out
#PBS -l mem=48g
#PBS -l vmem=48g
#PBS -l walltime=23:59:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

#Ref: https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php#--selectTypeToInclude

module load bedtools/2.25.0 java/1.8.0_91 picard-tools/2.5.0
ref_fasta=/hpf/largeprojects/adam/projects/lfs/resources/human_g1k_v37_decoy.fasta
gatk=/hpf/largeprojects/adam/projects/lfs/resources/gatk-4.0.1.0/gatk-package-4.0.1.0-local.jar
TMPDIR=/localhd/$PBS_JOBID

############################################################################################################################################
#### INPUT: input vcf
#### BED: bed file with coordinates of variants to extract
#### OUTPUT: output vcf
############################################################################################################################################

bedtools sort -i $bed > ${bed}.sorted
intersectBed -a $input -b ${bed}.sorted -header > $output

tabfile=$(echo "$output" | cut -f 1 -d '.').tab

java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $gatk VariantsToTable \
     -R ${ref_fasta} \
     -V ${output} \
     -F CHROM -F POS -F REF -F ALT -F QUAL -F FILTER -F ID -GF GT -GF AD -GF DP \
     -O ${tabfile}

