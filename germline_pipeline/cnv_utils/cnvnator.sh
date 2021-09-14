#!/bin/bash
#PBS -S /bin/bash
#PBS -N CNVator
#PBS -e log.err
#PBS -o log.out
#PBS -V
#PBS -l mem=16g
#PBS -l vmem=16g
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load cnvnator/0.3.3 ROOT/6.06.00

chrom=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
sample=$(echo $(basename $bam_file) | sed 's/.realigned-recalibrated.bam/\n/g')

echo "Computing cnvs for sample ${sample} using binsize of ${BINSIZE}"

echo "Readmapping sample ${sample}..."
#preprocessing
cnvnator -root "${sampledir}/${sample}.root.out" -unique -tree ${bam_file}
#generate histogram
echo "Generating histogram..."
cnvnator -root "${sampledir}/${sample}.root.out" -his ${BINSIZE} 
#calculate statistics
echo "Calculating statistics..."
cnvnator -root "${sampledir}/${sample}.root.out" -stat ${BINSIZE} 
#rd signal paritioning
echo "RD Signal Paritioning..."
cnvnator -root "${sampledir}/${sample}.root.out" -partition ${BINSIZE} 
#cnvnator
echo "Running CNVnator..."
cnvnator -root "${sampledir}/${sample}.root.out" -genome GRCh37 -call ${BINSIZE} >> "${sampledir}/${sample}.CNVnator.results"

cnvnator2VCF.pl "${sampledir}/${sample}.CNVnator.results" > "${sampledir}/${sample}.CNVnator.vcf"

