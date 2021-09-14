#!/bin/tcsh
#PBS -S /bin/tcsh
#PBS -l mem=32g
#PBS -l vmem=48g
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load tcl/8.6.0 bedtools/2.27.1


####### INPUT #####
set inBed=$1
###################

set outFile = ${inBed:r}.annotated.tsv
if ( ! -f $outFile ) then 
	echo $outFile
	$ANNOTSV/bin/AnnotSV -SVinputFile ${inBed} -SVinputInfo 1 -outputFile ${outFile} -svtBEDcol 4
endif

echo "Successfully annotated!"
