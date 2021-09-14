#!/bin/bash

#PBS -S /bin/bash
#PBS -N AnnotateCNVs
#PBS -e log/annotate.err
#PBS -o log/annotate.out
#PBS -l mem=32g
#PBS -l vmem=48g
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load annovar/2017.07.16 python/3.5.2

humandb=/hpf/largeprojects/adam/projects/lfs/resources/humandb/

table_annovar.pl --protocol refGene,dgvMerged,gwasCatalog,phastConsElements46way,tfbsConsSites,cytoBand,wgRna,targetScanS,genomicSuperDups,snp138 -operation g,r,r,r,r,r,r,r,r,r --buildver hg19 --outfile ${annovar_output} --remove --nastring . ${annovar_input} ${humandb}

#table_annovar.pl --minqueryfrac 0.5 --protocol refGene,dgvMerged,gwasCatalog,phastConsElements46way,tfbsConsSites,cytoBand,wgRna,targetScanS,genomicSuperDups,snp138 -operation g,r,r,r,r,r,r,r,r,r --buildver hg19 --outfile ${annovar_output} --remove --nastring . ${annovar_input} ${humandb} 

python /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/cnv_utils/annotateCNV.py -a ${annovar_output}
