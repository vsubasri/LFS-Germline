#!/bin/bash
#PBS -S /bin/bash
#PBS -N PLINK
#PBS -o PLINK.out
#PBS -e PLINK.err
#PBS -l nodes=1:ppn=4
#PBS -l gres=localhd:10

module load plink/1.90b3x
cd /hpf/largeprojects/adam/projects/lfs/lfs_germline/wt/plink
bfile=lfs
outdir=/hpf/largeprojects/adam/projects/lfs/lfs_germline/wt/plink/results

plink \
    --bfile ${bfile} \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump gwas_catalogue_cancerSNPs \
    --clump-snp-field ID \
    --clump-field PVALUE \
    --out lfs.clumped

plink --ped ${bfile}.ped --map ${bfile}.map --make-bed --out ${bfile}

plink \
    --bfile ${bfile} \
    --score gwas_catalogue_cancerSNPs 1 2 3 header \
    --q-score-range pval_range_list SNP.pvalue \
    --extract ${bfile}.valid.snp \
    --out ${outdir}/lfs.out

plink \
    --bfile ${bfile} \
    --score gwas_catalogue_cancerSNPs 1 2 3 header \
    --q-score-range pval_range_list SNP.pvalue \
    --extract bc.valid.snp \
    --out ${outdir}/bc.out

plink \
    --bfile ${bfile} \
    --score gwas_catalogue_cancerSNPs 1 2 3 header \
    --q-score-range pval_range_list SNP.pvalue \
    --extract os.valid.snp \
    --out ${outdir}/os.out

plink \
    --bfile ${bfile} \
    --score gwas_catalogue_cancerSNPs 1 2 3 header \
    --q-score-range pval_range_list SNP.pvalue \
    --extract ALL.valid.snp \
    --out ${outdir}/ALL.out

plink \
    --bfile ${bfile} \
    --score gwas_catalogue_cancerSNPs 1 2 3 header \
    --q-score-range pval_range_list SNP.pvalue \
    --extract leukemia.valid.snp \
    --out ${outdir}/leukemia.out

plink \
    --bfile ${bfile} \
    --score gwas_catalogue_cancerSNPs 1 2 3 header \
    --q-score-range pval_range_list SNP.pvalue \
    --extract nonmelanoma_skin.valid.snp \
    --out ${outdir}/nonmelanoma_skin.out

plink \
    --bfile ${bfile} \
    --score gwas_catalogue_cancerSNPs 1 2 3 header \
    --q-score-range pval_range_list SNP.pvalue \
    --extract glioma.valid.snp \
    --out ${outdir}/glioma.out



