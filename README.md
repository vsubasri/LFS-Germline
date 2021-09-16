# LFS Germline #

Li-Fraumeni syndrome (LFS) is an autosomal dominant cancer-predisposition syndrome associated with pathogenic germline variants in the TP53 tumour suppressor gene. However, 30% of patients that fit the clinical definition of LFS lack a germline pathogenic variant in TP53 and even among variant-TP53 carriers, approximately 20% remain cancer-free. We therefore characterized the germline genetic and/or epigenetic landscape of LFS patients harbouring pathogenic variants of (n=374) or wildtype TP53 (n=22). This directory includes the scripts used to perform these analyses.

## Variant Calling and Filtering ##

    germline_pipeline/

### SNV and indels ###
SNV and indel detection using GATK v.4.0.2.1 according to best practices for germline cohort data. All detected variants are further filtered to remove false positives using a set of filters designed for short read sequences as follows: read position 10-90, strandedness 1%-99%, distance to 3’ > 20, homopolymer <5, map quality difference <30, read length difference <25, and MMQS difference<100. Variants with less than 5 alternate reads detected using bam-readcount are removed. Population filters are applied to remove common variants found in non-cancer “normal” populations: gnomAD, ExAC and 1kGP, at an allele frequency > 0.01. InterVar, an implementation of the ACMG-AMP guidelines, is used to classify variants into 5 categories: pathogenic (P), likely pathogenic (LP), variant of uncertain significance (VUS), likely benign (LB) and benign (B). 

### Structural variants ###
Germline structural variant (SV) discovery is performed using DELLY v0.7.7 and Manta v1.0.3 in order to discover deletions, duplications and inversions. An ensemble approach was taken, retaining only the deletions, duplications and inversions called by both DELLY and Manta with a maximum allowed distance of 1kb between breakpoints. The following size-specific constraints are used for breakpoint consensus between the two tools: breakpoints within 100 bp for SVs ≤10 kbp, breakpoints within 1 kbp for SVs ≥10 kbp and ≤ 50 kbp and breakpoints within 10 kbp for SVs ≥ 50 kbp. A high-quality SV set was obtained by applying additional filtering criteria. A panel of SVs occurring as a result of sequencing artifacts was created to remove technical false positives using 72 germline genomes sequenced on the same HiSeqX platform (150bp paired end sequencing, minimum of 30X depth coverage). SVs present in ≥ 3 samples in the panel of normals are removed. Annotation is performed with AnnotSV and subsequently filtered to remove SVs with an AnnotSV score < 3 and present in publically available non-cancer “normal” databases: 1000 Genomes Project and gnomAD.

### Copy number variation ###
CNV is identified using read depth based CNV detection algorithms: CNVnator and ERDS. The intersection of the CNVnator and ERDS calls are subsequently filtered to remove repetitive and low-complexity regions using a predefined list formed by Trost et al. This is followed by annotation using ANNOVAR.

## Pathway Analysis ##

    pathway_analysis/

The hallmark gene sets from the MSigDB collections are utilized to perform a pathway analysis. Fisher's exact test is implemented for each of the 50 pathways, followed by FDR correction.

## Polygenic Risk Score ##

    prs/

145,583 GWAS rsSNPs from the EMBL-EBI GWAS Catalog42 are filtered for case-control studies that have identified rsSNPs associated with cancer development using the keywords: osteosarcoma, cancer, glioma, leukemia, carcinoma, malignancies, lymphoma, glioblastoma, neuroblastoma, tumor, and sarcoma; this resulted in 4,364 rsSNPs. In addition, rsSNPs evaluating the occurrence of specific LFS cancers in our dataset are also evaluated: osteosarcoma (6), glioma (31), breast cancer (490), non-melanoma skin cancers (53), and leukemia (117). PLINK 1.9 was used for LD clumping and p-value thresholding at 9 thresholds (1x10-2, 1x10-3, 1x10-4, 1x10-5, 1x10-6, 1x10-7, 1x10-10, 1x10-13, 1x10-16). A mixed effects logistic regression was fit for each PRS group using the lme4 package in R44, modelling family as a random effect and polygenic risk score, p53 variant status and sex as fixed effects. Bootstrapped 95% confidence intervals are calculated by resampling our data, taking 100 replicates and refitting the model to the resampled data.


## Epimutations ##

    cisCSCE/

###Preprocessing and addressing confounders in the methylation data
Raw beta values are corrected for dye-bias using ssNoob. The normalized beta values are then corrected for the batch effects between 96-well plates using ComBat. Probabilistic estimation of expression residuals (PEER)–a factor analysis method that infers hidden determinants and their effects on molecular profiles–was used to remove broad variance from known confounders (array type, batch, cancer status, age of sample collection, gender) as well as 100 hidden (latent) factors. Following PCA transformation of each cohort, samples within 3 standard deviations of the mean of PC1 (μPC1 ± 3σ) and PC2 (μPC2 ± 3σ) are retained. 

### Discovery of cancer-associated secondary constitutional epimutation (CSCE) ###

Using the preprocessed methylation data, each of the 452,497 methylation probes are tested for their association with cancer status. Diagnostic methylation probes are first identified in the discovery cohort using the Aziz Test (methylation ~ cancer status). P-values are adjusted using FDR and probes with FDRdiscovery < 0.1 are considered significant. Diagnostic methylation probes are then validated independently in the internal and external validation cohorts. Patients with both WGS and methylation are used to identify cancer-associated secondary constitutional epimutations (CSCE). A cis-CSCE was defined as a diagnostic probe significantly associated with a SNP within 10 kb upstream or downstream (cis-SNP). This association between a diagnostic probe and a cis-SNP was determined using Spearman correlation. Each SNP was coded additively as 1 (homozygous reference), 3 (heterozygous) or 4 (homozygous alternative). For each diagnostic probe the strongest cis-CSCE was retained to represent that probe. Subsequently, FDR correction was performed on all the cis-CSCE.

## Feature Enrichment ##

    feature_enrichment/

Feature enrichment analysis of cis-CSCE using 10 genomic properties, which included proximity to histone marks, open chromatin and expression levels from blood-derived samples. ChromImpute p-value signal tracks (bigwig files) are downloaded from the Roadmap Epigenomics Consortium (https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/) for the following genomic properties: H3K4me1, H3K4me3, H3K27me3, H3K9me3, H3K27ac, H3K36me3, DNase, H2A.Z, H3K79me2 and RNA-sequencing, in 23 blood-derived samples (E062, E034, E045, E033, E044, E043, E039, E041, E042, E040, E037, E048, E038, E047, E029, E031, E035, E051, E050, E036, E032, E046, E030). 100,000 methylation probes are sampled with replacement to calculate the baseline enrichment level. Average signal values are calculated at the genomic loci of the 100,000 probes for all features in each sample. Average signal values are calculated at the genomic loci of the significant cis-CSCE probes for all features in each sample. A paired Wilcoxon Test is used to calculate a p-value between the signal value of the null probes and the cis-CSCE probes for each feature. Similarly, effect size is calculated using Cohen’s distance (d).

## Telomere Analysis ##

    telomere/

Telseq (v0.0.1), a software package designed to calculate mean telomere length from WGS data, was run on the germline genomes. We used a read-length of 150 and a threshold of 12 telomeric repeats per read as parameters when running Telseq. 

## Other analyses ##

    analysis/

Collection of scripts used for analyses including:

- identification of disease segregating variants within families
- organize and filter kics patients
- extracting cancer associated risk SNPs in literature
- differential expression of patient's tumour

## Figures ##

    figures/

Figure 1: The LFS Cohort

Figure 2: TP53 Landscape of LFS

Figure 3: Association of genetic variants with phenotypic properties

Figure 4: Pathway analysis by cancer status and p53 status

Figure 5: Secondary constitutional epimutations in LFS

Figure 6: Polygenic Risk Score in LFS

