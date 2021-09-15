# LFS Germline #

Li-Fraumeni syndrome (LFS) is an autosomal dominant cancer-predisposition syndrome associated with pathogenic germline variants in the TP53 tumour suppressor gene. However, 30% of patients that fit the clinical definition of LFS lack a germline pathogenic variant in TP53 and even among variant-TP53 carriers, approximately 20% remain cancer-free. We therefore characterized the germline genetic and/or epigenetic landscape of LFS patients harbouring pathogenic variants of (n=374) or wildtype TP53 (n=22). This directory includes the scripts used to perform these analyses.

## Variant Calling and Filtering ##

germline_pipeline/

### SNV and indels ###
SNV and indel detection was performed using GATK v.4.0.2.1 according to best practices for germline cohort data.27 All detected variants were further filtered to remove false positives using a set of filters designed for short read sequences as follows: read position 10-90, strandedness 1%-99%, distance to 3’ > 20, homopolymer <5, map quality difference <30, read length difference <25, and MMQS difference<100.28 Variants with less than 5 alternate reads detected using bam-readcount were removed. Population filters were applied to remove common variants found in non-cancer “normal” populations: gnomAD, ExAC and 1kGP, at an allele frequency > 0.01.29 InterVar, an implementation of the ACMG-AMP guidelines,30 was used to classify variants into 5 categories: pathogenic (P), likely pathogenic (LP), variant of uncertain significance (VUS), likely benign (LB) and benign (B). 

### Structural variants ###
Germline structural variant (SV) discovery was performed using DELLY v0.7.7 and Manta v1.0.3 in order to discover deletions, duplications and inversions. An ensemble approach was taken, retaining only the deletions, duplications and inversions called by both DELLY and Manta with a maximum allowed distance of 1kb between breakpoints. Manta relies on read-pair, split-read and local-assembly support to call variants, while Delly uses read-pair and split-read support. The following size-specific constraints were used for breakpoint consensus between the two tools: breakpoints within 100 bp for SVs ≤10 kbp, breakpoints within 1 kbp for SVs ≥10 kbp and ≤ 50 kbp and breakpoints within 10 kbp for SVs ≥ 50 kbp. A high-quality SV set was obtained by applying additional filtering criteria. A panel of SVs occurring as a result of sequencing artifacts was created to remove technical false positives using 72 germline genomes sequenced on the same HiSeqX platform (150bp paired end sequencing, minimum of 30X depth coverage). SVs present in ≥ 3 samples in the panel of normals were removed. Annotation was performed with AnnotSV and subsequently filtered to remove SVs with an AnnotSV score < 3 and present in publically available non-cancer “normal” databases: 1000 Genomes Project and gnomAD.

### Copy number variation ###
CNV was identified using read depth based CNV detection algorithms: CNVnator and ERDS.35,36 The intersection of the CNVnator and ERDS calls were subsequently filtered to remove repetitive and low-complexity regions using a predefined list formed by Trost et al.37 This was followed by annotation using ANNOVAR.

## Pathway Analysis ##

pathway_analysis/

The hallmark gene sets from the MSigDB collections were utilized to perform a pathway analysis. Fisher's exact test was implemented for each of the 50 pathways, followed by FDR correction.

## Polygenic Risk Score ##

prs/

145,583 GWAS rsSNPs from the EMBL-EBI GWAS Catalog42 were filtered for case-control studies that have identified rsSNPs associated with cancer development using the keywords: osteosarcoma, cancer, glioma, leukemia, carcinoma, malignancies, lymphoma, glioblastoma, neuroblastoma, tumor, and sarcoma; this resulted in 4,364 rsSNPs. The LFS cohort primarily consists of individuals of European ancestry, as a result the rsSNPs from studies focusing primarily on non-European individuals were removed. The remaining 3189 rsSNPs were overlapped with the SNPs in the LFS cohort with high confidence genotype calls in all samples, leaving 1,568 rsSNPs to evaluate cancer risk, agnostic of cancer type. In addition, rsSNPs evaluating the occurrence of specific LFS cancers in our dataset were also evaluated: osteosarcoma (6), glioma (31), breast cancer (490), non-melanoma skin cancers (53), and leukemia (117). Due to the rarity of the other cancer types (e.g. adrenocortical carcinoma, choroid plexus carcinoma, rhabdomyosarcoma), there currently exists no studies that have identified associated GWAS SNPs in these patient populations. PLINK 1.9 was used for LD clumping and p-value thresholding at 9 thresholds (1x10-2, 1x10-3, 1x10-4, 1x10-5, 1x10-6, 1x10-7, 1x10-10, 1x10-13, 1x10-16). A mixed effects logistic regression was fit for each PRS group using the lme4 package in R44, modelling family as a random effect and polygenic risk score, p53 variant status and sex as fixed effects. We calculate bootstrapped 95% confidence intervals by resampling our data, taking 100 replicates and refitting the model to the resampled data.


## Epimutations ##

cisCSCE/

###Preprocessing and addressing confounders in the methylation data
Preprocessing and bias correction were performed on the discovery and each validation cohort separately. Raw beta values were corrected for dye-bias using ssNoob. The normalized beta values were then corrected for the batch effects between 96-well plates using ComBat. Probabilistic estimation of expression residuals (PEER)–a factor analysis method that infers hidden determinants and their effects on molecular profiles–was used to remove broad variance from known confounders (array type, batch, cancer status, age of sample collection, gender) as well as 100 hidden (latent) factors. Following PCA transformation of each cohort, samples within 3 standard deviations of the mean of PC1 (μPC1 ± 3σ) and PC2 (μPC2 ± 3σ) were retained. 

### Discovery of cancer-associated secondary constitutional epimutation (CSCE) ###

Using the preprocessed methylation data, each of the 452,497 methylation probes were tested for their association with cancer status as follows:
i. Diagnostic methylation probes were identified in the discovery cohort using the Aziz Test (methylation ~ cancer status). P-values were adjusted using false discovery rate (FDR) and probes with FDRdiscovery < 0.1 were considered significant.
ii. Diagnostic methylation probes were then identified independently in the internal and external validation cohorts using the Aziz Test (methylation ~ cancer status). P-values were adjusted using false discovery rate (FDR) and probes with FDRvalidation < 0.1 in both validation cohorts were considered significant. 
iii. The diagnostic methylation probes found in the discovery cohort with FDRdiscovery < 0.1 that were validated in both the internal and external validation cohorts with FDRvalidation < 0.1 resulted in 931 significant probes. 

The patients with both WGS and methylation (n=47) were used to identify cancer-associated secondary constitutional epimutations (CSCE). A cis-CSCE was defined as a diagnostic probe significantly associated with a SNP within 10 kb upstream or downstream (cis-SNP). This association between a diagnostic probe and a cis-SNP was determined using Spearman correlation. Each SNP was coded additively as 1 (homozygous reference), 3 (heterozygous) or 4 (homozygous alternative). For each diagnostic probe the strongest cis-CSCE was retained to represent that probe. Subsequently, FDR correction was performed on all the cis-CSCE.

## Feature Enrichment ##

feature_enrichment/

Feature enrichment analysis of cis-CSCE using 10 genomic properties, which included proximity to histone marks, open chromatin and expression levels from blood-derived samples. ChromImpute p-value signal tracks (bigwig files) were downloaded from the Roadmap Epigenomics Consortium (https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidatedImputed/) for the following genomic properties: H3K4me1, H3K4me3, H3K27me3, H3K9me3, H3K27ac, H3K36me3, DNase, H2A.Z, H3K79me2 and RNA-sequencing, in 23 blood-derived samples (E062, E034, E045, E033, E044, E043, E039, E041, E042, E040, E037, E048, E038, E047, E029, E031, E035, E051, E050, E036, E032, E046, E030). 100,000 methylation probes were sampled with replacement to calculate the baseline enrichment level. Average signal values were calculated at the genomic loci of the 100,000 probes for all features in each sample. Average signal values were calculated at the genomic loci of the 259 significant cis-CSCE probes for all features in each sample. A paired Wilcoxon Test was used to calculate a p-value between the signal value of the null probes and the cis-CSCE probes for each feature. Similarly, effect size was calculated using Cohen’s distance (d).

## Telomere Analysis ##

telomere/

Telseq (v0.0.1), a software package designed to calculate mean telomere length from WGS data, was run on the germline genomes. TelSeq estimates mean telomere length by counting the number of reads with telomeric repeats (TTAGGG), divided by GC-adjusted coverage multiplied by the mean chromosome size. We used a read-length of 150 and a threshold of 12 telomeric repeats per read as parameters when running Telseq. 

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

