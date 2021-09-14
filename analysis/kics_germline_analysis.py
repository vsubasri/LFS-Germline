import pandas as pd
import numpy as np
import os

kics_clin = pd.read_csv('~/research/KiCS/kics_lfsspectrumtumours.txt',sep='\t',skipinitialspace=True)
kics_id_map = pd.read_csv('~/research/KiCS/kics_ids.txt',sep='\t', skipinitialspace=True)
kics_id_map = dict(zip(kics_id_map.dna_id,kics_id_map.kics_id))

kics_file='~/research/KiCS/kics.hg19.pathogenic.txt'
thouG_file='~/research/KiCS/1000G.hg19.pathogenic.txt'



def get_pcpgs(file):
	prefix=os.path.basename(file).split('.')[0]
	psnvs = pd.read_csv(file,sep='\t',skipinitialspace=True, engine='python')
	psnvs['kics_id'] = psnvs['sample'].map(kics_id_map)
	psnvs = psnvs[psnvs['CPG']!="."]
	psnvs['gnomAD_genome_ALL'] = np.where(psnvs['gnomAD_genome_ALL'] == ".",None,psnvs['gnomAD_genome_ALL']).astype(float)
	psnvs['Freq_ExAC_ALL'] = np.where(psnvs['Freq_ExAC_ALL'] == ".",None,psnvs['Freq_ExAC_ALL']).astype(float)
	psnvs['Freq_1000g2015aug_all'] = np.where(psnvs['Freq_1000g2015aug_all'] == ".",None,psnvs['Freq_1000g2015aug_all']).astype(float)
	psnvs = psnvs[~((psnvs['gnomAD_genome_ALL'] > 0.001) | (psnvs['Freq_1000g2015aug_all'] > 0.001) | (psnvs['Freq_ExAC_ALL'] > 0.001))]
	psnvs.to_csv(prefix+'.hg19.pathogenic.cpg.txt',sep='\t',index=None)
	return(psnvs)

def get_kics_pvars(kics_file):
	## reformat snv file
	psnvs = get_pcpgs(kics_file)
	psnvs['kics_id'] = psnvs['sample'].map(kics_id_map)
	psnvs_lfsspectrum = psnvs[psnvs['kics_id'].isin(kics_clin['KiCS.ID'])]
	## remove individuals with cancer predisposition
	tp53_kics_ids = psnvs.kics_id[psnvs['Ref.Gene']=="TP53"].tolist() + ["220", "51", "167"]
	paraganglioma_ids = ["98", "222"]
	lynch_ids = ["63",  "83","120",  "141", "156", "171", "232", "219"]
	hereditary_syndrome_ids = tp53_kics_ids + paraganglioma_ids + lynch_ids
	psnvs = psnvs[~psnvs["kics_id"].isin(hereditary_syndrome_ids)]
	p_snvs_tocomb = psnvs[["Ref.Gene","kics_id"]]
	p_snvs_tocomb.columns = ["Gene name","kics_id"]
	## get pathogenic SVs
	svs = pd.read_csv('~/research/KiCS/kics.pathogenic.exonic.txt',sep='\t')
	svs['kics_id'] = svs['sample_id'].map(kics_id_map)
	svs_cpg = svs[svs.CPG == "CPG"]
	## remove INV that does not occur in exon
	svs_cpg = svs_cpg[~svs_cpg['AnnotSV ID'].isin(["9_97947446_98153304_INV","3_169474492_169510369_DUP"])]
	svs_cpg_tocomb = svs_cpg[["Gene name","kics_id"]]
	## combine snvs and svs
	all_p_cpg = pd.concat([p_snvs_tocomb,svs_cpg_tocomb])
	n_cpg = all_p_cpg["kics_id"].nunique() 

	## number of KiCS patients with class 3 variant
	freq_cpg = n_cpg/(208-len(hereditary_syndrome_ids))
	print('KiCS Class 3 Freq: {:f}'.format(freq_cpg))

	lfsspectrum_p_cpg = all_p_cpg[all_p_cpg.kics_id.isin(kics_clin['KiCS.ID'])]
	lfsspectrum_n_cpg = lfsspectrum_p_cpg["kics_id"].nunique() 
	## number of lfs spectrum KiCS patients with class 3 variant
	lfsspectrum_freq_cpg = lfsspectrum_n_cpg/(208-len(hereditary_syndrome_ids))
	print('LFS Spectrum KiCS Class 3 Freq: {:f}'.format(lfsspectrum_freq_cpg))
	

def get_1000G_pvars(thouG_file):
	psnvs = get_pcpgs(thouG_file)
	n_cpg = psnvs["sample"].nunique()
	freq_cpg = n_cpg/2504
	print('1000G Class 3 Freq: {:f}'.format(freq_cpg))

get_kics_pvars(kics_file)
get_1000G_pvars(thouG_file)


