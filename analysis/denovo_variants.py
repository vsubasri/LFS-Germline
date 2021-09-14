import pandas as pd
import os
import subprocess

dir="/Users/vsubasri/research/lfs_wt/structural_variations/final_svs"

def get_unaffs(dir,unaff_ids,filetype):
	unaff_vars = []
	for id in unaff_ids:
		unaff_path = os.path.join(dir,id,id+".annotated."+filetype)
		unaff = pd.read_csv(unaff_path,sep='\t')
		unaff_vars.append(unaff)
	unaff_vars = pd.concat(unaff_vars)
	return(unaff_vars)

def get_cancervars(dir,cancer_id,unaff_ids,filetype):
	os.chdir(os.path.join(dir,cancer_id))
	unaff_vars = get_unaffs(dir,unaff_ids,filetype)
	cancer_vars = pd.read_csv(os.path.join(cancer_id+'.annotated.'+filetype),sep='\t')
	cols = [ c for c in unaff_vars if c not in ['AnnotSV ID']]  + ['AnnotSV ID']
	unaff = unaff_vars[cols] ; unaff.to_csv('unaff_tmp.bed',sep='\t',index=False,header=False)
	cancer = cancer_vars[cols] ; cancer.to_csv('cancer_tmp.bed',sep='\t',index=False,header=False)
	bashCommand = ["bedtools","intersect","-f", "0.9","-F","0.9","-v","-a",'cancer_tmp.bed',"-b","unaff_tmp.bed"]
	outbed = open(cancer_id+".canceronly_"+filetype+"variants.txt", "w") ; outbed.write('\t'.join(cols)+'\n')
	outbed = open(cancer_id+".canceronly_"+filetype+"variants.txt", "a") ; subprocess.call(bashCommand, stdout=outbed)
	os.remove('unaff_tmp.bed') ; os.remove('cancer_tmp.bed')


def get_all_cancervars(dir,iddict):
	for cancer_id,unaff_ids in iddict.items():
		print("Getting cancer only variants for {}...".format(cancer_id))
		get_cancervars(dir,cancer_id,unaff_ids,"pathogenic")
		get_cancervars(dir,cancer_id,unaff_ids,"fpfiltered")


iddict = {	'1777' : ["1778","1792"] ,
			'1774' : ["1778","1792"],
			'735' : ["744"],
			'830_1' : ["833","835"],
			'832' : ["833","835"],
			'3024': ["3025"],
			'4093': ["3025"],
			'3021': ["3025"],
			'3024': ["3025"],
			'1971': ["1972","1973"],
			'1970':["1972","1973"]
			}

get_all_cancervars(dir,iddict)



