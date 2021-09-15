import pandas as pd
import os
import subprocess
import glob 

vars_dir="/hpf/largeprojects/adam/projects/lfs/lfs_germline/wt/gatk/filtered_variants/intervar"
svs_dir="/hpf/largeprojects/adam/projects/lfs/lfs_germline/wt/structural_variations/final_svs"

def get_exonic_CPG(path_outfile,prefix):
        path = pd.read_csv(path_outfile,sep='\t',dtype=str)
        path_exonic = path[path['Func.refGene'].isin(["exonic","splicing"])]
        path_exonic_CPG = path_exonic[path_exonic['CPG'] != "."]
        path_exonic = path_exonic[~path_exonic['ExonicFunc.refGene'].isin(['synonymous SNV','nonframeshift deletion','nonframeshift insertion'])]
        cpg_outpath = os.path.join(os.path.dirname(path_outfile),prefix+".exonic.cpg.txt")
        path_exonic_CPG.to_csv(cpg_outpath,sep='\t',index=False)
        exonic_outpath = os.path.join(os.path.dirname(path_outfile),prefix+".exonic.txt")
        path_exonic.to_csv(exonic_outpath,sep='\t',index=False)

def combine_vars(dir,pattern):
	all_files = []
	for file in glob.glob(os.path.join(dir,"*"+pattern)):
		vars = pd.read_csv(file,sep='\t')
		vars['Sample'] = os.path.basename(file).split('.')[0]
		all_files.append(vars)
	all_vars = pd.concat(all_files)
	all_vars.to_csv(os.path.join(dir,"all"+pattern),sep='\t',index=False)
	

def pathogenic_multisample(intervar_dir):
        m_intervar = os.path.join(intervar_dir,"*.hg19_multianno.txt.intervar")
        intervar_files = glob.glob(m_intervar)
        for file in intervar_files:
                print("Filtering {}".format(file))
                path_outfile = find_pathogenic(file)
                sample=os.path.splitext(path_outfile)[0]
                get_exonic_CPG(path_outfile,sample)
        combine_vars(intervar_dir,".hg19.pathogenic.vus.exonic.cpg.txt")
        combine_vars(intervar_dir,".hg19.pathogenic.vus.exonic.txt")

def find_pathogenic(m_path):
        search_status = ["Likely pathogenic", "Pathogenic", "Uncertain significance"]
        path_outfile = os.path.join(os.path.dirname(m_path), os.path.basename(m_path).split('.')[0] + ".hg19.pathogenic.vus.txt")

        with open(m_path, 'r') as in_f, open(path_outfile, 'w') as out_f:
                firstline=in_f.readline()
                indices=firstline.strip().split("\t")
                out_f.write(firstline)
                for line in in_f:
                        line_indices = line.strip().split("\t")
                        if any(x in line_indices[17] for x in search_status) and (not any(x in line_indices[13] for x in ["Benign","Likely benign"])):
                                out_f.write(line)
        return(path_outfile)

def get_cancer_vars(dir,cancer_id,unaff_ids,filetype):
	unaff_vars = get_unaffs(dir,unaff_ids,filetype)
	cancer_vars = pd.read_csv(os.path.join(cancer_id+'.hg19.pathogenic.vus.txt'),sep='\t')
	unaff_vars.to_csv('unaff_tmp.bed',sep='\t',index=False,header=False)
	cancer_vars.to_csv('cancer_tmp.bed',sep='\t',index=False,header=False)
	bashCommand = ["bedtools","intersect","-v","-a",'cancer_tmp.bed',"-b","unaff_tmp.bed"]
	outbed = open(cancer_id+".canceronly_variants.txt", "w") ; outbed.write('\t'.join(cancer_vars.columns)+'\n')
	outbed = open(cancer_id+".canceronly_variants.txt", "a") ; subprocess.call(bashCommand, stdout=outbed)
	os.remove('unaff_tmp.bed') ; os.remove('cancer_tmp.bed')

def get_all_cancer_vars(dir,iddict):
	os.chdir(dir)
	for cancer_id,unaff_ids in iddict.items():
		print("Getting cancer only variants for {}...".format(cancer_id))
		get_cancer_vars(dir,cancer_id,unaff_ids,".hg19.pathogenic.vus.txt")
		canceronly_file = os.path.join(dir,cancer_id+".canceronly_variants.txt")
		prefix = cancer_id + ".canceronly_variants"
		get_exonic_CPG(canceronly_file,prefix)	
	combine_vars(dir,".canceronly_variants.exonic.cpg.txt")
	combine_vars(dir, ".canceronly_variants.exonic.txt")

def get_unaffs(dir,unaff_ids,filetype):
        unaff_vars = []
        for id in unaff_ids:
                unaff_path = os.path.join(dir,id,id+filetype)
                unaff = pd.read_csv(unaff_path,sep='\t')
                unaff_vars.append(unaff)
        unaff_vars = pd.concat(unaff_vars)
        return(unaff_vars)


def get_cancer_svs(dir,cancer_id,unaff_ids,filetype):
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


def get_all_cancer_svs(dir,iddict):
	for cancer_id,unaff_ids in iddict.items():
		print("Getting cancer only variants for {}...".format(cancer_id))
		get_cancer_svs(dir,cancer_id,unaff_ids,".annotated.pathogenic")
		get_cancer_svs(dir,cancer_id,unaff_ids,".annotated.fpfiltered")


iddict = {
		## Family 3
		'1777' : ["1778","1792"] ,
		'1774' : ["1778","1792"],
		## Family 1
		'735' : ["744"],
		## Family 8
		'830_1' : ["833","835"],
		'832' : ["833","835"],
		## Family 6
		'3024': ["3025"],
		'4093': ["3025"],
		'3021': ["3025"],
		## Family 2
		'1971': ["1972","1973"],
		'1970':["1972","1973"],
		## Family 7
		'832': ["833","835"],
		'830_1': ["833","835"]

pathogenic_multisample(dir)
get_all_cancer_vars(vars_dir,iddict)			}

get_all_cancer_svs(svs_dir,iddict)


