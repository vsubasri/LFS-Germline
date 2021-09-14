import pandas as pd
import glob
import os
import itertools
import argparse 
import subprocess
import numpy as np
import scipy.stats

def format_cnvnator(file):
        cnv = pd.read_csv(file, sep ='\t', names=["CNV_type", "coordinates", "CNV_size", "normalized_RD", "e-val1", "e-val2", "e-val3", "e-val4", "q0"])
        cnv[["chrom","start"]]= cnv["coordinates"].str.split(":", expand=True)
        cnv[["start", "end"]] = cnv["start"].str.split("-", expand=True)
        cnv = cnv[["chrom", "start", "end", "CNV_type"]]
        return(cnv)

def preprocess_cnvnator(path):
        files = glob.glob(os.path.join(path, "*","*.CNVnator.results"))
        for file in files:
                cnvnator = format_cnvnator(file)
                cnvnator.to_csv(file + ".bed", sep='\t', index=None, header=None)
                dels = cnvnator[cnvnator.CNV_type == "deletion"]
                dels.to_csv(file + ".dels.bed", sep='\t', index=None, header=None)
                dups = cnvnator[cnvnator.CNV_type == "duplication"]
                dups.to_csv(file + ".dups.bed", sep='\t', index=None, header=None)

def combine(path, pattern, outdir, filename):
	files = glob.glob(os.path.join(path,"*", pattern))
	cnvs = []
	for file in files:
		print(file)
		if "erds" in path:
			cnv = pd.read_csv(file, sep='\t', header=None).iloc[:, 0:5]
			cnv.columns = ["chrom", "start", "end","read_depth","type"]
			cnv.drop('read_depth', axis=1, inplace=True)
		else:
			cnv = pd.read_csv(file, sep='\t',header=None).iloc[:, 0:4]
			cnv.columns = ["chrom", "start", "end", "type"]
		cnv["sample"] = os.path.basename(file).split('.')[0]
		cnvs.append(cnv)
	cnvs = pd.concat(cnvs)
	cnvs.to_csv(os.path.join(outdir,filename), sep='\t', header=None, index=None)


def format_chr(file):
	cnv = pd.read_csv(file, sep='\t', header=None)
	if not cnv[0].str.contains("chr").any():
		cnv[0] = "chr"+ cnv[0]
	cnv.to_csv(file, sep='\t', index=None, header=None)

def process_sample(outdir, cnvnatordels, cnvnatordups, erdsdels, erdsdups):
	RLCR="/hpf/largeprojects/adam/projects/lfs/resources/RLCRs.bed"
	sample = os.path.basename(outdir)
	for file in [cnvnatordels, cnvnatordups, erdsdels, erdsdups]:
		format_chr(file)
	with open(os.path.join(outdir,"cnvnator.dels.sorted.tmp"), "w") as f:
		subprocess.call(["bedtools","sort", "-i", cnvnatordels], stdout = f)
	with open(os.path.join(outdir,"cnvnator.dups.sorted.tmp"), "w") as f:
		subprocess.call(["bedtools","sort", "-i", cnvnatordups], stdout = f)
	with open(os.path.join(outdir,"erds.dels.sorted.tmp"), "w") as f:
               subprocess.call(["bedtools","sort", "-i", erdsdels], stdout = f)
	with open(os.path.join(outdir,"erds.dups.sorted.tmp"), "w") as f:
                subprocess.call(["bedtools","sort", "-i", erdsdups], stdout = f)
	with open(os.path.join(outdir,"dels.combined.tmp"), 'w') as f:
		subprocess.call(["bedtools", "intersect", "-a", os.path.join(outdir,"erds.dels.sorted.tmp"), "-b", os.path.join(outdir,"cnvnator.dels.sorted.tmp")] , stdout=f)
	with open(os.path.join(outdir, "dups.combined.tmp"), 'w') as f:
		subprocess.call(["bedtools", "intersect", "-a", os.path.join(outdir,"erds.dups.sorted.tmp"), "-b", os.path.join(outdir,"cnvnator.dups.sorted.tmp"), "-loj"], stdout=f)
	with open(os.path.join(outdir,"dels_noRLCRs.tmp"), 'w') as f:
		subprocess.call(["bedtools", "subtract", "-a", os.path.join(outdir,"dels.combined.tmp"), "-b", RLCR] , stdout=f)
	with open(os.path.join(outdir, "dups_noRLCRs.tmp"), 'w') as f:
		subprocess.call(["bedtools", "subtract", "-a", os.path.join(outdir,"dups.combined.tmp"), "-b", RLCR], stdout=f)
	dels = pd.read_csv(os.path.join(outdir,"dels_noRLCRs.tmp"), sep='\t', header=None).iloc[:, 0:3]
	dels.to_csv(os.path.join(outdir, "final_dels"), sep='\t', index=None, header=None)
	dels[3] = 0 ; dels[4] = 0
	dels.to_csv(os.path.join(outdir, "dels.annovar"), sep='\t', index=None, header=None)
	dups = pd.read_csv(os.path.join(outdir,"dups_noRLCRs.tmp"), sep='\t', header=None).iloc[:, 0:3]
	dups.to_csv(os.path.join(outdir, "final_dups"), sep='\t', index=None, header=None)
	dups[3] = 0 ; dups[4] = 0
	dups.to_csv(os.path.join(outdir, "dups.annovar"), sep='\t', index=None, header=None)
	for tmp in glob.glob(os.path.join(outdir,"*.tmp")):
		os.remove(tmp)
	vars_del="annovar_input="+os.path.join(outdir,"dels.annovar")+",annovar_output="+os.path.join(outdir, sample+".dels.anno")
	vars_dup="annovar_input="+os.path.join(outdir,"dups.annovar")+",annovar_output="+os.path.join(outdir,sample+".dups.anno")
	print(vars_del)
	print(vars_dup)
	del_id = subprocess.check_output(["qsub", "-v", vars_del, "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/cnv_utils/annotate_cnvr.sh"]).rstrip().decode('utf-8')
	dup_id = subprocess.check_output(["qsub", "-v", vars_dup, "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/cnv_utils/annotate_cnvr.sh"]).rstrip().decode('utf-8')
	return([del_id, dup_id])

def combine_sample(outdir):
	if not os.path.isdir(os.path.join(outdir,"cnvnator_erds")):
		os.mkdir(os.path.join(outdir,"cnvnator_erds"))
	deps = []
	for subdirs, dirs, files in os.walk(os.path.join(outdir,"cnvnator")):
		for sample in dirs:	
			sampledir = os.path.join(outdir,"cnvnator_erds",sample)
			erds_del= os.path.join(outdir,"erds",sample, sample+".del.events")
			erds_dup=os.path.join(outdir,"erds",sample,sample+".dup.events")
			cnvnator_del= os.path.join(outdir,"cnvnator",sample,sample+".CNVnator.results.dels.bed")
			cnvnator_dup= os.path.join(outdir,"cnvnator",sample,sample+".CNVnator.results.dups.bed")
			if not os.path.isdir(sampledir):
				print(sample)
				os.mkdir(sampledir)
				sample_id = process_sample(sampledir, cnvnator_del, cnvnator_dup, erds_del, erds_dup)
				deps.append(sample_id)
	return(list(itertools.chain.from_iterable(deps)))

def preprocess(cnvdir):
	erds_dir =  os.path.join(cnvdir, "erds")
	cnvnator_dir = os.path.join(cnvdir, "cnvnator")
	preprocess_cnvnator(cnvnator_dir)
	deps = combine_sample(cnvdir)
	with open(os.path.join(cnvdir,"deps.tmp"), "w") as f:
                f.write(":".join(deps))
	combine(erds_dir, "*del.events", cnvdir, "erds.dels.tmp")
	combine(erds_dir, "*dup.events", cnvdir, "erds.dups.tmp")
	combine(cnvnator_dir, "*results.dels.bed", cnvdir, "cnvnator.dels.tmp")
	combine(cnvnator_dir, "*results.dups.bed", cnvdir, "cnvnator.dups.tmp")	

def main():

	parser = argparse.ArgumentParser(
		description='Process cnvs',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)

	parser.add_argument('-c','--cnvdir', required=False , dest='cnvdir', help='CNV directory containing CNVnator and ERDS calls')

	args = parser.parse_args()

	if args.cnvdir:
		preprocess(args.cnvdir)

main()


