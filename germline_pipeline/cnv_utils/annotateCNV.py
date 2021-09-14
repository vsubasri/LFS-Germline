import pandas as pd
import glob
import os
import itertools
import argparse 
import subprocess
import numpy as np
import scipy.stats

def check_overlap(x1, x2, y1, y2):
	overlap= max(0, min(x2, y2) - max(x1, y1))
	total=x2-x1
	return(overlap/total)

def mean_ci(data, confidence=0.95):
	a = 1.0 * np.array(data)
	n = len(a)
	m, se = np.mean(a), scipy.stats.sem(a)
	h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
	return(str(round(m,2)), " [" + str(round(m-h,2)) + "," + str(round(m+h,2)) + "]")

def check_cpg(gene):
	cp_genes_list= open("/hpf/largeprojects/adam/projects/lfs/resources/actionable_cpgs.txt", "r")
	cp_genes = cp_genes_list.read().split('\n')
	if gene in cp_genes:
		cpg='CPG'
	else:
		cpg='.'
	return(cpg)

def annotate_file(file_prefix):
	file = file_prefix + ".hg19_multianno.txt"
	if "dups" in file:
		type_ind = "observedgains"
	else:
		type_ind = "observedlosses"
	print("Annotating {}".format(file))
	dgv = pd.read_csv("/hpf/largeprojects/adam/projects/lfs/resources/DGV/GRCh37_hg19_variants_2016-05-15.txt", sep='\t')
	outfile= os.path.splitext(file)[0]+".final"
	with open(file,'r') as in_f, open(outfile, 'w') as out_f:		
		fl_indices=in_f.readline().strip().split("\t") + ['DGV_freq','DGV_freqfrac','DGVoverlap_mean','DGVoverlap_ci','CPG']
		firstline='\t'.join(fl_indices)+"\n"
		out_f.write(firstline)
		for line in in_f:
			line_indices = line.strip().split("\t")
			cpg = check_cpg(line_indices[6])
			acc_ids = line_indices[10].replace("Name=","").split(',')
			overlap_all = []
			if (acc_ids[0] != '.'):
				total_obs = 0
				for acc in acc_ids:
					entry = dgv[dgv['variantaccession'] == acc]
					overlap = check_overlap(int(line_indices[1]), int(line_indices[2]), int(entry['start']), int(entry['end']))
					overlap_all.append(overlap)
					total_obs = total_obs + entry[[type_ind]].values[0][0]
				mean, ci = mean_ci(overlap_all, confidence=0.95)
				freq_frac = str(round(total_obs/54946,2))
				newline = '\t'.join(line_indices + [str(total_obs), freq_frac, mean, ci, cpg]) + '\n'
			else:
				newline = '\t'.join(line_indices + ['.','.','.', '.', cpg]) + '\n'
			out_f.write(newline)

def generate_bedfiles(outdir):
	files = glob.glob(os.path.join(outdir, "*.annovar"))
	for file in files:
		sample=os.path.basename(file).split('.')[0]
		print("Generating bed file for {}...".format(sample))
		filtered = pd.read_csv(file, sep='\t')
		bed = cfiltered[cfiltered.columns[0:2]]
		if "del" in file:
			bed['type'] = "DEL"
		else:
			bed['type'] = "DUP"
		bed['sample'] = sample
		outfile = os.path.join(outdir,os.path.splitext(os.path.basename(file))[0]+".bed")
		bed.to_csv(outfile, sep='\t',index=None,header=None)

def main():

	parser = argparse.ArgumentParser(
		description='Process cnvs',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)

	parser.add_argument('-a','--annotate', required=False , dest='annotate', help='CNV ANNOVAR file to annotate')
	parser.add_argument('-b','--bedfile', required=False , dest='bedfile', help='directory with CNV annovar input files')

	args = parser.parse_args()

	if args.bedfile:
		generate_bedfiles(args.bedfile)	
	if args.annotate:
		annotate_file(args.annotate)	
main()


