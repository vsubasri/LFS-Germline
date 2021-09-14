#!/usr/bin/env python

import glob
import argparse
import csv
import pandas as pd
import numpy as np
import os
import subprocess
import re 
import io
import math
from optparse import OptionParser
from functools import reduce

pd.options.mode.chained_assignment = None  # default='warn'

def read_largefile(fname, chunk_size=10000):
    reader = pd.read_csv(fname, header=0, delimiter='\t',low_memory=False, iterator=True)
    chunks = []
    loop = True
    while loop:
        try:
            chunk = reader.get_chunk(chunk_size)
            chunks.append(chunk)
        except StopIteration:
            loop = False
            print("Iteration is stopped")
    df_ac = pd.concat(chunks, ignore_index=True)
    return(df_ac)

def split_chrom(in_filename):
	command=['awk \'NR==1{h=$0; next};!seen[$1]++{f=$1".tsv"; print h > f};{f=$1".tsv"; print >> f; close(f)}\'', in_filename]
	subprocess.call(command, shell=True)

def filter_all(in_filename,out_filename):
	genotypes=read_largefile(geno_filename)
	with open(in_filename, 'r') as in_f, open(out_filename, 'w') as out_f:
		indices=in_f.readline().strip().split("\t")
		first_line=pd.concat((pd.DataFrame([indices]),pd.DataFrame([genotypes.columns[7:]])), axis=1)
		s1='\t'.join(str(v) for v in first_line.ix[0]) + '\n'
		out_f.write(s1)
		gnomad_genome_index=indices.index('gnomAD_genome_ALL')
		thoug_index=indices.index('1000g2015aug_all')
		exac_index=indices.index('ExAC_ALL')
		gnomad_exome_index=indices.index('gnomAD_exome_ALL')
		exonic_index=indices.index('Func.refgGene')
		func_index=indices.index('ExonicFunc.refGene')
		for line in in_f:
			funcs = line.strip().split("\t")
			elements = pd.to_numeric(line.strip().split("\t"),errors='coerce')
			if ((funcs[exonic_index] == 'exonic') & 
				(funcs[func_index] == 'nonsynonymous SNV') & 
				((elements[gnomad_genome_index] < 0.01) | (np.isnan(elements[gnomad_genome_index]))) & 
				((elements[thoug_index] < 0.01) | (np.isnan(elements[thoug_index]))) & 
				((elements[exac_index] < 0.01) | (np.isnan(elements[exac_index]))) &
				((elements[gnomad_exome_index] < 0.01) | (np.isnan(elements[gnomad_exome_index])))):
				geno=genotypes.loc[(genotypes['CHROM'] == elements[0]) & (genotypes['POS'] == elements[1])]
				elements=pd.DataFrame([line.strip().split("\t")])
				new_line=pd.concat((elements, geno[geno.columns[7:]].reset_index(drop=True)),axis=1)
				new_line[new_line.columns[0:3]] = new_line[new_line.columns[0:3]].astype(np.int64)
				s2='\t'.join(str(v) for v in new_line.ix[0]) + '\n'
				out_f.write(s2)

def combine_pathogenic(pathdir, path_outfile, vus):
	if vus:
		path_outfiles=glob.glob(os.path.join(pathdir,"*.hg19.pathogenic.vus.txt"))
	else:
		path_outfiles=glob.glob(os.path.join(pathdir,"*.hg19.pathogenic.txt"))
	for s in path_outfiles:
		s_intervar = read_largefile(s)
		sample=os.path.basename(s).split('.')[5]
		#sample=os.path.basename(s).split('.')[0]
		s_intervar["sample"] = sample
		if not os.path.isfile(path_outfile):
			with open(path_outfile, 'w') as f:
				f.write('\t'.join(s_intervar.columns) + '\n')		
		print('Concatenating sample {}...'.format(sample))
		s_intervar.to_csv(path_outfile, sep='\t', index=False, header=False, mode='a')

def filter_func(intervar_file):
	func_outfile = intervar_file + ".ingene"
	with open(intervar_file, 'r') as in_f, open(func_outfile, 'w') as out_f:
		firstline=in_f.readline()
		indices=firstline.strip().split("\t")
		out_f.write(firstline)
		for line in in_f:
			line_indices = line.strip().split("\t")
			if (line_indices[7] not in ["intergenic","upstream","downstream"]) and any(x in line_indices[17] for x in ["Likely pathogenic", "Pathogenic", "Uncertain significance"]):
				 out_f.write(line)

def pathogenic_multisample(intervar_dir, vus):
	m_intervar = os.path.join(intervar_dir,"*.hg19_multianno.txt.intervar")
	intervar_files = glob.glob(m_intervar)
	path_outfiles = []
	for file in intervar_files:
		print("Filtering {}".format(file))
		int_split = file
#		int_split = split_intervar(file)
		summ_outfile, path_outfile = find_pathogenic(int_split, vus)
		path_outfiles.append(path_outfile)
	
	if vus:
		path_outfile = os.path.join(intervar_dir,"all.hg19.pathogenic.vus.txt")
		summ_outfile = os.path.join(intervar_dir,"all.hg19.pathogenic.vus.summary.txt")
		multihits_outfile = os.path.join(intervar_dir,"all.hg19.pathogenic.vus.multiplehits.txt")
	else:
		path_outfile = os.path.join(intervar_dir,"all.hg19.pathogenic.txt")
		summ_outfile = os.path.join(intervar_dir,"all.hg19.pathogenic.summary.txt")
		multihits_outfile = os.path.join(intervar_dir,"all.hg19.pathogenic.multiplehits.txt")

	combine_pathogenic(intervar_dir, path_outfile, vus)

	path_intervar = pd.read_csv(path_outfile, delimiter='\t', header=0, dtype={'Start': str, 'sample': str})
	grouped = path_intervar[["#Chr","Start","Ref.Gene", "sample"]].groupby("Ref.Gene", as_index=False).agg({'#Chr':'first','Start':[lambda y: ','.join(y.unique())],'sample': ['count',lambda x: ','.join(x.unique())]})

	multiple_hits(path_intervar, multihits_outfile)

	cp_genes_list= pd.read_csv("/hpf/largeprojects/adam/projects/lfs/resources/CPG_NEJM1508054.csv")
	cp_genes = cp_genes_list['Gene']
	grouped.reset_index()
	grouped.columns = ['Gene', 'Count', 'Samples', 'Chr', 'Pos']
	grouped["cp_genes"] = np.where(grouped["Gene"].isin(cp_genes), "T", "F")
	grouped.to_csv(summ_outfile, sep='\t', index=False) 
	return(summ_outfile)

def find_pathogenic(m_path, vus):
	if (vus):
		search_status = ["Likely pathogenic", "Pathogenic", "Uncertain significance"]
		path_outfile = os.path.join(os.path.dirname(m_path), m_path.replace(".hg19_multianno.txt.intervar","") + ".hg19.pathogenic.vus.txt")
		summ_outfile =  os.path.join(os.path.dirname(m_path), m_path.replace(".hg19_multianno.txt.intervar","") + ".hg19.pathogenic.vus.summary.txt") 
	else:
		search_status = ["Likely pathogenic", "Pathogenic"] 
		path_outfile =  os.path.join(os.path.dirname(m_path), m_path.replace(".hg19_multianno.txt.intervar","") + ".hg19.pathogenic.txt")
		summ_outfile =  os.path.join(os.path.dirname(m_path), m_path.replace(".hg19_multianno.txt.intervar","") + ".hg19.pathogenic.summary.txt")
	if os.path.isfile(path_outfile):
		return(summ_outfile,path_outfile)
	with open(m_path, 'r') as in_f, open(path_outfile, 'w') as out_f:
		firstline=in_f.readline()
		indices=firstline.strip().split("\t")
		out_f.write(firstline)
		for line in in_f:
			line_indices = line.strip().split("\t")
			if any(x in line_indices[17] for x in search_status) and (not any(x in line_indices[13] for x in ["Benign","Likely benign"])):
				out_f.write(line)
	return(summ_outfile,path_outfile)

def separate_snvs_indels(path_outfile):
	path_intervar = pd.read_csv(path_outfile, delimiter='\t', header=0)
	snvs =  path_intervar[~(path_intervar.Ref.str.contains('-') | path_intervar.Alt.str.contains('-'))]
	snvs.to_csv(os.path.splitext(path_outfile)[0] + ".snvs.txt", sep='\t', index=False)
	indels = path_intervar[(path_intervar.Ref.str.contains('-') | path_intervar.Alt.str.contains('-'))]
	indels.to_csv(os.path.splitext(path_outfile)[0] + ".indels.txt", sep='\t', index=False)

def multiple_hits(path_intervar, multihits_outfile):
	if not path_intervar.filter(like='sample').empty:
		multiplehits = path_intervar[["#Chr","Start","Ref.Gene", "sample"]].groupby(["Ref.Gene", "sample"], as_index=False).agg({'#Chr':'first','Start':[lambda y: ','.join(y.unique()), 'count']})
		cols = ['Gene','Sample', 'Pos', 'Count','Chr']
	else:
		multiplehits = path_intervar[["#Chr","Start","Ref.Gene"]].groupby(["Ref.Gene"], as_index=False).agg({'#Chr':'first','Start':[lambda y: ','.join(y.unique()), 'count']})
		cols = ['Gene','Pos', 'Count','Chr']
	multiplehits = multiplehits[multiplehits['Start']['count'] > 1]
	multiplehits.reset_index()
	multiplehits.columns = cols
	multiplehits.to_csv(multihits_outfile, sep='\t', index=False)

def pathogenic_singlesample(int_presplit, vus):
	m_path = int_presplit
#	m_path = split_intervar(int_presplit)
	summ_outfile, path_outfile = find_pathogenic(m_path, vus)
	path_intervar = pd.read_csv(path_outfile, delimiter='\t', header=0, dtype={'Start': str, 'sample': str})
	## summary grouped by gene 
	grouped = path_intervar[["#Chr","Start","Ref.Gene"]].groupby("Ref.Gene", as_index=False).agg({'#Chr':'first','Start':[lambda y: ','.join(y.unique())]})
	multiple_hits(path_intervar, summ_outfile)
	cp_genes_list= pd.read_csv("/hpf/largeprojects/adam/projects/lfs/resources/CPG_NEJM1508054.csv")
	cp_genes = cp_genes_list['Gene']
	grouped.reset_index()
	grouped.columns = ['Gene','Chr', 'Pos']
	grouped["cp_genes"] = np.where(grouped["Gene"].isin(cp_genes), "", "F")
	grouped.to_csv(summ_outfile, sep='\t', index=False)  
	return(summ_outfile)

def extract_genes(m_path, pathtogenes,ext):
	cpgenes_list= open(pathtogenes, "r")
	cpgenes = cpgenes_list.read().split('\n')
	intervar_files = glob.iglob(os.path.join(m_path, "*"+ext))
	cpgenes_outfiles = []
	for file in intervar_files:
		cpgenes_outfile = os.path.splitext(file)[0]+".genes"
		print("Writing file {} ...".format(cpgenes_outfile))
		cpgenes_outfiles.append(cpgenes_outfile) 
		if not os.path.exists(cpgenes_outfile):
			with open(file, 'r') as in_f, open(cpgenes_outfile, 'w') as out_f:
				firstline=in_f.readline()
				indices=firstline.strip().split("\t")
				out_f.write(firstline)
				for line in in_f:
					line_indices = line.strip().split("\t")
					if (line_indices[5] in cpgenes) and (line_indices[7] not in ["intergenic", "upstream", "downstream"]) and not any(x in line_indices[17] for x in ["Benign", "Likely benign"]):
						out_f.write(line)
	cpgenes_multisample = m_path.split('*')[0] +".genes"
	## multisample
	if not os.path.isfile(m_path):
		for s in cpgenes_outfiles:
			s_cpgenes = read_largefile(s)
			sample=s.split('.hg19')[0].split('.')[-1]
			s_cpgenes["sample"] = sample
			if not os.path.isfile(cpgenes_multisample):
				with open(cpgenes_multisample, 'w') as f:
					f.write('\t'.join(s_cpgenes.columns) + '\n')
			with open(cpgenes_multisample, 'a') as f:
				s_cpgenes.to_csv(f, sep='\t', index=False, header=False)

def extract_vaf(pathfile, input_vcf):
	os.system("""sed '1d' """ + pathfile +""" | cut -f1,2,3 | bedtools sort  > """ + os.path.splitext(pathfile)[0] + ".bed")
	vars="input="+input_vcf+",bed="+ os.path.splitext(pathfile)[0] + ".bed"+",output="+ os.path.splitext(pathfile)[0]+".vcf"
	submit=["qsub","-v",vars, "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/extract_vaf.sh"]
	job_id=subprocess.check_output(submit).rstrip().decode('utf-8')

def clinical_analyses(summ_outfile, clinical_file):
	summ = pd.read_csv(summ_outfile, delimiter='\t', header=0)
	clinical = pd.read_csv(clinical_file, delimiter='\t', header=0)
	tt = dict(zip(clinical['sample'], clinical['tumortype']))
	summ['tumortype'] = summ['Samples'].replace(tt, regex=True)
	ageofonset = dict(zip(clinical['sample'], clinical['ageofonset'].astype(str)))
	summ['ageofonset'] = summ['Samples'].replace(ageofonset, regex=True).str.split(',').apply(lambda x: [ i for i in np.array(x, dtype=np.float64) if not math.isnan(i) ])
	summ['mean_ageofonset'] = [ np.round(np.mean(x),decimals=2) if len(x) != 0 else np.nan for x in summ['ageofonset'] ]
	summ['sd_ageofonset'] = [ np.round(np.std(x),decimals=2) if len(x) != 0 else np.nan for x in summ['ageofonset'] ]
	summ.to_csv(os.path.splitext(summ_outfile)[0] + ".clinical.summary.txt", sep='\t', index=False)

def summarize_snvs(input):
	outfile=os.path.splitext(input)[0] + ".input"
	with open(input,'r') as in_f, open(outfile, 'w') as out_f:
		out_f.write(in_f.readline())
		for line in in_f:
			line_indices = line.strip().split("\t")
			ref = line_indices[2]
			alts = line_indices[3].split(',')
			samples = line_indices[7:] 
			samples_gts = []
			for sample in samples:
				sample_alleles = sample.split('/')
				if set(sample_alleles) & set(alts):
					gt = "1"
				else: 
					gt = "0"
				samples_gts.append(gt)
			newline = '\t'.join(line_indices[0:7] + samples_gts) + '\n'
			out_f.write(newline)

def split_intervar(file):
	outfile = file + ".final"
	if os.path.exists(outfile):
		return(outfile)
	else:
		with open(file, 'r') as in_f, open(outfile, 'w') as out_f:
			firstline=in_f.readline()
			indices=firstline.strip().split("\t") + ['genotype','allele_depth','depth_of_coverage','genotype_quality','phred_likelihood']
			out_f.write('\t'.join(indices)+'\n')
			for line in in_f:
				linesplit = line.rstrip().split('\t')
				fi = linesplit[45].split(':')
				if len(fi) > 5:
					fi = fi[1:4] + [fi[len(fi)-1]]
				linesplit = linesplit[0:44]+fi
				out_f.write('\t'.join(linesplit)+'\n')
		return(outfile)

def main():
	parser = argparse.ArgumentParser(
		description='Filter variants',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	)
	parser.add_argument("-m", '--summarize', action="store", dest='summarize', help='file/directory of single sample gatk tab files to summarize into bed files', required=False)
	parser.add_argument("-v", '--vcf', action="store", dest='vcf', help='vcf containing all variants', required=False)
	parser.add_argument('--vus', action="store_true", dest='vus', help='whether to include vus in pathogenic file', required=False)
	parser.add_argument("-i", '--intervar', action="store", dest="intervar", help='single intervar file or directory containing multiple intervar files', required=False)
	parser.add_argument("-c", '--clinical', action="store", dest="clinical", help='clinical file', required=False)
	parser.add_argument("-g", '--genes', action="store", dest="genes", help='extract variants from list of genes', required=False)
	parser.add_argument("-b", '--bed', action="store", dest="bed", help='bed file with list of coordinates to extract minor allele frequency', required=False)
	parser.add_argument("-s", '--separate', action="store", dest="separate", help='separate snvs and indels', required=False)
	parser.add_argument("-f", '--func', action="store", dest="func", help='intervar file to filter for ingene/splicing variants', required=False)	
	parser.add_argument("-e", '--ext', action="store", dest="ext", help='extension for amalgamating files',required=False)

	args = parser.parse_args() 

	## identify pathogenic variants from intervar output
#	if args.intervar:
#		if os.path.isfile(args.intervar):
#			summ_outfile = pathogenic_singlesample(args.intervar, args.vus)
#		else:
#			summ_outfile = pathogenic_multisample(args.intervar, args.vus)

	## separate snvs and indels -- must be tab file with "Ref" and "Alt" columns
	if args.separate:
		separate_snvs_indels(args.separate)

	## create summarized ml-input formatted file of snvs/indels
	if args.summarize:
		summarize_snvs(args.summarize)

	## extract variants for select genes from intervar file 
	if args.genes and args.intervar and args.ext:
		extract_genes(args.intervar, args.genes, args.ext)

	## extract vaf of variants at select coordinates 
	if args.bed and args.vcf:
		extract_vaf(args.bed, args.vcf)
	
	if args.func:
		filter_func(args.func)

if __name__ == '__main__':
	main()
