#!/usr/bin/env python
import subprocess
import pandas as pd
import numpy as np
import os
import glob
import argparse
import re
from functools import reduce
from collections import Counter

def check_cpg(gene):
	cp_genes_list= open("/hpf/largeprojects/adam/projects/lfs/resources/actionable_cpgs.txt", "r")
	cp_genes = cp_genes_list.read().split('\n')
	if gene in cp_genes:
		cpg='CPG'
	else:
		cpg='.'
	return(cpg)

def pathogenic(annotdir):
	pd.options.mode.chained_assignment = None
	pathsvs_list = []
	files = glob.glob(os.path.join(annotdir, "*/*.annotated.final"))	
	pathsvs_all = [] ; filtered_all = []
	for file in files:
		print(file)
		outfile= os.path.splitext(file)[0]+".pathogenic"
		filtered_outfile = os.path.splitext(file)[0]+".fpfiltered"
		cpg_outfile=os.path.splitext(file)[0]+".pathogenic.cpg"
		allsvs = pd.read_csv(file, delimiter='\t')
		filteredsvs = allsvs[~allsvs['Panel_Status'].str.contains("GT3", na=False)]
		filteredsvs.columns.values[14] = "SAMPLE_INFO"
		filteredsvs['sample_id'] = os.path.basename(file).split('.')[0]
		filtered_all.append(filteredsvs)
		filteredsvs.to_csv(filtered_outfile,sep='\t', index=False)
		pathsvs = filteredsvs[filteredsvs['GD_POPMAX_AF'] < 0.001]
		pathsvs = pathsvs[~((pathsvs['SV type'] == "DEL") & (pathsvs['DGV_LOSS_Frequency'] > 0.01))]
		pathsvs = pathsvs[~((pathsvs['SV type'] == "DUP") & (pathsvs['DGV_GAIN_Frequency'] > 0.01))]
		pathsvs = pathsvs[pathsvs['AnnotSV ranking'] > 2]
		fullsvs = pathsvs[pathsvs['AnnotSV type'] == "full"] ; fullids = fullsvs['AnnotSV ID'].tolist(); 
		pathsvs = pathsvs[pathsvs['AnnotSV ID'].isin(fullids)]
		pathsvs_all.append(pathsvs)
		pathsvs_cpg = pathsvs[pathsvs['CPG'] == "CPG"]
		pathsvs.to_csv(outfile,sep='\t', index=False)
		pathsvs_cpg.to_csv(cpg_outfile,sep='\t', index=False)
	final_pathsvs = pd.concat(pathsvs_all, ignore_index=True)
	final_filtered = pd.concat(filtered_all, ignore_index=True)
	final_filtered_exonic = final_filtered[~final_filtered['CDS length'].isin([0,np.nan])]
	final_pathsvs_exonic = final_pathsvs[~final_pathsvs['CDS length'].isin([0,np.nan])]
	final_pathsvs.to_csv(os.path.join(annotdir,"all.pathogenic.txt"), sep='\t', index=False)
	final_filtered.to_csv(os.path.join(annotdir,"all.fpfiltered.txt"), sep='\t', index=False)
	final_filtered_cpg = final_filtered[final_filtered['CPG'] == "CPG"]
	final_filtered_cpg.to_csv(os.path.join(annotdir,"all.fpfiltered.cpg.txt"), sep='\t', index=False)
	final_pathsvs_exonic.to_csv(os.path.join(annotdir,"all.pathogenic.exonic.txt"), sep='\t', index=False)
	final_filtered_exonic.to_csv(os.path.join(annotdir,"all.fpfiltered.exonic.txt"), sep='\t', index=False)

def recurrent_svs(annotdir):
	gene_path = os.path.join(annotdir, "all.fpfiltered.txt")
	svs = pd.read_csv(gene_path, delimiter='\t', low_memory=False)
	svs = svs[svs['AnnotSV type'] == "split"]
	svs['ID'] = svs['Gene name'] + '_' + svs['SV type'] + '_' + svs['location']
	indivsvs = pd.melt(svs[['ID','sample_id']], id_vars='sample_id').drop(['variable'], axis=1).drop_duplicates(keep='first')
	ml_output = indivsvs.pivot_table(index='sample_id',columns='value',aggfunc=np.size).fillna(0)
	ml_output[ml_output != 0] = 1
#	byid= indivsvs.groupby(['value']).agg({'sample_id': lambda x: ','.join(x.unique())})
#	byid['ID'] = byid.index
	same_effect = ml_output.groupby(level=0).max().fillna(0)
	same_effect.drop([col for col, val in same_effect.sum().iteritems() if val < 2], axis=1, inplace=True)
	same_effect.to_csv(os.path.join(annotdir,'svs.fpfiltered.bygene.inputfile.txt'),sep='\t')
#	byid.to_csv(os.path.join(annotdir,'svs.pathogenic.recurrent_svs.txt'),sep='\t', index=False)

def check_cpg(gene):
        cp_genes_list= open("/hpf/largeprojects/adam/projects/lfs/resources/actionable_cpgs.txt", "r")
        cp_genes = [x for x in cp_genes_list.read().split('\n') if x]
        if gene in cp_genes:
                cpg='CPG'
        else:
                cpg='.'
        return(cpg)

def annotate_cpg(annotateddir):
        files = glob.glob(os.path.join(annotateddir, "*/*.annotated.tsv"))
        for file in files:
                outfile= os.path.splitext(file)[0]+".final"
                if not os.path.exists(outfile):
                        print(outfile)
                        with open(file,'r') as in_f, open(outfile, 'w') as out_f:
                                fl_indices=in_f.readline().strip().split("\t") + ['CPG']
                                firstline='\t'.join(fl_indices)+"\n"
                                out_f.write(firstline)
                                for line in in_f:
                                        line_indices = line.strip().split("\t")
                                        genes = line_indices[16].split('/')
                                        cpgs = []
                                        for gene in genes:
                                                cpg_tmp = check_cpg(gene)
                                                cpgs.append(cpg_tmp)
                                        if "CPG" in cpgs:
                                                newline = '\t'.join(line_indices + ['CPG']) + '\n'
                                        else:
                                                newline = '\t'.join(line_indices + ['.']) + '\n'
                                        out_f.write(newline)

def main():
	parser = argparse.ArgumentParser(
		description='Filter structural variants',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	)
	parser.add_argument("-s", '--svdir', action="store", dest="svdir", help='path to structural variants')

	args = parser.parse_args()
	annotate_cpg(args.svdir)
	pathogenic(args.svdir)
	recurrent_svs(args.svdir)

if __name__ == '__main__':
	main()

