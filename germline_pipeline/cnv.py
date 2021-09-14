#!/usr/bin/python 
import subprocess
import pandas as pd
import glob 
import os

def format(file):
	cnv = pd.read_csv(file, sep ='\t', names=["CNV_type", "coordinates", "CNV_size", "normalized_RD", "e-val1", "e-val2", "e-val3", "e-val4", "q0"])
	cnv[["chrom","start"]]= cnv["coordinates"].str.split(":", expand=True)
	cnv[["start", "end"]] = cnv["start"].str.split("-", expand=True)
	cnv = cnv[["chrom", "start", "end", "CNV_type"]]
	return(cnv)

def preprocess(path):
	files = glob.glob(os.path.join(path, "*.CNVnator.results"))
	for file in files:
		cnvnator = format(file)
		cnvnator.to_csv(file + ".bed", sep='\t', index=None, header=None)
		dels = cnvnator[cnvnator.CNV_type == "deletion"]
		dels.to_csv(file + ".dels.bed", sep='\t', index=None, header=None)
		dups = cnvnator[cnvnator.CNV_type == "duplication"]
		dups.to_csv(file + ".dups.bed", sep='\t', index=None, header=None)

def run_cnvnator(bam, outdir, binsize):
	sample=os.path.basename(bam).split('.')[0]
	sampledir = os.path.join(outdir, sample)
	if glob.glob(os.path.join(sampledir,"*CNVnator.vcf")):
		return
	if not os.path.isdir(sampledir):
		os.mkdir(sampledir)
	print("Running CNVnator for {}".format(sample))
	os.chdir(sampledir)
	vars ="bam_file="+bam+",sampledir="+sampledir+",BINSIZE="+str(binsize)
	comm = ["qsub", "-v", vars, "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/cnv_utils/cnvnator.sh"]
	job_id=subprocess.check_output(comm).rstrip().decode('utf-8')	
	return(job_id)

def run_erds(bam, outdir, erdsdir, erds_deps):
	sample=os.path.basename(bam).split('.')[0]
	vcf = os.path.join(os.path.dirname(outdir), "gatk", sample, sample+".vcf")
	sampledir = os.path.join(erdsdir, sample)
	if glob.glob(os.path.join(sampledir,"*erds.vcf")):
		return
	if not os.path.isdir(sampledir):
		os.mkdir(sampledir)
	print("Running ERDS for {}".format(sample))
	os.chdir(sampledir)	
	vars ="erds_output="+sampledir+",bam_file="+bam+",vcf="+vcf
	if erds_deps:
		deps = "depend=afterok:" + ":".join(erds_deps)
		comm = ["qsub", "-W", deps,"-v", vars, "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/cnv_utils/erds.sh"]
	else:
		comm = ["qsub","-v", vars, "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/cnv_utils/erds.sh"]
	job_id=subprocess.check_output(comm).rstrip().decode('utf-8')
	return(job_id)

def cnv_pipeline(bams, outdir, binsize, erds_deps=None):
	os.chdir(outdir)
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	cnvnatordir = os.path.join(outdir, "cnvnator")
	if not os.path.isdir(cnvnatordir):
		os.mkdir(cnvnatordir)
	erdsdir = os.path.join(outdir, "erds")
	if not os.path.isdir(erdsdir):
		os.mkdir(erdsdir)
	if not os.path.isdir(os.path.join(outdir,"log")):
		os.mkdir(os.path.join(outdir,"log"))
	deps = []
	for bam in bams:
		cnv_id = run_cnvnator(bam, cnvnatordir, binsize)
		erds_id = run_erds(bam, outdir, erdsdir, erds_deps)
		if erds_id:
			deps.append(erds_id)
		if cnv_id:
			deps.append(cnv_id)
	vars="cnvnator_erds="+outdir
	if deps:
		alldeps="depend=afterok:" + ":".join(deps)
		comm=["qsub","-W", alldeps,"-v",vars, "-d", outdir, "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/cnv_utils/run_CNVanalysis.sh"]
	else:
		comm=["qsub","-v",vars, "-d", outdir, "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/cnv_utils/run_CNVanalysis.sh"]
	job_id=subprocess.check_output(comm).rstrip().decode('utf-8')	


