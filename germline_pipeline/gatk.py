import argparse 
import os
import subprocess
import glob
import numpy as np
from math import floor
import csv

def parse_version(v3):
        steps=["haplotypecaller", "genotypegvcfs", "recalibration", "filter"]
        if v3:
                scripts= ["/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/haplotype_caller3.sh",
                "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/genotype_gvcfs3.sh",
                "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/recalibration3.sh"]
        else:
                scripts= ["/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/haplotype_caller4.sh",
                "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/genotype_gvcfs4.sh",
                "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/recalibration4.sh"]
        scripts.append("/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/run_variantfiltering.sh")
        script_dict=dict(zip(steps, scripts))
        return(script_dict)

def haplotypecaller(bams, outdir, intervallist, allelespecific, script):
	dependencies=[]
	for bam in bams:
		sample=os.path.basename(bam).split('.')[0]
		sample_dir=os.path.join(outdir, sample)
		if not os.path.isdir(sample_dir):
			os.makedirs(sample_dir)
			os.makedirs(os.path.join(sample_dir,"log"))
		os.chdir(sample_dir)
		if not os.path.isfile(os.path.join(sample_dir,sample+".g.vcf")):
			print("Running Haplotype Caller for {}".format(sample))
			vars="bam="+bam+",INTERVALLIST="+str(intervallist)+",ALLELESPECIFIC="+str(allelespecific)
			comm=["qsub","-v",vars,"-d",sample_dir, script]
			job_id=subprocess.check_output(comm).rstrip().decode('utf-8')
			dependencies.append(job_id)
	return(dependencies)

def parse_controls(n, outdir):
	controls = []
	with open('/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/gatk_utils/control.list','r') as f:
		reader=csv.reader(f,delimiter='\n')
		for control in reader:
			sample_dir= os.path.join(outdir,os.path.basename(control[0]).rstrip(".g.vcf"))
			if not os.path.isdir(sample_dir):
				os.mkdir(sample_dir)
			if len(controls) == n:
				return(controls)
			controls.extend(control)

def genotypegvcfs(geno_dependencies, bams, outdir, intervallist, allelespecific, script):
	control_num=0
	##find optimal number to combine
	batches = floor(len(bams)/22)
	remainder= len(bams) % 22
	if (remainder > 12):
		batch_size = np.concatenate([np.repeat(22, batches), [remainder]])
	elif ((remainder < 12) and (len(bams) < 12)):
		print("WARNING: less than 13 samples running at one time, running with LMS controls for increased accuracy of calls")
		control_num = 13-remainder
		control_gvcfs = parse_controls(control_num,outdir)
		batch_size = np.array([remainder])
	else:
		alloc = floor(remainder/batches) + 22
		alloc_rem = remainder % batches + 22
		batch_size = np.concatenate([np.repeat(alloc, batches),[alloc_rem]])
	##call script
	deps="depend=afterok:" + ":".join(geno_dependencies)
	geno_paths = []
	input_lists = []
	bam_lists = []
	job_ids= []
	for idx, batch in enumerate(batch_size):
		batch_bams = bams[:batch+1]
		bams = bams[batch+1:]
		idx=idx+1
		inputlist=os.path.join(outdir, "input.batch"+str(idx)+".list")
		gvcf_list = [ os.path.join(outdir,os.path.basename(x).rstrip(".realigned-recalibrated.bam"), os.path.basename(x).rstrip(".realigned-recalibrated.bam") + ".g.vcf") for x in batch_bams]
		if control_num != 0:
			gvcf_list = gvcf_list + control_gvcfs
		bamlist=os.path.join(outdir, "bam.batch"+str(idx)+".list")
		with open(bamlist, 'w', newline='') as f_bam:
			tsv_output = csv.writer(f_bam, delimiter='\n')
			tsv_output.writerow(batch_bams)
		bam_lists.append(bamlist)
		geno_paths.append(os.path.join(outdir, "genotyped_variants."+str(idx)+".g.vcf"))
		with open(inputlist, 'w', newline='') as f_output:
			tsv_output = csv.writer(f_output, delimiter='\n')
			tsv_output.writerow(gvcf_list)
		input_lists.append(inputlist)
		os.chdir(outdir)
		vars="INPUTLIST="+inputlist+",INTERVALLIST="+str(intervallist)+",ALLELESPECIFIC="+str(allelespecific)
		comm=["qsub","-W",deps,"-v",vars,"-d",outdir, script]
		if deps == "depend=afterok:":
			comm=["qsub","-v",vars,"-d",outdir, script]
		job_id=subprocess.check_output(comm).rstrip().decode('utf-8')
		job_ids.append(job_id)
		if len(bams) == 0:
			break
	return(job_ids, geno_paths, input_lists, bam_lists)

def recalibrate(recal_dependencies, geno_paths, outdir, intervallist, allelespecific, script):
	filter_deps = []
	recal_paths = []
	os.chdir(outdir)
	deps="depend=afterok:" + ':'.join(recal_dependencies)
	for path in geno_paths:
		batch = path.split('.')[-3]
		recal_paths.append(os.path.join(outdir,"indel.recalibrated."+batch+".vcf"))
		vars="INPUTFILE="+path+",INTERVALLIST="+str(intervallist)+",ALLELESPECIFIC="+str(allelespecific)
		comm= ["qsub","-W",deps, "-v",vars,"-d",outdir,script]
		job_id=subprocess.check_output(comm).rstrip().decode('utf-8')
		filter_deps.append(job_id)
	return(filter_deps, recal_paths)

def gatk(bams, outdir, intervallist, allelespecific, v3):
	scripts = parse_version(v3)
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	os.chdir(outdir)
	if not os.path.isdir(os.path.join(outdir, "log")):
		os.mkdir(os.path.join(outdir, "log"))
	geno_deps = haplotypecaller(bams, outdir, intervallist, allelespecific, scripts['haplotypecaller'])
	recal_deps, geno_paths, input_lists, bam_lists  = genotypegvcfs(geno_deps, bams, outdir, intervallist, allelespecific, scripts['genotypegvcfs'])
	filter_deps, recal_paths = recalibrate(recal_deps, geno_paths, outdir, intervallist, allelespecific, scripts['recalibration'])
	erds_deps = run_filtering(filter_deps, recal_paths, bam_lists,  outdir, scripts["filter"])
	return(erds_deps)

def run_filtering(filter_dependencies, recal_paths, inputlists, outdir, script):
	deps="depend=afterok:" + ":".join(filter_dependencies)
	erds_deps = []
	for inputlist, recal_path in zip(inputlists, recal_paths):
		vars="GATKDIR="+outdir+",BAMDIR="+inputlist+",MULTISAMPLE="+recal_path
		comm=["qsub","-W",deps,"-v",vars,"-d",outdir,script]
		job_id=subprocess.check_output(comm).rstrip().decode('utf-8')
		erds_deps.append(job_id)	
	return(erds_deps)

