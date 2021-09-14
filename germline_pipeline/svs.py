import os
import subprocess
import pandas as pd
import numpy as np
import sys
import glob

def create_bamtable(outdir, baminput):
	bamtable_path=os.path.join(outdir, "bam_table.tab")
	samples =  pd.Series([ os.path.basename(x).split(".")[0] for x in baminput])
	dummy = pd.Series(np.repeat("/hpf/largeprojects/adam/projects/lfs/lfs_germline/aligned_bams/NA12878/NA12878_20k.b37.bam", len(baminput), axis = 0))
	bamtable = pd.concat([samples,dummy,baminput], axis=1, keys = ["Names","Tumor_bams", "Normal_bams"])
	bamtable.to_csv(bamtable_path, sep='\t', index=False)	

def run_cfilter(bams, outdir, cfilter_deps=None):
	deps_list=[]
	types = ["del", "dup", "tra", "inv"]
	pype_dir = os.path.join(outdir,"pype_output")
	if not os.path.isdir(pype_dir):
		os.mkdir(pype_dir)
	for bam in bams:
		sample=os.path.basename(bam).split('.')[0]
		sampledir=os.path.join(outdir, "data", sample)
		pype_sample = os.path.join(pype_dir,sample)
		if not os.path.isdir(pype_sample):
			os.mkdir(pype_sample)
		if not os.path.isdir(os.path.join(pype_sample, "log")):
			os.mkdir(os.path.join(pype_sample, "log"))
		os.chdir(pype_sample)
		coverage_deps = cfilter_deps
		if not glob.glob(bam.replace(".bam","") +"*gatk.depthofcoverage*"):
			print("Need to create gatk coverage file for {}".format(sample))
			vars="bam="+bam+",sample="+sample
			if coverage_deps:
				comm=["qsub","-W", cfilter_deps,"-v", vars, "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/delly_utils/gatk_coverage.sh"]
				coverage_id=subprocess.check_output(comm).rstrip().decode('utf-8')
				coverage_deps = coverage_deps + ":"+coverage_id
			else:
				comm=["qsub","-v", vars, "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/delly_utils/gatk_coverage.sh"]
				coverage_id=subprocess.check_output(comm).rstrip().decode('utf-8')	
				coverage_deps = coverage_id
		for type in types:
			tab=os.path.join(sampledir,"normal_"+sample+"."+type+".tab")
			output=os.path.join(pype_sample,sample+"."+type+".Cfilter.tab")
			if os.path.isfile(output):
				continue
			vars="output="+output+",tab="+tab+",bam="+bam
			if coverage_deps:
				comm=["qsub","-W", coverage_deps,"-v", vars, "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/delly_utils/cfilter.sh"]
			else:
				comm=["qsub","-v", vars, "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/delly_utils/cfilter.sh"]
			job_id=subprocess.check_output(comm).rstrip().decode('utf-8')
			deps_list.append(job_id)
	deps="depend=afterany:" + ":".join(deps_list)
	return(deps)
			
def run_delly(baminput, outdir, yaml):
	dellydir=os.path.join(outdir,"delly")
	if not os.path.isdir(dellydir):
		os.mkdir(dellydir)
	os.chdir(dellydir)
	create_bamtable(dellydir, baminput)
	os.system("bash /hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/runDelly.sh")
	with open (os.path.join(dellydir,"cfilterdeps.txt"), "r") as myfile:
		deps=myfile.readlines()[0].rstrip('\n') + ":"
	return(deps)	

def run_manta(bams, outdir):
	mantadir=os.path.join(outdir,"manta")
	if not os.path.isdir(mantadir):
		os.mkdir(mantadir)
	os.chdir(mantadir)
	manta_deps=[]
	for bam in bams:
		cmd = ["bash", "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/setupWorkflowManta.sh",mantadir, bam]
		job_id=subprocess.check_output(cmd).rstrip().decode('utf-8')
		if job_id:
			manta_deps.append(job_id)
	manta_deps = ":".join(manta_deps)
	return(manta_deps)

def combine_delly_manta(bams,outdir,apply_deps):
	deps="depend=afterok:"
	if not os.path.isdir(os.path.join(outdir,"final_svs")):
		os.mkdir(os.path.join(outdir,"final_svs"))
	for bam in bams:
		vars="bam_file="+bam+",sv_dir="+outdir
		comm = [ "qsub","-W", apply_deps,"-v", vars, "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/combinePrefiltered.sh" ]
		job_id=subprocess.check_output(comm).rstrip().decode('utf-8')
		deps=deps+job_id+":"
	return(deps)

def sv_pipeline(bams, outdir, yaml, cfilteronly):
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	if not os.path.isdir(os.path.join(outdir,"log")):
                os.mkdir(os.path.join(outdir,"log"))
	delly_deps = run_delly(bams, outdir, yaml)
	manta_deps = run_manta(bams,outdir)
	all_deps = delly_deps + manta_deps
	os.chdir(outdir)
	deps=combine_delly_manta(bams,outdir,all_deps)
	vars ="svdir="+os.path.join(outdir, "final_svs")
	comm = ["qsub", "-W", deps, "-v", vars, "/hpf/largeprojects/adam/projects/lfs/code/germline_pipeline/sv_utils/runApplyFilter.sh"]
	job_id=subprocess.check_output(comm).rstrip().decode('utf-8')
	

