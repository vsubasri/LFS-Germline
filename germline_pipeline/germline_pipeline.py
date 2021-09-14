import pandas as pd
import argparse
import gatk
import svs
import cnv 
import os
import subprocess
import glob

def parse_bams(bams_input):
	if os.path.isfile(bams_input):
		bams=pd.read_csv(bams_input, sep='\n', header=None)[0]
	elif os.path.isdir(bams_input):
		bams=glob.glob(os.path.join(bams_input, "*.bam"))
	else:
		raise IOError("Wrong bam input")
	return(bams)

def main():
	
	parser = argparse.ArgumentParser(
		description='Run GATK Pipeline.',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)

	parser.add_argument('-o','--outdir', required=True , dest='outdir', help='Output destination for variants')
	parser.add_argument('-i','--intervallist', required=False, action='store_false', dest='intervallist', help='Use of interval list')
	parser.add_argument('-a','--allelespecific', required=False, action='store_true', dest='allelespecific', help='Use of interval list')
	parser.add_argument('-c', '--cfilteronly',required=False, action='store_true', dest='cfilteronly', help='Skip SV calling')
	parser.add_argument('-b','--bams', required=True, dest='bams', help='directory containing bams or line delimite file containing bams')
	parser.add_argument('-v','--v3', required=False, action='store_true', dest='v3', help='GATK Version')
	parser.add_argument('-s','--binsize', required=False, default="500", type=int, dest='binsize', help='CNVnator binsize')
	parser.add_argument('-y','--yaml', required=False, default="/hpf/largeprojects/adam/projects/lfs/lfs_germline/delly/lfs.yaml", dest='yaml', help='path to yaml for delly')
	parser.add_argument('--sv', dest='sv', action='store_true')
	parser.add_argument('--snv', dest='snv', action='store_true')
	parser.add_argument('--cnv', dest='cnv', action='store_true')
	
	args = parser.parse_args()

	bams = parse_bams(args.bams)
	if args.snv or (args.cnv and not os.path.exists(os.path.join(args.outdir,"gatk"))):
		erds_deps = gatk.gatk(bams, os.path.join(args.outdir,"gatk"), args.intervallist, args.allelespecific, args.v3)
		if args.cnv:
			cnv.cnv_pipeline(bams, os.path.join(args.outdir, "cnvs"), args.binsize, erds_deps)
	if args.cnv and os.path.exists(os.path.join(args.outdir,"gatk")):
		cnv.cnv_pipeline(bams, os.path.join(args.outdir, "cnvs"), args.binsize)
	if args.sv:
		svs.sv_pipeline(bams, os.path.join(args.outdir, "structural_variations"), args.yaml, args.cfilteronly)

main()

