#!/bin/bash

module purge; module load python/3.4.0

YAML=/hpf/largeprojects/adam/projects/lfs/lfs_germline/structural_variations/delly/lfs.yaml

python /hpf/largeprojects/adam/local/sw/pype/2.3.1/pypegermline.py --c $YAML --tab bam_table.tab --delly

pwd

