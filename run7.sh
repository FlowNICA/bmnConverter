#!/bin/bash

#SBATCH -p fast
#SBATCH -t 01:00:00
#SBATCH -J bmn_converter
#SBATCH -o /home/ovgol/log/%A_%a.log
#SBATCH -a 1-1 

shft=4000

id=$((${SLURM_ARRAY_TASK_ID}+${shft}))
in=~/nica/bmndata/run7/plotnikov
out=~/nica/bmndata/run7/tree
outqa=~/nica/bmndata/run7/tree/qa
inFile=${in}/${id}.root
outFile=${out}/${id}.tree.root
qaFile=${outqa}/${id}.qa.root
mkdir -pv $out
mkdir -pv $outqa

. /home/ovgol/soft/bmnroot/build_run7_plotnikov/config.sh

time root -b -l -q $VMCWORKDIR/macro/run/bmnloadlibs.C \
	           /home/ovgol/soft/bmnConverter/convertBmn_run7.C"(\"${inFile}\", \"${inFile}\", \"${outFile}\")"
#time root -b -l -q /home/ovgol/soft/bmnConverter/qa.C"(\"${outFile}\", \"${qaFile}\")"
rm -v ${outFile}_tmp
