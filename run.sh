#!/bin/bash

#SBATCH -p fast
#SBATCH -t 01:00:00
#SBATCH -J bmn_converter
#SBATCH -o /home/ovgol/log/%A_%a.log

in=$1
shft=$2
[ -z $shft ] && shft=0
id=$SLURM_ARRAY_TASK_ID
out=/mnt/pool/nica/7/ovgol/mc/bmnsim/XeCs
mkdir -pv $out
mkdir -pv $out/qa

. /mnt/pool/nica/7/mam2mih/soft/basov/bmnroot-mamaev/build/config.sh

mkdir $out
time root -b -l -q /home/ovgol/soft/bmn_converter/convertBmn.C"(\"${in}/${id}/dst2_geant_output.root\", \"${in}/${id}/geant_output.root\", \"${out}/$((${id}+${shft})).tree.root\")"
time root -b -l -q /home/ovgol/soft/bmn_converter/qa.C"(\"${out}/$((${id}+${shft})).tree.root\", \"${out}/qa/$((${id}+${shft})).qa.root\")"
