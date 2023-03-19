#!/bin/bash

#$ -wd /scratch1/mmamaev/data/DCM-QGSM-SMM/3.0
#$ -cwd
#$ -N Convert
#$ -q all.q
#$ -l h=!(ncx182.jinr.ru|ncx211.jinr.ru)
#$ -l h_rt=01:00:00
#$ -l s_rt=01:00:00
#
#$ -o /scratch1/mmamaev/data/DCM-QGSM-SMM/3.0/log
#$ -e /scratch1/mmamaev/data/DCM-QGSM-SMM/3.0/log

in=$1
shft=$2
[ -z $shft ] && shft=0
id=$SLURM_ARRAY_TASK_ID
out=/scratch1/mmamaev/data/DCM-QGSM-SMM/3.0
mkdir -pv $out
mkdir -pv $out/qa

. /scratch1/mmamaev/bmn_environment.sh

mkdir $out
time root -b -l -q /scratch1/mmamaev/bmnConverter/convertBmn.C"(\"${in}/${id}/dst2_geant_output.root\", \"${in}/${id}/geant_output.root\", \"${out}/$((${id}+${shft})).tree.root\")"
time root -b -l -q /scratch1/mmamaev/bmnConverter/qa.C"(\"${out}/$((${id}+${shft})).tree.root\", \"${out}/qa/$((${id}+${shft})).qa.root\")"
