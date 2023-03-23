#!/bin/bash

#$ -wd /scratch1/mmamaev/data/DCM-QGSM-SMM/3.0
#$ -cwd
#$ -N Convert
#$ -q all.q
#$ -l h=!(ncx182.jinr.ru|ncx211.jinr.ru)
#$ -l h_rt=01:00:00
#$ -l s_rt=01:00:00
#
#$ -o /scratch1/mmamaev/data/log
#$ -e /scratch1/mmamaev/data/log

in=$1
out=$2
id=$SGE_TASK_ID

mkdir -pv $out
#mkdir -pv $out/qa

source /cvmfs/nica.jinr.ru/sw/os/login.sh
module add GCC-Toolchain/
source /scratch1/mmamaev/bmn_environment.sh

time root -b -l -q /scratch1/mmamaev/bmnConverter/convertBmn.C"(\"${in}/${id}/dst_geant_output.root\", \"${in}/${id}/geant_output.root\", \"${in}/${id}/full_geometry .root\", \"${out}/${id}.tree.root\")"
#time root -b -l -q /scratch1/mmamaev/bmnConverter/qa.C"(\"${out}/${id}.tree.root\", \"${out}/qa/${id}.qa.root\")"
