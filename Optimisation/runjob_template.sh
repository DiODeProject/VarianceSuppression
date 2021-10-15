#!/bin/bash
#$ -N jobname
#$ -l mem=4G 
#$ -l rmem=4G
#$ -l h_rt=hrs:min:00
#$ -e errors/
#$ -o outputs/

module load apps/python/conda
source activate python35
 
python script_scipy_optimise.py -t ${1} -z ${2} -r ${3} -o ${4}

