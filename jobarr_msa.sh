#!/bin/bash -l
#$ -t 1-10
#$ -N myjobarr_msa

module load miniconda
mamba activate me700-final
python3 MSA_main.py $SGE_TASK_ID
