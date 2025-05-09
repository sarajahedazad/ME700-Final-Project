#!/bin/bash -l
#$ -t 1-10
#$ -N myjobarr_fea

module load miniconda
mamba activate me700-final
python3 FEA_main.py $SGE_TASK_ID
