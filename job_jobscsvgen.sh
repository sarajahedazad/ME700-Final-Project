#!/bin/bash -l
#$ -N myjob_jobscsvgen

module load miniconda
mamba activate me700-final
python3 archstruct_jobscsvgen_main.py

