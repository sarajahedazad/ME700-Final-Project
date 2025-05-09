#!/bin/bash -l
#$ -m b
#$ -N myjob_meshgen

module load miniconda
mamba activate me700-final
python3 archstruct_meshgen_main.py



