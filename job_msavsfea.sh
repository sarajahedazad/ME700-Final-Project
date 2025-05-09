#!/bin/bash -l
#$ -m aes
#$ -N myjob_msavsfea

module load miniconda
mamba activate me700-final
python3 results_MSAvsFEA_main.py
