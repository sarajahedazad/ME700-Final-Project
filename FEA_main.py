from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import pickle
import h5py
import time
import pandas as pd
from reading_saving_functions import save_FEA_results, get_params_dict
# Dealing with ufl legacy
import sys
try:
    from ufl import diag, Jacobian, shape
except:
    from ufl_legacy import diag, Jacobian, shape

from FEA_functions import *

#---------input--------------

parent_dir_inputs = 'inputs/archstruct'
jobparams_csv_folder_dir = os.path.join( parent_dir_inputs, 'jobparams_csv')


''' Parameters that one can read from a csv input file '''
jobtype = 'samplecheck'
perturb_type = 'nan'
job_id = int(sys.argv[1]) - 1
#job_id = 3

params_dict = get_params_dict( job_id, jobparams_csv_folder_dir , jobtype, perturb_type )


#---------running fea-------------

results_dict, _ = run_fea( perturb_type, params_dict )

#-------------output--------------
outputs_root_dir = 'outputs/archstruct/results2D/FEA'
save_FEA_results( outputs_root_dir, params_dict, results_dict)