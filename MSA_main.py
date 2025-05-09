# importing modules
from MSA_demo import *
from reading_saving_functions import *
import sys

# reading: params
parent_dir_inputs = 'inputs/archstruct'
jobparams_csv_folder_dir = os.path.join( parent_dir_inputs, 'jobparams_csv')


''' Parameters that one can read from a csv input file '''
jobtype = 'samplecheck'
perturb_type = 'nan'
job_id = int(sys.argv[1]) - 1
#job_id = 3
mesh_folder_root = 'inputs/archstruct/mesh/mesh2D'
outputs_root_dir = 'outputs/archstruct/results2D/MSA'

params_dict = get_params_dict( job_id, jobparams_csv_folder_dir , jobtype, perturb_type )
frame = get_frame_from_params( mesh_folder_root, params_dict)
structure_length = params_dict['ny'] * params_dict['cell_size']

# disp critical, eigenvector
disp_critical, eigenvector = calc_dispcritical_eigmode_fixed_to_displaced_end( frame, structure_length)

# force critical
force_critical = calc_forcecritical_fixed_to_displaced_end( frame, structure_length, disp_critical )

# save disp critical, force critical, eigenvector
save_MSA_results( outputs_root_dir, params_dict, disp_critical, force_critical, eigenvector)

