import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import h5py
import pickle
import time
import pandas as pd
import math

def get_params_dict(job_id, folder_path, jobtype, perturb_type):
    """
    Load all parameters from CSV into a single dict.
    CSV must contain columns for variable, fixed, and solver settings.
    """
    csv_file = os.path.join(folder_path, f"jobs_{jobtype}_perturb{perturb_type}.csv")
    df = pd.read_csv(csv_file)
    params = df.iloc[job_id].to_dict()
    # normalize types
    params['strkey'] = str(int(params['strkey']))
    params['ID']     = int(params['ID'])
    return params

def get_params_dict_from_df(job_id, df):
    params = df.iloc[job_id].to_dict()
    # normalize types
    params['strkey'] = str(int(params['strkey']))
    params['ID']     = int(params['ID'])
    return params

def get_name_extensions(perturb_type, params):
    # variable parameters
    strkey        = params['strkey']
    ID            = params['ID']
    a             = params['a']
    lcar          = params['lcar']
    tstep         = params['tstep']
    pf_up, pf_dn  = params['pf_force_up'], params['pf_force_down']
    H_up, H_dn    = params['H_force_up'],  params['H_force_down']
    pf_disp       = params['pf_dispfield']

    # fixed + solver settings
    cell_size      = params['cell_size']
    cellvalues_str = params['cellvalues_str']
    ny, nx         = params['ny'], params['nx']
    E, nu          = params['E'], params['nu']
    stop_thresh    = params['stop_thresh']
    rel_param      = params['rel_param']
    abs_tol, rel_tol = params['abs_tol'], params['rel_tol']
    max_bs           = params['max_bs_iter_cnt']
    find_load        = params['find_critical_load']

    keys_ext = f'ny{ny}nx{nx}celltypes{cellvalues_str}'
    mesh_ext = f'{keys_ext}cellsize{cell_size}lcar{lcar}'
    mat_ext  = f'r{a}'
    sol_ext  = f'step{tstep}thresh{stop_thresh}abstol{abs_tol}reltol{rel_tol}'
    pert_ext = f'perturbtype{perturb_type}forcepfup{pf_up}forcepfdown{pf_dn}locup{H_up}locdown{H_dn}disppf{pf_disp}'

    save_ext = f'{mesh_ext}{mat_ext}_{sol_ext}_{pert_ext}'
    save_shortened_ext = f'lcar{lcar}_{sol_ext}_{pert_ext}'

    return {
        'keys': keys_ext,
        'mesh': mesh_ext,
        'mat':  mat_ext,
        'solve': sol_ext,
        'perturb': pert_ext,
        'save': save_ext,
        'save_shortened_ext': save_shortened_ext
    }

def save_MSA_dispcritical( outputs_root_dir, params, disp_critical):
    saving_folder_dir = os.path.join( outputs_root_dir, f'resultslcar{params['lcar']}' )
    os.makedirs( saving_folder_dir, exist_ok = True)
    name_ext_dict = get_name_extensions('nan', params)
    saving_file_name = f'MSA_ID{params['ID']}_dispcritical_{name_ext_dict['mesh']}{name_ext_dict['mat']}.txt'
    saving_file_dir = os.path.join( saving_folder_dir, saving_file_name)
    np.savetxt( saving_file_dir, np.array([disp_critical]) )

def save_MSA_forcecritical( outputs_root_dir, params, force_critical):
    saving_folder_dir = os.path.join( outputs_root_dir, f'resultslcar{params['lcar']}' )
    os.makedirs( saving_folder_dir, exist_ok = True)
    name_ext_dict = get_name_extensions('nan', params)
    saving_file_name = f'MSA_ID{params['ID']}_forcecritical_{name_ext_dict['mesh']}{name_ext_dict['mat']}.txt'
    saving_file_dir = os.path.join( saving_folder_dir, saving_file_name)
    np.savetxt( saving_file_dir, np.array([force_critical]))

def save_MSA_eigenvector( outputs_root_dir, params, eigenvector):
    saving_folder_dir = os.path.join( outputs_root_dir, f'resultslcar{params['lcar']}' )
    os.makedirs( saving_folder_dir, exist_ok = True)
    name_ext_dict = get_name_extensions('nan', params)
    saving_file_name = f'MSA_ID{params['ID']}_eigenvector_{name_ext_dict['mesh']}{name_ext_dict['mat']}.txt'
    saving_file_dir = os.path.join( saving_folder_dir, saving_file_name)
    np.savetxt( saving_file_dir, eigenvector)

def save_MSA_results( outputs_root_dir, params, disp_critical, force_critical, eigenvector):
    save_MSA_dispcritical( outputs_root_dir, params, disp_critical)
    save_MSA_forcecritical( outputs_root_dir, params, force_critical)
    save_MSA_eigenvector( outputs_root_dir, params, eigenvector)

def save_FEA_dispcritical( outputs_root_dir, params, results_dict ):
    disp_critical = results_dict['disp_critical']
    saving_folder_dir = os.path.join( outputs_root_dir, f'resultslcar{params['lcar']}' )
    os.makedirs( saving_folder_dir, exist_ok = True)
    name_ext_dict = get_name_extensions('nan', params)
    saving_file_name = f'FEA_ID{params['ID']}_dispcritical_{name_ext_dict['save']}.txt'
    saving_file_dir = os.path.join( saving_folder_dir, saving_file_name)
    np.savetxt( saving_file_dir, np.array([disp_critical]) )

def save_FEA_forcecritical( outputs_root_dir, params, results_dict):
    force_critical = results_dict['force_critical']
    saving_folder_dir = os.path.join( outputs_root_dir, f'resultslcar{params['lcar']}' )
    os.makedirs( saving_folder_dir, exist_ok = True)
    name_ext_dict = get_name_extensions('nan', params)
    saving_file_name = f'FEA_ID{params['ID']}_forcecritical_{name_ext_dict['save']}.txt'
    saving_file_dir = os.path.join( saving_folder_dir, saving_file_name)
    np.savetxt( saving_file_dir, np.array([force_critical]))

def save_FEA_eigenvector( outputs_root_dir, params, results_dict):
    eigenvector = results_dict['eigmode_critical1']
    saving_folder_dir = os.path.join( outputs_root_dir, f'resultslcar{params['lcar']}' )
    os.makedirs( saving_folder_dir, exist_ok = True)
    name_ext_dict = get_name_extensions('nan', params)
    saving_file_name = f'FEA_ID{params['ID']}_eigenvector_{name_ext_dict['save']}.txt'
    saving_file_dir = os.path.join( saving_folder_dir, saving_file_name)
    np.savetxt( saving_file_dir, eigenvector)

def save_FEA_initconf( outputs_root_dir, params, results_dict):
    initconf = results_dict['initconf']
    saving_folder_dir = os.path.join( outputs_root_dir, f'resultslcar{params['lcar']}' )
    os.makedirs( saving_folder_dir, exist_ok = True)
    name_ext_dict = get_name_extensions('nan', params)
    saving_file_name = f'FEA_ID{params['ID']}_initconf_{name_ext_dict['save']}.txt'
    saving_file_dir = os.path.join( saving_folder_dir, saving_file_name)
    np.savetxt( saving_file_dir, initconf)

def save_FEA_results( outputs_root_dir, params, results_dict):
    save_FEA_dispcritical( outputs_root_dir, params, results_dict)
    save_FEA_forcecritical( outputs_root_dir, params, results_dict)
    save_FEA_eigenvector( outputs_root_dir, params, results_dict)
    save_FEA_initconf( outputs_root_dir, params, results_dict)
