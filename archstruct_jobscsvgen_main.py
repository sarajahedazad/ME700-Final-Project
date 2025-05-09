import pickle
from pathlib import Path
from itertools import product
import os

import pandas as pd

# ─── Configuration Lists ─────────────────────────────────────────────────────
cell_size_lst           = [20] # Length of each cell
cellvalues_str_lst      = ['456'] # Types of cells
ny_lst                  = [4] # Number of cells in the y direction
nx_lst                  = [2] # Number of cells in the x direction
E_lst                   = [1000] # Young's modulus
nu_lst                  = [0.33] # Poisson's ratio
a_lst                = [0.25] # Radius of cross sectional area
stop_thresh_lst         = [0.2] # The threshhold in which the FEA analysis should be stopped (in case it is not finding the critical point)
rel_param_lst           = [0.5] # Parameter for solver
rel_param_flag_lst      = [True] # Flag to whether or not use the previous parameter
max_bs_iter_cnt_lst     = [120] # maximum number of iteratrions allowed in bisection method (in case of finding the critical point)
find_critical_load_lst  = [True] # whether or not to find the critical disp and load

# Parameter lists
# The numbers like '44444444' account for cell type arrangements
strkey_lst           = ['44444444', '44444455', '55664466', '66666666', '55445555',
                        '54444445', '55655655', '54466445', '46644664', '46555564']

abs_tol_lst             = [1e-10] # absolute tolerance for solver
rel_tol_lst             = [1e-10] # relative tolerance for solver
lcar_lst             = [1] # critical length for mesh
tstep_lst            = [0.001] # the displacement step that is used in FEA

pf_force_updown_lst  = [[0.0001,0.0001]] # Force perturbation factors: In force perturbation, two small forces are applied. One is placed uppper, the other lower. These are the values for them.
H_force_updown_lst   = [[45,35]] # the locations of applying force perturbation
pf_dispfield_lst     = [1] # Displacement field perturbation factor

# ─── Original parameter lists ────────────────────────────────────────────────
cell_values             = [4, 5, 6] 

base_dir                = Path('inputs/archstruct')
keys_dir                = base_dir / 'keys'
name_ext                = f'ny{ny_lst[0]}nx{nx_lst[0]}celltypes{cellvalues_str_lst[0]}'
jobparams_dir           = base_dir / 'jobparams_csv'
jobtype                 = 'samplecheck' # this goes in the name of the csv file

# ─── Load key dictionaries ─────────────────────────────────────────────────────
# example: ID 1 has the strkey of '44444444'
keyscsv_name = f'keys_{ name_ext }.csv'
keyscsv_dir = os.path.join( keys_dir, keyscsv_name )
keyscsv_df = pd.read_csv( keyscsv_dir )
ID2strkey_dict  = dict(zip(keyscsv_df['ID'], keyscsv_df['strkey'].astype(str)))
strkey2ID_dict  = dict(zip(keyscsv_df['strkey'].astype(str), keyscsv_df['ID']))


# ─── DataFrame Builder ───────────────────────────────────────────────────────
def get_parameters_dataframe(params: dict) -> pd.DataFrame:
    cols = [
        'strkey', 'ID', 'a', 'lcar', 'tstep',
        'pf_force_up', 'pf_force_down',
        'H_force_up', 'H_force_down', 'pf_dispfield',
        'cell_size', 'cellvalues_str', 'ny', 'nx',
        'E', 'nu', 'stop_thresh', 'rel_param', 'rel_param_flag',
        'abs_tol', 'rel_tol', 'max_bs_iter_cnt', 'find_critical_load'
    ]

    product_iter = product(
        params['strkey_lst'], params['a_lst'], params['lcar_lst'],
        params['tstep_lst'], params['pf_force_updown_lst'],
        params['H_force_updown_lst'], params['pf_dispfield_lst'],
        params['cell_size_lst'], params['cellvalues_str_lst'],
        params['ny_lst'], params['nx_lst'],
        params['E_lst'], params['nu_lst'], params['stop_thresh_lst'],
        params['rel_param_lst'], params['rel_param_flag_lst'],
        params['abs_tol_lst'], params['rel_tol_lst'], params['max_bs_iter_cnt_lst'],
        params['find_critical_load_lst']
    )

    records = []
    for (
        strkey, a, lcar, tstep,
        pf_updown, H_updown, pf_disp,
        cell_size, cellvalues_str, ny, nx,
        E, nu, stop_thresh, rel_param,
        rel_param_flag, abs_tol, rel_tol, max_bs_iter_cnt,
        find_crit
    ) in product_iter:
        records.append({
            'strkey':           strkey,
            'ID':               strkey2ID_dict[strkey],
            'a':                a,
            'lcar':             lcar,
            'tstep':            tstep,
            'pf_force_up':      pf_updown[0],
            'pf_force_down':    pf_updown[1],
            'H_force_up':       H_updown[0],
            'H_force_down':     H_updown[1],
            'pf_dispfield':     pf_disp,
            'cell_size':        cell_size,
            'cellvalues_str':   cellvalues_str,
            'ny':               ny,
            'nx':               nx,
            'E':                E,
            'nu':               nu,
            'stop_thresh':      stop_thresh,
            'rel_param':        rel_param,
            'rel_param_flag':   rel_param_flag,
            'abs_tol':          abs_tol,
            'rel_tol':          rel_tol,
            'max_bs_iter_cnt':  max_bs_iter_cnt,
            'find_critical_load': find_crit
        })

    return pd.DataFrame.from_records(records, columns=cols)

# ─── Helpers ─────────────────────────────────────────────────────────────────
def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def make_param_dict(
    strkey_lst, a_lst, lcar_lst, tstep_lst,
    pf_force_updown_lst, H_force_updown_lst, pf_dispfield_lst,
    cell_size_lst, cellvalues_str_lst, ny_lst, nx_lst,
    E_lst, nu_lst, stop_thresh_lst, rel_param_lst,
    rel_param_flag_lst, abs_tol_lst, rel_tol_lst, max_bs_iter_cnt_lst,
    find_critical_load_lst
) -> dict:
    return {
        'strkey_lst':            strkey_lst,
        'a_lst':                 a_lst,
        'lcar_lst':              lcar_lst,
        'tstep_lst':             tstep_lst,
        'pf_force_updown_lst':   pf_force_updown_lst,
        'H_force_updown_lst':    H_force_updown_lst,
        'pf_dispfield_lst':      pf_dispfield_lst,
        'cell_size_lst':         cell_size_lst,
        'cellvalues_str_lst':    cellvalues_str_lst,
        'ny_lst':                ny_lst,
        'nx_lst':                nx_lst,
        'E_lst':                 E_lst,
        'nu_lst':                nu_lst,
        'stop_thresh_lst':       stop_thresh_lst,
        'rel_param_lst':         rel_param_lst,
        'rel_param_flag_lst':    rel_param_flag_lst,
        'abs_tol_lst':           abs_tol_lst,
        'rel_tol_lst':           rel_tol_lst,
        'max_bs_iter_cnt_lst':   max_bs_iter_cnt_lst,
        'find_critical_load_lst': find_critical_load_lst
    }

# ─── CSV Saver ────────────────────────────────────────────────────────────────
def save_jobs_csv(
    save_dir: Path, perturb_type: str, jobtype: str,
    strkey_lst, a_lst, lcar_lst, tstep_lst,
    pf_force_updown_lst, H_force_updown_lst, pf_dispfield_lst,
    cell_size_lst, cellvalues_str_lst, ny_lst, nx_lst,
    E_lst, nu_lst, stop_thresh_lst, rel_param_lst,
    rel_param_flag_lst, abs_tol_lst, rel_tol_lst, max_bs_iter_cnt_lst,
    find_critical_load_lst
):
    if perturb_type == 'nan':
        pf_force_updown_lst = [[0,0]]
        H_force_updown_lst  = [['nan','nan']]
        pf_dispfield_lst    = [0]
    elif perturb_type == 'smallforce':
        pf_dispfield_lst    = [0]
    elif perturb_type == 'dispfield':
        pf_force_updown_lst = [[0,0]]
        H_force_updown_lst  = [['nan','nan']]

    params = make_param_dict(
        strkey_lst, a_lst, lcar_lst, tstep_lst,
        pf_force_updown_lst, H_force_updown_lst,
        pf_dispfield_lst,
        cell_size_lst, cellvalues_str_lst,
        ny_lst, nx_lst,
        E_lst, nu_lst,
        stop_thresh_lst, rel_param_lst,
        rel_param_flag_lst,
        abs_tol_lst, rel_tol_lst,
        max_bs_iter_cnt_lst,
        find_critical_load_lst
    )

    df = get_parameters_dataframe(params)
    out_file = save_dir / f'jobs_{jobtype}_perturb{perturb_type}.csv'
    df.to_csv(out_file, index=False)
    print(f"→ Saved {len(df)} jobs to {out_file}")

# ─── Main ────────────────────────────────────────────────────────────────────

ensure_dir(jobparams_dir)

for perturb in ['nan','smallforce','dispfield']:
    save_jobs_csv(
        jobparams_dir, perturb, jobtype,
        strkey_lst, a_lst, lcar_lst, tstep_lst,
        pf_force_updown_lst, H_force_updown_lst, pf_dispfield_lst,
        cell_size_lst, cellvalues_str_lst,
        ny_lst, nx_lst,
        E_lst, nu_lst,
        stop_thresh_lst, rel_param_lst,
        rel_param_flag_lst,
        abs_tol_lst, rel_tol_lst,
        max_bs_iter_cnt_lst,
        find_critical_load_lst
    )
