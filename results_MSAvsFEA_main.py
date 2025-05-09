import numpy as np
import matplotlib.pyplot as plt
from reading_saving_functions import *
import os
from MSA_demo import get_frame_from_params
from MSA_shapefunctions_2D import ShapeFunctions2D

# Define directories for saving results
MSA_parent_output_folder_dir = 'outputs/archstruct/results2D/MSA'
FEA_parent_output_folder_dir = 'outputs/archstruct/results2D/FEA'
mesh_folder_root = 'inputs/archstruct/mesh/mesh2D'
wrappedup_folder = 'outputs/archstruct/results2D/wrappedup'

# Create the folder for wrapped-up results if it doesn't exist
os.makedirs(wrappedup_folder, exist_ok=True)

# Directory for input parameters
parent_dir_inputs = 'inputs/archstruct'
jobparams_csv_folder_dir = os.path.join(parent_dir_inputs, 'jobparams_csv')

# Define job type and perturbation type for the analysis
jobtype = 'samplecheck'
perturb_type = 'nan'

# Read the job parameter CSV file
csv_file = os.path.join(jobparams_csv_folder_dir, f"jobs_{jobtype}_perturb{perturb_type}.csv")
df = pd.read_csv(csv_file)

# Initialize lists for storing critical displacement and force values for MSA and FEA
job_id_lst = list(range(len(df)))
dispcritical_MSA_lst, dispcritical_FEA_lst = [], []
forcecritical_MSA_lst, forcecritical_FEA_lst = [], []
dispcritical_percentageerror_lst = []
forcecritical_percentageerror_lst = []
ID_lst = []
strkey_lst = []

# Loop through each job ID in the parameter file
for job_id in job_id_lst:
    # Extract parameters for the current job
    params = get_params_dict_from_df(job_id, df)
    name_ext_dict = get_name_extensions(perturb_type, params)
    lcar = params['lcar']

    # Define output directories for MSA and FEA results
    MSA_output_folder_dir = os.path.join(MSA_parent_output_folder_dir, f'resultslcar{lcar}')
    FEA_output_folder_dir = os.path.join(FEA_parent_output_folder_dir, f'resultslcar{lcar}')

    # Load critical displacement and force values for MSA
    dispcritical_MSA_name = f'MSA_ID{params["ID"]}_dispcritical_{name_ext_dict["mesh"]}{name_ext_dict["mat"]}.txt'
    forcecritical_MSA_name = f'MSA_ID{params["ID"]}_forcecritical_{name_ext_dict["mesh"]}{name_ext_dict["mat"]}.txt'
    dispcritical_MSA = np.loadtxt(os.path.join(MSA_output_folder_dir, dispcritical_MSA_name))
    forcecritical_MSA = np.loadtxt(os.path.join(MSA_output_folder_dir, forcecritical_MSA_name))

    # Load critical displacement and force values for FEA
    dispcritical_FEA_name = f'FEA_ID{params["ID"]}_dispcritical_{name_ext_dict["save"]}.txt'
    forcecritical_FEA_name = f'FEA_ID{params["ID"]}_forcecritical_{name_ext_dict["save"]}.txt'
    dispcritical_FEA = np.loadtxt(os.path.join(FEA_output_folder_dir, dispcritical_FEA_name))
    forcecritical_FEA = np.loadtxt(os.path.join(FEA_output_folder_dir, forcecritical_FEA_name))

    # Append values to the respective lists
    dispcritical_MSA_lst.append(dispcritical_MSA)
    forcecritical_MSA_lst.append(forcecritical_MSA)
    dispcritical_FEA_lst.append(dispcritical_FEA)
    forcecritical_FEA_lst.append(forcecritical_FEA)

    # Calculate percentage errors for displacement and force
    dispcritical_percentageerror = 100 * np.abs(dispcritical_MSA - dispcritical_FEA) / dispcritical_MSA
    dispcritical_percentageerror_lst.append(dispcritical_percentageerror)

    forcecritical_percentageerror = 100 * np.abs(forcecritical_MSA - forcecritical_FEA) / forcecritical_MSA
    forcecritical_percentageerror_lst.append(forcecritical_percentageerror)

    # Append job ID and string key for tracking
    ID_lst.append(params['ID'])
    strkey_lst.append(params['strkey'])

# Create a dictionary for saving results in a DataFrame
msavsfea_dict = {
    'ID': ID_lst, 
    'strkey': strkey_lst,
    'dispcritical_MSA': dispcritical_MSA_lst, 
    'forcecritical_MSA': forcecritical_MSA_lst,
    'dispcritical_FEA': dispcritical_FEA_lst, 
    'forcecritical_FEA': forcecritical_FEA_lst,
    'dispcritical_percentageerror': dispcritical_percentageerror_lst, 
    'forcecritical_percentageerror': forcecritical_percentageerror_lst
}

# Save results to a CSV file
df_results = pd.DataFrame.from_dict(msavsfea_dict)
df_results.to_csv(os.path.join(wrappedup_folder, f'results_{jobtype}.csv'))

# Plotting critical displacement values for MSA and FEA
plt.figure()
plt.plot(list(range(len(ID_lst))), dispcritical_MSA_lst, marker='.', label='MSA')
plt.plot(list(range(len(ID_lst))), dispcritical_FEA_lst, marker='.', label='FEA')
plt.xticks(list(range(len(ID_lst))), ID_lst)
plt.xlabel('Sample ID')
plt.ylabel('Critical Displacement')
plt.legend()
plt.savefig(os.path.join(wrappedup_folder, f'jobs{jobtype}_dispcritical.png'), dpi=500)
plt.close()

# Plotting critical force values for MSA and FEA
plt.figure()
plt.plot(list(range(len(ID_lst))), forcecritical_MSA_lst, marker='.', label='MSA')
plt.plot(list(range(len(ID_lst))), forcecritical_FEA_lst, marker='.', label='FEA')
plt.xticks(list(range(len(ID_lst))), ID_lst)
plt.xlabel('Sample ID')
plt.ylabel('Critical Force')
plt.legend()
plt.savefig(os.path.join(wrappedup_folder, f'jobs{jobtype}_forcecritical.png'), dpi=500)
plt.close()

# Plotting percentage errors for displacement and force
plt.figure()
plt.plot(list(range(len(ID_lst))), dispcritical_percentageerror_lst, marker='.', label='Critical Displacement Error')
plt.plot(list(range(len(ID_lst))), forcecritical_percentageerror_lst, marker='.', label='Critical Force Error')
plt.xticks(list(range(len(ID_lst))), ID_lst)
plt.xlabel('Sample ID')
plt.ylabel('Percentage Error ( 100 * abs(FEA - MSA)/ MSA )')
plt.legend()
plt.savefig(os.path.join(wrappedup_folder, f'jobs{jobtype}_errors.png'), dpi=500)
plt.close()

# Comparison of eigen modes for each job
for job_id in job_id_lst:
    # Load eigenvectors for MSA and FEA
    params = get_params_dict_from_df(job_id, df)
    name_ext_dict = get_name_extensions(perturb_type, params)
    
    eigenvector_MSA_name = f'MSA_ID{params["ID"]}_eigenvector_{name_ext_dict["mesh"]}{name_ext_dict["mat"]}.txt'
    eigenvector_MSA = np.loadtxt(os.path.join(MSA_output_folder_dir, eigenvector_MSA_name))

    eigenvector_FEA_name = f'FEA_ID{params["ID"]}_eigenvector_{name_ext_dict["save"]}.txt'
    eigenvector_FEA = np.loadtxt(os.path.join(FEA_output_folder_dir, eigenvector_FEA_name))
    
    # Load initial configuration for FEA
    initconf_FEA_name = f'FEA_ID{params["ID"]}_initconf_{name_ext_dict["save"]}.txt'
    initconf_FEA = np.loadtxt(os.path.join(FEA_output_folder_dir, initconf_FEA_name))

    # Generate and save MSA eigenmode visualization
    frame = get_frame_from_params(mesh_folder_root, params)
    saving_dir_with_name = os.path.join(wrappedup_folder, f'MSA_ID{params["ID"]}_eigenmode.png')
    sf = ShapeFunctions2D(eigenvector_MSA, frame, n=2, scale=100)
    sf.plot_element_interpolation(saving_dir_with_name, f'MSA Eigenmode, ID = {params["ID"]}')

    # Plot and save FEA eigenmode visualization
    conf = initconf_FEA + 140 * eigenvector_FEA
    plt.figure()
    plt.scatter(initconf_FEA[:, 0], initconf_FEA[:, 1], s=1, c='k', label='Init Conf')
    plt.scatter(conf[:, 0], conf[:, 1], s=1, c='r', label='Final Conf')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.axis('equal')
    plt.title(f'FEA Eigenmode, ID = {params["ID"]}')
    plt.savefig(os.path.join(wrappedup_folder, f'FEA_ID{params["ID"]}_eigenmode.png'), dpi=500)
    plt.close()

    # Print progress
    print(f'Jobs done: {job_id + 1} / {len(job_id_lst)}')
