import numpy as np
from MSA_geometry_2D import *
from MSA_boundaryconditions_2D import *
from MSA_stiffnessmatrices_2D import *
from MSA_shapefunctions_2D import *
from MSA_solver import *
import os
from dolfin import * 
import matplotlib.pyplot as plt
from reading_saving_functions import *


def read_mesh( mesh_file_dir ):
    mesh = Mesh()
    with XDMFFile(mesh_file_dir) as infile:
        infile.read(mesh)
    return mesh

def frame_from_mesh( mesh, E, A, I ):
    frame = Frame2D()
    frame.generate_frame_from_mesh( mesh, E, A, I )
    return frame

def get_frame_from_params( mesh_folder_root, params_dict):
    lcar = params_dict['lcar']
    name_ext_dict = get_name_extensions('nan', params_dict)
    mesh_folder_dir = os.path.join( mesh_folder_root, f'meshlcar{lcar}' )
    mesh_file_name = f'mesh{name_ext_dict['mesh']}ID{ params_dict['ID'] }.xdmf'
    mesh_file_dir = os.path.join( mesh_folder_dir, mesh_file_name)
    E = params_dict['E']
    a = params_dict['a']
    A = np.pi * a **2
    I = np.pi * a **4/4
    mesh = read_mesh( mesh_file_dir )
    frame = frame_from_mesh( mesh, E, A, I )
    return frame


def apply_fixed_end_and_displaced_end_bcs( frame, structure_length, applied_disp, direction = 'y'):
    bcs = BoundaryConditions2D(frame)
    bcs.apply_boundary_conditions_by_direction(direction, 0, (0, 0), (0), tol=1e-8)
    bcs.apply_boundary_conditions_by_direction(direction, structure_length, (0, -applied_disp), (0), tol=1e-8)
    bcs.apply_free_conditions()
    bcs.set_up_bounds()
    return bcs

def calc_dispcritical_eigmode_fixed_to_displaced_end( frame, structure_length, arbitrary_disp = 0.01, direction= 'y'):
    '''-------------------------------'''
    '''-----------BCs-----------------'''
    '''-------------------------------'''
    bcs = apply_fixed_end_and_displaced_end_bcs( frame, structure_length, arbitrary_disp, direction)
    '''-------------------------------------------------------------'''
    '''-----------Build the Global Stiffness Matrix-----------------'''
    '''-------------------------------------------------------------'''
    stiffmat = StiffnessMatrices2D(frame)
    K = stiffmat.get_global_elastic_stiffmatrix()
    '''-------------------------------------------------------------'''
    '''-------------------Solving for unknowns----------------------'''
    '''-------------------------------------------------------------'''
    Delta, F = solve_stiffness_system(K, bcs)
    '''-------------------------------------------------------------'''
    '''----------------Critical Displacements-------------'''
    '''-------------------------------------------------------------'''
    K_g = stiffmat.get_global_geometric_stiffmatrix(Delta)

    smallest_eigvalue, eigenvectors_allstructure = compute_critical_load(K, K_g, bcs)

    disp_critical = arbitrary_disp * smallest_eigvalue

    return disp_critical, eigenvectors_allstructure

def calc_forcecritical_fixed_to_displaced_end( frame, structure_length, disp_critical, direction= 'y' ):
    '''-------------------------------'''
    '''-----------BCs-----------------'''
    '''-------------------------------'''
    bcs = apply_fixed_end_and_displaced_end_bcs( frame, structure_length, disp_critical, direction)
    '''-------------------------------------------------------------'''
    '''-----------Build the Global Stiffness Matrix-----------------'''
    '''-------------------------------------------------------------'''
    stiffmat = StiffnessMatrices2D(frame)
    K = stiffmat.get_global_elastic_stiffmatrix()
    '''-------------------------------------------------------------'''
    '''-------------------Solving for unknowns----------------------'''
    '''-------------------------------------------------------------'''
    Delta, F = solve_stiffness_system(K, bcs)
    '''-------------------------------------------------------------'''
    '''----------------Critical Forces-------------'''
    '''-------------------------------------------------------------'''
    K_g = stiffmat.get_global_geometric_stiffmatrix(Delta)

    f_bound = bcs.get_force_on_boundary( direction, 0, F)
    force_critical = np.sum( f_bound )

    return force_critical