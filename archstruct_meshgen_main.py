from archstruct_meshgen_functions import save_archstruct_mesh_fromstrkeylst
import gmsh
import pygmsh
import meshio
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import itertools
SEED = 5
np.random.seed(seed=SEED)

"""Information about the Structures"""
ny, nx = 4, 2
cell_size, lcar, z_prune, element_type = 20, 1, True, 'line'
cellvalues = [4, 5, 6] # cell types
cellvalues_str = ''.join(map(str, cellvalues))
if z_prune:
  dims = '2D'
else:
  dims = '3D'

"""Saving Directories"""
root_dir = 'inputs/archstruct'
mesh_folder_dir = os.path.join( root_dir, f'mesh/mesh{dims}/meshlcar{lcar}' )

keys_folder_name = 'keys'
keys_folder_dir = os.path.join( root_dir, keys_folder_name )

os.makedirs( mesh_folder_dir, exist_ok=True )
os.makedirs( keys_folder_dir, exist_ok=True )

"""Making the dataset"""
strkey_lst           = ['44444444', '44444455', '55664466', '66666666', '55445555',
                        '54444445', '55655655', '54466445', '46644664', '46555564']
save_archstruct_mesh_fromstrkeylst( strkey_lst, ny, nx, cell_size, cellvalues_str, lcar, mesh_folder_dir, element_type, z_prune, keys_folder_dir)