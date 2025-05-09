import gmsh
import pygmsh
import meshio
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import pickle
import itertools
SEED = 5
np.random.seed(seed=SEED)

def save_mesh(file_path, mesh_object, element_type, z_prune):
  # Element Type
  element_dict = {}
  for element in mesh_object.cells:
    if element.type == element_type:
      element_dict[ element_type ] = element.data

  # Cell Points
  if z_prune:
    mesh_points = mesh_object.points[ :, :2]
  else:
    mesh_points = mesh_object.points

  # Save
  meshio.write( file_path + ".xdmf", meshio.Mesh(points = mesh_points, cells = element_dict ))

# Each cell type, can be represented with a number
def get_cellkeys_func_dict():
    cellkeys_function_dict = { 0: cell_empty, 
                               4: cell_horizental, 5: cell_vertical, 6: cell_cross}

    return  cellkeys_function_dict

def cell_empty( geom,x_coord, y_coord, cell_size, lcar): # Right to left, diagonal
    p1 = geom.add_point( [x_coord, y_coord, 0], lcar)
    p2 = geom.add_point( [x_coord + cell_size, y_coord, 0], lcar)
    p3 = geom.add_point( [x_coord + cell_size, y_coord + cell_size, 0], lcar)
    p4 = geom.add_point( [x_coord, y_coord + cell_size, 0], lcar)

    l1 = geom.add_line( p1, p2 )
    l2 = geom.add_line( p2, p3 )
    l3 = geom.add_line( p3, p4 )
    l4 = geom.add_line( p1, p4 )

    return [l1, l2, l3, l4]


def cell_horizental(geom,x_coord, y_coord, cell_size, lcar): # Left ro right, diagonal
    p1 = geom.add_point( [x_coord, y_coord, 0], lcar)
    p2 = geom.add_point( [x_coord + cell_size, y_coord, 0], lcar)
    p3 = geom.add_point( [x_coord + cell_size, y_coord + cell_size, 0], lcar)
    p4 = geom.add_point( [x_coord, y_coord + cell_size, 0], lcar)

    p5 = geom.add_point( [x_coord, y_coord + cell_size / 2, 0], lcar)
    p6 = geom.add_point( [x_coord + cell_size, y_coord + cell_size / 2, 0], lcar)

    l1 = geom.add_line( p1, p2 )
    l2 = geom.add_line( p2, p3 )
    l3 = geom.add_line( p3, p4 )
    l4 = geom.add_line( p4, p1 )
    l5 = geom.add_line( p5, p6 )

    return [l1, l2, l3, l4, l5]

def cell_vertical(geom,x_coord, y_coord, cell_size, lcar): # Left ro right, diagonal
    p1 = geom.add_point( [x_coord, y_coord, 0], lcar)
    p2 = geom.add_point( [x_coord + cell_size, y_coord, 0], lcar)
    p3 = geom.add_point( [x_coord + cell_size, y_coord + cell_size, 0], lcar)
    p4 = geom.add_point( [x_coord, y_coord + cell_size, 0], lcar)

    p5 = geom.add_point( [x_coord + cell_size / 2, y_coord, 0], lcar)
    p6 = geom.add_point( [x_coord + cell_size / 2, y_coord + cell_size, 0], lcar)

    l1 = geom.add_line( p1, p2 )
    l2 = geom.add_line( p2, p3 )
    l3 = geom.add_line( p3, p4 )
    l4 = geom.add_line( p4, p1 )
    l5 = geom.add_line( p5, p6 )

    return [l1, l2, l3, l4, l5]

def cell_cross(geom,x_coord, y_coord, cell_size, lcar): # Left ro right, diagonal
    p1 = geom.add_point( [x_coord, y_coord, 0], lcar)
    p2 = geom.add_point( [x_coord + cell_size, y_coord, 0], lcar)
    p3 = geom.add_point( [x_coord + cell_size, y_coord + cell_size, 0], lcar)
    p4 = geom.add_point( [x_coord, y_coord + cell_size, 0], lcar)

    p5 = geom.add_point( [x_coord + cell_size / 2, y_coord, 0], lcar)
    p6 = geom.add_point( [x_coord + cell_size / 2, y_coord + cell_size, 0], lcar)

    p7 = geom.add_point( [x_coord, y_coord + cell_size / 2, 0], lcar)
    p8 = geom.add_point( [x_coord + cell_size, y_coord + cell_size / 2, 0], lcar)

    p9 = geom.add_point( [x_coord + cell_size / 2, y_coord + cell_size / 2, 0], lcar)

    l1 = geom.add_line( p1, p5 )
    l2 = geom.add_line( p5, p2 )
    l3 = geom.add_line( p2, p8 )
    l4 = geom.add_line( p8, p3 )

    l5 = geom.add_line( p3, p6 )
    l6 = geom.add_line( p6, p4 )
    l7 = geom.add_line( p4, p7 )
    l8 = geom.add_line( p7, p1 )

    l9 = geom.add_line( p5, p9 )
    l10 = geom.add_line( p6, p9 )
    l11 = geom.add_line( p7, p9 )
    l12 = geom.add_line( p8, p9 )

    return [l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12]

def generate_archstruct_mesh( cellkeys_arr, ny, nx, cell_size, lcar ):
    cellkeys_func_dict  = get_cellkeys_func_dict()
    lst = []
    with pygmsh.occ.Geometry() as geom:
        for j in range( ny ):
            for i in range( nx ):
                x_coord = i * cell_size
                y_coord = j * cell_size
                cell_key = cellkeys_arr[j][i]
                lst_cell = cellkeys_func_dict[ cell_key ]( geom, x_coord, y_coord, cell_size, lcar )
                lst = lst + lst_cell
        archstruct = geom.boolean_union( lst )
        gmsh.option.setNumber("General.ExpertMode", 1)
        archstruct_mesh = geom.generate_mesh(  )
    return archstruct_mesh

def cellkeys_array2string(arr):
    return ''.join(str(x) for row in arr for x in row)

def cellkeys_string2array(strkey, ny, nx):
    if len(strkey) != ny * nx:
        raise ValueError(f"String length ({len(strkey)}) does not match requested shape {(ny, nx)}")
    flat = [int(ch) for ch in strkey]
    return np.array(flat, dtype=int).reshape((ny, nx))
  
def save_archstruct_mesh( cellkeys_arr, ny, nx, cell_size, lcar, mesh_file_name, saving_folder_dir, element_type, z_prune ):
    lattice_mesh = generate_archstruct_mesh( cellkeys_arr, ny, nx, cell_size, lcar )
    mesh_file_path = os.path.join( saving_folder_dir, mesh_file_name )
    save_mesh( mesh_file_path, lattice_mesh, element_type, z_prune )

def save_archstruct_mesh_fromstrkeylst(strkey_lst, ny, nx, cell_size, cellvalues_str, lcar, mesh_folder_dir, element_type, z_prune, keys_folder_dir):
    name_extention1 = f'ny{ny}nx{nx}celltypes{cellvalues_str}'
    name_extention2 = f'cellsize{cell_size}lcar{lcar}'
    name_extention = f'{name_extention1}{name_extention2}'
    ID_lst = []
    strkeys_lst = []
    for k, strkey in enumerate( strkey_lst ):
      cellkeys_arr = cellkeys_string2array( strkey, ny, nx)
      mesh_file_name = f'mesh{name_extention}ID{ k + 1 }'
      save_archstruct_mesh( cellkeys_arr, ny, nx, cell_size, lcar, mesh_file_name, mesh_folder_dir, element_type, z_prune )
      ID_lst.append( k + 1 )
      strkeys_lst.append( cellkeys_array2string( cellkeys_arr ) )
    ID_strkey = {'ID': ID_lst, 'strkey': strkeys_lst}
    df_ID_strkey = pd.DataFrame.from_dict(ID_strkey)
    keys_file_name = f'keys_{name_extention1}.csv'
    keys_file_dir = os.path.join( keys_folder_dir, keys_file_name)
    df_ID_strkey.to_csv( keys_file_dir )
