from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import h5py
import pickle
import time
import pandas as pd
import math

# Dealing with ufl legacy
try:
    from ufl import diag, Jacobian, shape
except ImportError:
    from ufl_legacy import diag, Jacobian, shape


# ---------------------- FEA Utility Functions ----------------------
def read_mesh(folder_dir, file_name):
    mesh = Mesh()
    with XDMFFile(os.path.join(folder_dir, file_name)) as infile:
        infile.read(mesh)
    return mesh


def get_functionspace( mesh ):
    Ue = VectorElement("CG", mesh.ufl_cell(), 2, dim=2) # displacement
    Te = FiniteElement("CG", mesh.ufl_cell(), 1) # rotation
    V = FunctionSpace(mesh, MixedElement([Ue, Te]))   

    v_ = TestFunction(V)
    u_, theta_ = split(v_)
    dv = TrialFunction(V)
    v = Function(V, name="Generalized displacement")
    u, theta = split(v)
    return V, v_, u_, theta_, dv, v, u, theta


def compute_tangents(mesh):
    Jac  = Jacobian(mesh)
    gdim = mesh.geometry().dim()
    Jac  = as_vector([Jac[i,0] for i in range(gdim)])
    g01  = Jac/sqrt(dot(Jac,Jac))
    g02  = as_vector([-g01[1], g01[0]])
    R0   = outer(g01, as_vector([1,0])) + outer(g02, as_vector([0,1]))
    return g01, R0


def tgrad(u, g01):
    return dot(grad(u), g01)


def rotation_matrix(theta):
    return as_tensor([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])


def strains(R0, Rot, g01, u, theta):
    defo = dot(R0.T, dot(Rot.T, g01 + tgrad(u, g01)) - g01)
    curv = tgrad(theta, g01)
    return defo, curv


def constitutive_matrices(a, E, nu):
    S = Constant(pi * a**2)
    I = Constant(pi * a**4 / 4)
    kappa = Constant(6/7)
    G = E/(2*(1+nu))
    ES   = E*S
    EI   = E*I
    GS2  = G*S*kappa
    C_N  = diag(as_vector([ES, GS2]))
    return C_N, EI


def eigvalvect(tangent_form, bcs, N_solve, V):
    A = PETScMatrix()
    assemble(tangent_form, tensor=A)
    for bc in bcs:
        bc.apply(A)
    solver = SLEPcEigenSolver(A)
    solver.parameters['spectral_transform'] = 'shift-and-invert'
    solver.parameters['spectral_shift']     = 0.0
    solver.parameters['tolerance']          = 1e-16
    solver.solve(N_solve)

    vals, vecs = [], []
    for i in range(N_solve):
        r, _, rx, _ = solver.get_eigenpair(i)
        vals.append(r)
        mode = Function(V)
        mode.vector()[:] = rx
        vecs.append(mode)
    return vals, vecs


def not_same_sign(a, b):
    return a * b < 0

# ---------------------- Name Extension ----------------------
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

# ---------------------- Core Solver ----------------------
def run_fea(perturb_type, params_dict):
    parameters["form_compiler"]["cpp_optimize"] = True
    parameters["form_compiler"]["quadrature_degree"] = 3
    parameters['reorder_dofs_serial'] = False

    ffc_options = {"optimize": True, \
                "eliminate_zeros": True, \
                "precompute_basis_const": True, \
                "precompute_ip_const": True}
    #-----------------Inputs--------------------------
    ''' Variable Parameters'''
    strkey = params_dict['strkey'] # A string of numbers that shows that corresponds 
    ID = int( params_dict['ID'] ) # A string of numbers that shows that corresponds
    a = params_dict['a'] # radius of the cross-sectional area
    lcar = params_dict['lcar'] # characteristic length of the mesh
    tstep = params_dict['tstep']
    # two forces are applied as perturbation. One is applied on a higher height, the other to a lower one.
    pf_force_up = params_dict['pf_force_up'] # Perturbation factor: force applied to location1, that is placed higher
    pf_force_down = params_dict['pf_force_down'] # Perturbation factor: force applied to location2, that is placed lower
    H_force_up = params_dict['H_force_up'] # Height: force applied to location1, that is placed higher
    H_force_down = params_dict['H_force_down'] # Height: force applied to location2, that is placed lower
    pf_dispfield = params_dict['pf_dispfield' ] # Perturbation factor in case the perturbation is applied to the displacement field

    if perturb_type != 'smallforce':
        H1, H2 = 45, 35
    else:
        H1, H2 = H_force_up, H_force_down

    ''' Fixed Parameters'''
    cell_size = params_dict['cell_size'] # Size of each side of the square-shaped cells that are the building blocks of the structure
    cellvalues_str = params_dict['cellvalues_str']
    ny, nx = params_dict['ny'], params_dict['nx'] # ny, nx: number of cells on the horizontal and vertical axes respectively
    E = params_dict['E']
    nu = params_dict['nu']

    ''' Solver Settings'''
    stop_thresh = params_dict['stop_thresh'] # the threshhold at which the simulation should stop
    rel_param = params_dict['rel_param']
    abs_tol = params_dict['abs_tol']
    rel_tol = params_dict['rel_tol']
    max_bs_iter_cnt = params_dict['max_bs_iter_cnt']

    iters = int( stop_thresh / tstep )
    t_lst = [tstep*(i) for i in range(iters + 1 )]
    #----------------Directories----------------------
    # Reading

    mesh_file_name   = f"mesh{get_name_extensions(perturb_type, params_dict)['mesh']}ID{params_dict['ID']}.xdmf"

    parent_dir_inputs = 'inputs/archstruct'
    mesh_folder_dir = os.path.join( parent_dir_inputs, f'mesh/mesh2D/meshlcar{lcar}')

    # Saving

    #----------------Reading Mesh and Setting Spaces and Domains-----------------
    mesh = read_mesh( mesh_folder_dir, mesh_file_name )
    V, v_, u_, theta_, dv, v, u, theta = get_functionspace( mesh )

    def load_center_up(x, on_boundary):
        return near(x[1], H1 , 1e-6) and near(x[0], cell_size * nx / 2 , 1e-6)
    def load_center_down(x, on_boundary):
        return near(x[1], H2 , 1e-6) and near(x[0], cell_size * nx / 2 , 1e-6)


    facets = MeshFunction("size_t", mesh, 0)
    facets.set_all(0)
    AutoSubDomain(load_center_up).mark(facets,1)
    AutoSubDomain(load_center_down).mark(facets,2)

    #------------------- Boundary Conditions
    def bottom(x, on_boundary):
        return near(x[1], 0, 1e-6) 
    def top(x, on_boundary):
        return near(x[1], ny * cell_size, 1e-6) 

    BC_bot = DirichletBC(V, Constant((0.0,0.0,0.0)), bottom, method='pointwise') # fixed displacement and rotation
    BC_top_x = DirichletBC(V.sub(0).sub(0), Constant(0.0), top, method='pointwise') # fix displacement
    BC_top_rot = DirichletBC(V.sub(1), Constant(0.0), top, method='pointwise') # fix rotation

    apply_disp = Expression("-t", t=0.0, degree = 0) # Create expression to compress the top
    BC_top_y = DirichletBC(V.sub(0).sub(1),apply_disp,top, 'pointwise') # incrementally compress the top 

    bcs = [BC_bot, BC_top_y, BC_top_rot, BC_top_x]

    #--------Kinematic and Weak Forms-----------
    g01, R0 = compute_tangents( mesh )
    Rot = rotation_matrix(theta)
    defo, curv = strains(R0, Rot, g01, u, theta)
    C_N, EI = constitutive_matrices( a, E , nu )

    F_max_up = Constant((pf_force_up, 0.0))
    F_max_down = Constant((pf_force_down, 0.0))

    dS = Measure("dS", domain=mesh, subdomain_data=facets)
    dx = Measure("dx", domain=mesh)

    elastic_energy = 0.5 * (dot(defo, dot(C_N, defo)) + (EI*curv**2))*dx

    F_int = derivative(elastic_energy, v, v_)
    F_ext = avg( dot(F_max_up, u_) ) * dS(1) + avg( dot(F_max_down, u_) ) * dS(2)

    residual = F_int - F_ext
    tangent_form = derivative(residual, v, dv)


    # boundary conditions & solver
#-------Solving----------
    problem = NonlinearVariationalProblem(residual, v, bcs, tangent_form)
    solver  = NonlinearVariationalSolver(problem)

    prm = solver.parameters

    prm['newton_solver']['linear_solver'] = "default"
    prm['newton_solver']['absolute_tolerance'] = abs_tol
    prm['newton_solver']['relative_tolerance'] = rel_tol
    if params_dict['rel_param_flag']:
        prm['newton_solver']['relaxation_parameter'] = params_dict['rel_param']

    N_eig_solve = 3
    eigval_lst_dict = {}
    for i in range( N_eig_solve ):
        eigval_lst_dict[ f'eigval{i + 1}' ] = []
    eigmode_critical_dict = {}
    for i in range( N_eig_solve ):
        eigmode_critical_dict[ f'eigmode_critical{i + 1}' ] = None
    v_reac = Function(V)
    bcRy = DirichletBC(V.sub(0).sub(1), Constant(1.0), bottom)
    force_lst = []
    disp_lst = []
    dispfield_lst = []

    #-----------------------Iterations-----------------------------------#
    force_critical, disp_critical, eigval_critical= 1000, 1000, 1000 # saves a big value for these in case those cannot be calculated
    eigmode_critical_vect = None

    findcriticalstate_flag = params_dict['find_critical_load']
    start_time = time.time()
    for iter, t_n in enumerate( t_lst ):
        # Performing iterations
        print( f'This is iteration number {iter + 1}')
        print(f'Horizontal dispalcement is {t_n}')
        apply_disp.t = t_n
        solver.solve()
        print( f'Iteration number {iter + 1} is done')

        # Appending data
        dispfield_lst.append(v.vector()[:])
        disp_lst.append( t_n )
        bcRy.apply(v_reac.vector())
        force = assemble(action(residual, v_reac))
        force_lst.append( force )

        # Analyzing eigenvalues
        eigenvalues_1toN_lst, eigenvectors_1toN_lst = eigvalvect( tangent_form, bcs, N_eig_solve, V )


        for i in range( N_eig_solve ):
            eigval_lst_dict[ f'eigval{i + 1}' ].append( eigenvalues_1toN_lst[ i ] )
        if iter > 0 and perturb_type != 'smallforce' and findcriticalstate_flag:
            a = disp_lst[ -2 ] # lower bound in the bisection method
            b = disp_lst[ -1 ] # upper bound in the bisection method
            f_a = eigval_lst_dict[ 'eigval1' ][ -2 ] # function at a (in this case, first eigenvalue)
            f_b = eigval_lst_dict[ 'eigval1' ][ -1 ] # function at b 
            if not_same_sign(f_a, f_b):
                c = (b + a ) / 2 
                apply_disp.t = c
                solver.solve()
                eigenvalues_1toN_lst, eigenvectors_1toN_lst = eigvalvect( tangent_form, bcs, N_eig_solve, V )
                f_c = eigenvalues_1toN_lst[ 0 ]
                bs_iter_cnt = 1 # bisection iteration count
                # Performing bisection method to find the critical point
                while ( np.abs(b - a ) >  100 * rel_tol and np.abs( f_c ) > 100 * abs_tol ) and max_bs_iter_cnt > bs_iter_cnt: 
                    if  a * c > 0:
                        f_a = f_c
                        a = c
                    else:
                        f_b = f_c
                        b = c
                    c = (b + a ) / 2
                    apply_disp.t = c
                    solver.solve()
                    eigenvalues_1toN_lst, eigenvectors_1toN_lst = eigvalvect( tangent_form, bcs, N_eig_solve, V )
                    f_c = eigenvalues_1toN_lst[ 0 ]
                    
                    print( f'Eigenvalue is {f_c}, for t { c }' )
                    bs_iter_cnt += 1 

                findcriticalstate_flag = False

                disp_critical = a
                apply_disp.t = disp_critical
                solver.solve()
                bcRy.apply(v_reac.vector())
                force_critical = assemble(action(residual, v_reac))
                eigenvalues_1toN_lst, eigenvectors_1toN_lst = eigvalvect( tangent_form, bcs, N_eig_solve, V )
                eigval_critical = eigenvalues_1toN_lst[0]
                eigmode_critical = eigenvectors_1toN_lst[0]

                eigmode_critical_dict = {}
                for i in range( N_eig_solve ):
                    eigmode_critical_dict[ f'eigmode_critical{i + 1}' ] =  eigenvectors_1toN_lst[ i ].vector()[:]

                # adding the eigenmode at the critical point as a perturbation to the last stable state
                if perturb_type == 'dispfield':
                    #last_stable_step = eigval_lst_dict[ 'eigval1' ][ -2 ] # The last step that converged before getting to negative eigenvalues
                    last_stable_step = disp_lst[ -2 ] 
                    apply_disp.t = last_stable_step
                    solver.solve()
                    v.vector()[:] += pf_dispfield * eigmode_critical.vector()[:]
                elif perturb_type == 'nan':
                    #apply_disp.t = eigval_lst_dict[ 'eigval1' ][ -1 ] 
                    apply_disp.t = disp_lst[-1]
                    solver.solve()

    end_time = time.time()
    run_time = end_time - start_time

    #-------Outputs--------
    # getting initial conf and final conf
    x_dofs = V.sub(0).sub(0).dofmap().dofs()
    y_dofs = V.sub(0).sub(1).dofmap().dofs()
    theta_dofs = V.sub(1).dofmap().dofs()
    dofs = V.tabulate_dof_coordinates()
    dof_coords = dofs.reshape((-1, 2))

    x_nodal_coord = dof_coords[x_dofs][:,0]
    y_nodal_coord = dof_coords[y_dofs][:,1]
    initconf_xy = np.concatenate( (x_nodal_coord.reshape(-1, 1),  y_nodal_coord.reshape(-1, 1)), axis = 1 )

    finalconf_x = x_nodal_coord + dispfield_lst[-1][x_dofs]
    finalconf_y = y_nodal_coord + dispfield_lst[-1][y_dofs]
    finalconf_xy = np.concatenate( (finalconf_x.reshape(-1, 1),  finalconf_y.reshape(-1, 1)), axis = -1)
    for i in range( N_eig_solve ):
        if eigmode_critical_dict[ f'eigmode_critical{i + 1}' ] is not None:
            eigmode_critical_vect = eigmode_critical_dict[ f'eigmode_critical{i + 1}' ].copy()
            eigmode_critical_dispx = eigmode_critical_vect[x_dofs]
            eigmode_critical_dispy = eigmode_critical_vect[y_dofs]
            eigmode_conf_xy = np.concatenate( (eigmode_critical_dispx.reshape(-1, 1),  eigmode_critical_dispy.reshape(-1, 1)), axis = -1)
            eigmode_critical_dict[ f'eigmode_critical{i + 1}' ] = eigmode_conf_xy.copy()


    dispfield_dict = { 'x_nodal_coord':x_nodal_coord, 'y_nodal_coord':y_nodal_coord,
                       'dispfield_lst':dispfield_lst, 'x_dofs':x_dofs, 'y_dofs':y_dofs }

    results_dict = {}
    # setting up a dictionary that contains all solution arrays
    results_dict['force_y'] = np.array( force_lst )
    results_dict['disp'] = np.array( disp_lst )
    for i in range( N_eig_solve ):
        results_dict[ f'eigval{ i + 1 }' ] = eigval_lst_dict[ f'eigval{i + 1}' ]
    for i in range( N_eig_solve ):
        results_dict[ f'eigmode_critical{ i + 1 }' ] = eigmode_critical_dict[ f'eigmode_critical{i + 1}' ]
    eigvals = np.concatenate( tuple([ np.array( eigval_lst_dict[ f'eigval{ i + 1 }' ]).reshape(1, -1) for i in range( N_eig_solve ) ]), axis = 0) 
    results_dict['eigvals'] = eigvals
    forcedispeigval_critic_tuple = (np.array(force_critical).reshape(1), np.array(disp_critical).reshape(1), np.array(eigval_critical).reshape(1) )
    forcedispeigval_critic = np.concatenate(forcedispeigval_critic_tuple )
    results_dict[ 'forcedispeigval_critic' ] = forcedispeigval_critic # an array that conatins critical force, displacemnt and eigen value
    results_dict[ 'force_critical' ] = np.array( [force_critical] )
    results_dict[ 'disp_critical' ] = np.array( [disp_critical] )
    results_dict[ 'eigval_critical' ] = np.array( [eigval_critical] )
    results_dict['initconf'] = initconf_xy
    results_dict['finalconf'] = finalconf_xy
    results_dict[ 'run_time' ] = np.array( [run_time] )

    print( 'disp field list size is', len( dispfield_lst ) )
    print( f" criticla force is { results_dict[ 'force_critical' ]  }")
    print( f" criticla disp is { results_dict[ 'disp_critical' ]  }")

    return results_dict, dispfield_dict
