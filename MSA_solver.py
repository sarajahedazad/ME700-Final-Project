import numpy as np
from scipy.linalg import eig
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix

def partition(matrix, known_DoF_idx):
    known_idx = known_DoF_idx
    n = matrix.shape[0]
    
    # Identify free DOF indices (those not prescribed)
    free_idx = np.setdiff1d(np.arange(n), known_idx)
    
    # Partition the stiffness matrix into submatrices:
    # K_ss: supported-supported, K_sf: supported-free,
    # K_fs: free-supported, K_ff: free-free.
    matrix_ss = matrix[np.ix_(known_idx, known_idx)]
    matrix_sf = matrix[np.ix_(known_idx, free_idx)]
    matrix_fs = matrix[np.ix_(free_idx, known_idx)]
    matrix_ff = matrix[np.ix_(free_idx, free_idx)]

    return matrix_ss, matrix_sf, matrix_fs, matrix_ff

def sparse_solve_stiffness_system(K, bcs):
    known_disp_idx = bcs.BCs_supported_indices 
    Delta_prescribed = bcs.BCs_Delta_supported_values

    known_F_idx = bcs.BCs_free_indices 
    F_prescribed = bcs.BCs_F_free_values
     
    F_full = np.zeros((bcs.n_DoFs))
    F_full[known_F_idx] = F_prescribed

    n = K.shape[0]
    K_ss, K_sf, K_fs, K_ff = partition(K, known_disp_idx)
    
    # Partition the displacement and force vectors:
    Delta_s = np.array(Delta_prescribed)
    F_f = F_full[known_F_idx]
    
    # For the free DOFs, the governing equations are:
    # F_f = K_fs * X_s + K_ff * X_f
    # Solve for the unknown displacements X_f:
    K_ff_sparse = csr_matrix(K_ff)
    K_fs_sparse = csr_matrix(K_fs)

    # Compute the right-hand side
    rhs = F_f - K_fs_sparse @ Delta_s

    # Solve the sparse system
    Delta_f = spsolve(K_ff_sparse, rhs)

    
    # Compute reaction forces at the supported DOFs:
    F_s = K_ss @ Delta_s + K_sf @ Delta_f
    
    # Assemble the complete displacement and force vectors:
    Delta = np.zeros(n)
    F = np.copy(F_full)
    Delta[known_disp_idx] = Delta_s
    Delta[known_F_idx] = Delta_f
    F[known_disp_idx] = F_s  # reaction forces
    
    return Delta, F


def solve_stiffness_system(K, bcs):
    known_disp_idx = bcs.BCs_supported_indices 
    Delta_prescribed = bcs.BCs_Delta_supported_values

    known_F_idx = bcs.BCs_free_indices 
    F_prescribed = bcs.BCs_F_free_values
     
    F_full = np.zeros((bcs.n_DoFs))
    F_full[known_F_idx] = F_prescribed

    n = K.shape[0]
    K_ss, K_sf, K_fs, K_ff = partition(K, known_disp_idx)
    
    # Partition the displacement and force vectors:
    Delta_s = np.array(Delta_prescribed)
    F_f = F_full[known_F_idx]
    
    # For the free DOFs, the governing equations are:
    # F_f = K_fs * X_s + K_ff * X_f
    # Solve for the unknown displacements X_f:
    Delta_f = np.linalg.solve(K_ff, F_f - K_fs @ Delta_s)
    
    # Compute reaction forces at the supported DOFs:
    F_s = K_ss @ Delta_s + K_sf @ Delta_f
    
    # Assemble the complete displacement and force vectors:
    Delta = np.zeros(n)
    F = np.copy(F_full)
    Delta[known_disp_idx] = Delta_s
    Delta[known_F_idx] = Delta_f
    F[known_disp_idx] = F_s  # reaction forces
    
    return Delta, F

def compute_critical_load(K_elastic, K_geometric, bcs):
    n = K_elastic.shape[0] 

    known_DoF_idx = bcs.BCs_supported_indices 
    k_ss, _, _, k_ff = partition(K_elastic, known_DoF_idx) 
    k_g_ss, _, _, k_g_ff = partition(K_geometric, known_DoF_idx) 

    # Solve the generalized eigenvalue problem: (K_elastic - lambda * K_geometric) = 0
    eigenvalues, eigenvectors = eig(k_ff, -k_g_ff)

    #print( eigenvalues )

    # Define a tolerance for checking the "realness" of an eigenvalue.
    tol = 1e-8

    # Filter indices for eigenvalues whose imaginary part is negligible
    real_indices = np.where(np.abs(np.imag(eigenvalues)) < tol)[0]

    # Extract the real parts of the eigenvalues and the corresponding eigenvectors.
    real_eigenvalues = np.real(eigenvalues[real_indices])
    real_eigenvectors = eigenvectors[:, real_indices]

    # Further filter to select only the positive eigenvalues.
    positive_indices = np.where(real_eigenvalues > 0)[0]
    if positive_indices.size == 0:
        raise ValueError("No positive real eigenvalue found.")

    positive_eigenvalues = real_eigenvalues[positive_indices]
    positive_eigenvectors = real_eigenvectors[:, positive_indices]

    # Identify the index of the smallest positive eigenvalue.
    min_index = np.argmin(positive_eigenvalues)

    # Retrieve the smallest positive eigenvalue and its corresponding eigenvector.
    smallest_positive_eigenvalue = positive_eigenvalues[min_index]
    corresponding_eigenvector = positive_eigenvectors[:, min_index]

    eigenvectors_allstructure = np.zeros(n)
    eigenvectors_allstructure[bcs.BCs_free_indices] = corresponding_eigenvector

    # # Sort indices based on the eigenvalues
    # sorted_indices = np.argsort(positive_eigenvalues)

    # # Ensure that there are at least two positive eigenvalues
    # if sorted_indices.size < 2:
    #     raise ValueError("Fewer than two positive eigenvalues are available.")

    # # Retrieve the index for the second smallest eigenvalue
    # second_min_index = sorted_indices[1]

    # # Use the second smallest index to retrieve eigenvalue and eigenvector
    # second_smallest_positive_eigenvalue = positive_eigenvalues[second_min_index]
    # corresponding_eigenvector = positive_eigenvectors[:, second_min_index]

    # # Optionally print properties of the imaginary part of the eigenvector
    # print(np.max(np.imag(corresponding_eigenvector)), 
    #     np.min(np.imag(corresponding_eigenvector)), 
    #     np.linalg.norm(np.imag(corresponding_eigenvector)))

    # # Assign the eigenvector values to the full structure
    # eigenvectors_allstructure = np.zeros(n)
    # eigenvectors_allstructure[bcs.BCs_free_indices] = corresponding_eigenvector


    return smallest_positive_eigenvalue, eigenvectors_allstructure
    #return second_smallest_positive_eigenvalue, eigenvectors_allstructure
