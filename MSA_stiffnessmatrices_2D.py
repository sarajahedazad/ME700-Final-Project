import numpy as np

def local_elastic_stiffness_matrix_2D_beam(E: float, A: float, L: float, I: float) -> np.ndarray:
    k_e = np.zeros((6, 6))
    
    # Axial stiffness terms
    axial = E * A / L
    k_e[0, 0] = axial
    k_e[0, 3] = -axial
    k_e[3, 0] = -axial
    k_e[3, 3] = axial

    # Bending stiffness terms
    bending_factor = E * I / L**3
    k_e[1, 1] = 12 * bending_factor
    k_e[1, 4] = -12 * bending_factor
    k_e[4, 1] = -12 * bending_factor
    k_e[4, 4] = 12 * bending_factor

    k_e[1, 2] = 6 * L * bending_factor
    k_e[2, 1] = 6 * L * bending_factor
    k_e[1, 5] = 6 * L * bending_factor
    k_e[5, 1] = 6 * L * bending_factor

    k_e[2, 4] = -6 * L * bending_factor
    k_e[4, 2] = -6 * L * bending_factor
    k_e[4, 5] = -6 * L * bending_factor
    k_e[5, 4] = -6 * L * bending_factor

    k_e[2, 2] = 4 * L**2 * bending_factor
    k_e[5, 5] = 4 * L**2 * bending_factor
    k_e[2, 5] = 2 * L**2 * bending_factor
    k_e[5, 2] = 2 * L**2 * bending_factor

    return k_e

def check_unit_vector(vec: np.ndarray):
    if np.isclose(np.linalg.norm(vec), 1.0):
        return
    else:
        raise ValueError("Expected a unit vector for reference vector (2D).")

def check_parallel(vec_1: np.ndarray, vec_2: np.ndarray):
    if np.isclose(np.cross(vec_1, vec_2), 0.0):
        raise ValueError("Reference vector is parallel to beam axis (2D).")
    else:
        return

def rotation_matrix_2D(x1: float, y1: float, x2: float, y2: float):
    L = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    c = (x2 - x1) / L
    s = (y2 - y1) / L
    
    gamma = np.array([
        [ c,  s, 0],
        [-s,  c, 0],
        [ 0,  0, 1]
    ])
    
    return gamma

def transformation_matrix_2D(gamma: np.ndarray) -> np.ndarray:
    Gamma = np.zeros((6, 6))
    Gamma[0:3, 0:3] = gamma
    Gamma[3:6, 3:6] = gamma
    return Gamma

def local_geometric_stiffness_matrix_2D_beam(L, Fx):
    k_g = np.zeros((6, 6))

    # Axial terms
    k_g[0, 0] = Fx / L
    k_g[0, 3] = -Fx / L
    k_g[3, 0] = -Fx / L
    k_g[3, 3] = Fx / L

    # Bending terms
    k_g[1, 1] = 6 * Fx / (5 * L)
    k_g[1, 2] = Fx / 10
    k_g[1, 4] = -6 * Fx / (5 * L)
    k_g[1, 5] = Fx / 10

    k_g[2, 1] = Fx / 10
    k_g[2, 2] = 2 * Fx * L / 15
    k_g[2, 4] = -Fx / 10
    k_g[2, 5] = -Fx * L / 30

    k_g[4, 1] = -6 * Fx / (5 * L)
    k_g[4, 2] = -Fx / 10
    k_g[4, 4] = 6 * Fx / (5 * L)
    k_g[4, 5] = -Fx / 10

    k_g[5, 1] = Fx / 10
    k_g[5, 2] = -Fx * L / 30
    k_g[5, 4] = -Fx / 10
    k_g[5, 5] = 2 * Fx * L / 15

    return k_g


#-------------------------------
#-------------------------------
#-------------------------------
class StiffnessMatrices2D:
    def __init__(self, frame_obj):
        self.frame = frame_obj
        self.n_points = self.frame.points.shape[0]
        self.n_connectivities = self.frame.connectivities.shape[0]
        self.n_DoFs = self.n_points * 3

    def get_element_parameters(self, element_idx):
        connection = self.frame.connectivities[element_idx]
        E = self.frame.E_array[element_idx]
        A = self.frame.A_array[element_idx]
        L = self.frame.L_array[element_idx]
        I = self.frame.I_array[element_idx]
        return connection, E, A, L, I

    def get_element_points(self, connection):
        p0_idx = connection[0]
        p0 = self.frame.points[p0_idx]
        p1_idx = connection[1]
        p1 = self.frame.points[p1_idx]
        return p0_idx, p0, p1_idx, p1

    def get_transformation_matrix_2D(self, p0, p1):
        gamma = rotation_matrix_2D(p0[0], p0[1], p1[0], p1[1])
        Gamma = transformation_matrix_2D(gamma)
        return Gamma

    def get_global_elastic_stiffmatrix(self):
        K = np.zeros((self.n_DoFs, self.n_DoFs))
        for element_idx in range(self.n_connectivities):
            connection, E, A, L, I = self.get_element_parameters(element_idx)
            p0_idx, p0, p1_idx, p1 = self.get_element_points(connection)

            k_element_local = local_elastic_stiffness_matrix_2D_beam(E, A, L, I)
            Gamma = self.get_transformation_matrix_2D(p0, p1)
            k_element_global = Gamma.T @ k_element_local @ Gamma

            p0_DoF_idx = p0_idx * 3
            p1_DoF_idx = p1_idx * 3

            K[p0_DoF_idx:p0_DoF_idx + 3, p0_DoF_idx:p0_DoF_idx + 3] += k_element_global[0:3, 0:3]
            K[p0_DoF_idx:p0_DoF_idx + 3, p1_DoF_idx:p1_DoF_idx + 3] += k_element_global[0:3, 3:6]
            K[p1_DoF_idx:p1_DoF_idx + 3, p0_DoF_idx:p0_DoF_idx + 3] += k_element_global[3:6, 0:3]
            K[p1_DoF_idx:p1_DoF_idx + 3, p1_DoF_idx:p1_DoF_idx + 3] += k_element_global[3:6, 3:6]

        return K

    def get_element_local_internal_F(self, element_idx, Delta):
        connection, E, A, L, I = self.get_element_parameters(element_idx)
        p0_idx, p0, p1_idx, p1 = self.get_element_points(connection)

        Delta_el_global = np.concatenate((Delta[p0_idx * 3:p0_idx * 3 + 3], Delta[p1_idx * 3:p1_idx * 3 + 3]))
        Gamma = self.get_transformation_matrix_2D(p0, p1)
        Delta_el_local = Gamma @ Delta_el_global

        K_el_local = local_elastic_stiffness_matrix_2D_beam(E, A, L, I)
        F_el_local = K_el_local @ Delta_el_local

        return F_el_local

    def get_element_local_geometric_stiffness_matrix(self, element_idx, Delta):
        connection, E, A, L, I = self.get_element_parameters(element_idx)
        F_el_local = self.get_element_local_internal_F(element_idx, Delta)
        Fx = F_el_local[3]  # axial force at node 2
        k_geom_element_local = local_geometric_stiffness_matrix_2D_beam(L, Fx)
        return k_geom_element_local

    def get_global_geometric_stiffmatrix(self, Delta):
        K_g = np.zeros((self.n_DoFs, self.n_DoFs))
        for element_idx in range(self.n_connectivities):
            connection, E, A, L, I = self.get_element_parameters(element_idx)
            p0_idx, p0, p1_idx, p1 = self.get_element_points(connection)

            k_g_element_local = self.get_element_local_geometric_stiffness_matrix(element_idx, Delta)
            Gamma = self.get_transformation_matrix_2D(p0, p1)
            k_g_element_global = Gamma.T @ k_g_element_local @ Gamma

            p0_DoF_idx = p0_idx * 3
            p1_DoF_idx = p1_idx * 3

            K_g[p0_DoF_idx:p0_DoF_idx + 3, p0_DoF_idx:p0_DoF_idx + 3] += k_g_element_global[0:3, 0:3]
            K_g[p0_DoF_idx:p0_DoF_idx + 3, p1_DoF_idx:p1_DoF_idx + 3] += k_g_element_global[0:3, 3:6]
            K_g[p1_DoF_idx:p1_DoF_idx + 3, p0_DoF_idx:p0_DoF_idx + 3] += k_g_element_global[3:6, 0:3]
            K_g[p1_DoF_idx:p1_DoF_idx + 3, p1_DoF_idx:p1_DoF_idx + 3] += k_g_element_global[3:6, 3:6]

        return K_g