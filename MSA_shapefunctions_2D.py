import numpy as np
import matplotlib.pyplot as plt
import sympy
from MSA_geometry_2D import *
from MSA_stiffnessmatrices_2D import * 

def interpolate_two_points(p0, p1, n):
    points_lst = [tuple(p0 + (p1 - p0) * i / n) for i in range(n + 1)]
    return np.array(points_lst)

class ShapeFunctions2D:
    def __init__(self, eigenvector, frame_obj, n=20, scale=5):
        self.eigenvector = eigenvector
        self.frame = frame_obj
        self.n = n
        self.scale = scale
        self.x_local = sympy.symbols(f'x:{self.n}')

    def transformation_66_matrix_2D(self, gamma):
        Gamma = np.zeros((6, 6))
        Gamma[0:3, 0:3] = gamma
        Gamma[3:6, 3:6] = gamma
        return Gamma

    def transformation_nn_matrix_2D(self, gamma):
        Gamma = np.zeros((2 * self.n, 2 * self.n))
        for i in range(self.n):
            Gamma[2 * i: 2 * i + 2, 2 * i: 2 * i + 2] = gamma[:2, :2]
        return Gamma

    def evaluate(self, expr, x0):
        func = sympy.lambdify(self.x_local, expr, 'numpy')
        func_val = func(*x0)
        return func_val.reshape(-1)

    def linear_N1(self, length):
        expr = sympy.Matrix([1 - self.x_local[i] / length for i in range(self.n)])
        return expr

    def linear_N2(self, length):
        expr = sympy.Matrix([self.x_local[i] / length for i in range(self.n)])
        return expr

    def hermite_N1(self, length):
        expr = sympy.Matrix([1 - 3 * (self.x_local[i] / length) ** 2 + 2 * (self.x_local[i] / length) ** 3 for i in range(self.n)])
        return expr

    def hermite_N2(self, length):
        expr = sympy.Matrix([3 * (self.x_local[i] / length) ** 2 - 2 * (self.x_local[i] / length) ** 3 for i in range(self.n)])
        return expr 

    def hermite_N3(self, length):
        expr = sympy.Matrix([self.x_local[i] * (1 - self.x_local[i] / length) ** 2 for i in range(self.n)])
        return expr

    def hermite_N4(self, length):
        expr = sympy.Matrix([self.x_local[i] * ((self.x_local[i] / length) ** 2 - self.x_local[i] / length) for i in range(self.n)])
        return expr

    def get_element_info(self, element_idx):
        connection = self.frame.connectivities[element_idx]
        p0_idx, p1_idx = connection
        p0 = self.frame.points[p0_idx]
        p1 = self.frame.points[p1_idx]
        length = self.frame.L_array[element_idx]
        return connection, p0_idx, p0, p1_idx, p1, length

    def get_eigenvector_element_global(self, p0_idx, p1_idx):
        return np.concatenate((self.eigenvector[3 * p0_idx: 3 * p0_idx + 3], self.eigenvector[3 * p1_idx: 3 * p1_idx + 3]))

    def calc_element_interpolation(self, element_idx):
        connection, p0_idx, p0, p1_idx, p1, length = self.get_element_info(element_idx)
        gamma = rotation_matrix_2D(p0[0], p0[1], p1[0], p1[1])
        Gamma = self.transformation_66_matrix_2D(gamma)
        eigenvector_el_global = self.get_eigenvector_element_global(p0_idx, p1_idx)
        eigenvector_el_local = Gamma @ eigenvector_el_global

        u_p0, v_p0, theta_p0, u_p1, v_p1, theta_p1 = eigenvector_el_local
        x_local_val = np.linspace(0, length, num=self.n)

        u_local = u_p0 * self.evaluate(self.linear_N1(length), x_local_val) + u_p1 * self.evaluate(self.linear_N2(length), x_local_val)
        v_local_part1 = v_p0 * self.evaluate(self.hermite_N1(length), x_local_val) + v_p1 * self.evaluate(self.hermite_N2(length), x_local_val)
        v_local_part2 = theta_p0 * self.evaluate(self.hermite_N3(length), x_local_val) + theta_p1 * self.evaluate(self.hermite_N4(length), x_local_val)
        v_local = v_local_part1 + v_local_part2

        uv_local = np.stack((self.scale * u_local, self.scale * v_local)).flatten('F')
        Gamma_nn = self.transformation_nn_matrix_2D(gamma)
        uv_global_element = Gamma_nn.T @ uv_local

        interpolated_points = interpolate_two_points(p0, p1, self.n - 1)
        return interpolated_points + uv_global_element.reshape(-1, 2)

    def plot_element_interpolation(self, saving_dir_with_name, plot_title = '', dpi=500):
        points = self.frame.points
        connectivities = self.frame.connectivities

        plt.figure()
        for connection in connectivities:
            p0 = points[connection[0]]
            p1 = points[connection[1]]
            plt.plot([p0[0], p1[0]], [p0[1], p1[1]], 'k--', label='Elements' if connection is connectivities[0] else "")

        for j, connection in enumerate(connectivities):
            element_interpolated = self.calc_element_interpolation(j)
            plt.plot(element_interpolated[:, 0], element_interpolated[:, 1], 'r', label='Interpolated Shape' if j == 0 else "")

        plt.xlabel('X')
        plt.ylabel('Y')
        plt.legend()
        plt.axis('equal')
        #plt.grid(True)
        plt.title( plot_title )
        plt.savefig(saving_dir_with_name, dpi=dpi)
        plt.show()
