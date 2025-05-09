import numpy as np

class OverDefinedError(Exception):
    pass

class UnderDefinedError(Exception):
    pass

class BoundaryConditions2D:
    def __init__(self, frame_obj):
        self.frame = frame_obj
        self.n_points = self.frame.points.shape[0]
        self.n_DoFs = self.n_points * 3

        # displacement bounds
        self.BCs_disp_indices = []
        self.BCs_disp_values = []

        # rotation bounds
        self.BCs_rot_indices = []
        self.BCs_rot_values = []

        # force bounds
        self.BCs_force_indices = []
        self.BCs_force_values = []

        # moment bounds
        self.BCs_moment_indices = []
        self.BCs_moment_values = []

        # consolidated boundary conditions
        self.BCs_supported_indices = []
        self.BCs_Delta_supported_values = []
        self.BCs_free_indices = []
        self.BCs_F_free_values = []

    def add_disp_bound_xy(self, point_coords, value_x, value_y):
        idx = self.frame.points_unique[tuple(point_coords)]
        self.BCs_disp_indices += [idx * 3, idx * 3 + 1]
        self.BCs_disp_values += [value_x, value_y]

    def add_disp_bound_x(self, point_coords, value: float):
        idx = self.frame.points_unique[tuple(point_coords)]
        self.BCs_disp_indices += [idx * 3]
        self.BCs_disp_values += [value]

    def add_disp_bound_y(self, point_coords, value: float):
        idx = self.frame.points_unique[tuple(point_coords)]
        self.BCs_disp_indices += [idx * 3 + 1]
        self.BCs_disp_values += [value]

    def add_rot_bound_z(self, point_coords, value: float):
        idx = self.frame.points_unique[tuple(point_coords)]
        self.BCs_rot_indices += [idx * 3 + 2]
        self.BCs_rot_values += [value]

    def add_force_bound_xy(self, point_coords, value_x, value_y):
        idx = self.frame.points_unique[tuple(point_coords)]
        self.BCs_force_indices += [idx * 3, idx * 3 + 1]
        self.BCs_force_values += [value_x, value_y]

    def add_force_bound_x(self, point_coords, value):
        idx = self.frame.points_unique[tuple(point_coords)]
        self.BCs_force_indices += [idx * 3]
        self.BCs_force_values += [value]

    def add_force_bound_y(self, point_coords, value):
        idx = self.frame.points_unique[tuple(point_coords)]
        self.BCs_force_indices += [idx * 3 + 1]
        self.BCs_force_values += [value]

    def add_moment_bound_z(self, point_coords, value):
        idx = self.frame.points_unique[tuple(point_coords)]
        self.BCs_moment_indices += [idx * 3 + 2]
        self.BCs_moment_values += [value]

    def validate_bounds(self):
        intersect_len = len(np.intersect1d(self.BCs_supported_indices, self.BCs_free_indices))
        union_len = len(np.union1d(self.BCs_supported_indices, self.BCs_free_indices))
        if intersect_len > 0:
            raise OverDefinedError("The problem is overdefined.")
        elif union_len < self.n_DoFs:
            raise UnderDefinedError("The problem is underdefined.")
        else:
            print("Bounds are good to go!")

    def set_up_bounds(self):
        # Known supported indices (Dirichlet BCs)
        BCs_supported_indices = np.concatenate((self.BCs_disp_indices, self.BCs_rot_indices)).astype(int)
        BCs_Delta_supported_values = np.concatenate((self.BCs_disp_values, self.BCs_rot_values))
        order = np.argsort(BCs_supported_indices)
        self.BCs_supported_indices = BCs_supported_indices[order]
        self.BCs_Delta_supported_values = BCs_Delta_supported_values[order]

        # Known free indices (Neumann BCs)
        BCs_free_indices = np.concatenate((self.BCs_force_indices, self.BCs_moment_indices)).astype(int)
        BCs_F_free_values = np.concatenate((self.BCs_force_values, self.BCs_moment_values))
        order = np.argsort(BCs_free_indices)
        self.BCs_free_indices = BCs_free_indices[order]
        self.BCs_F_free_values = BCs_F_free_values[order]

        self.validate_bounds()

    # --- New generalized methods for 2D ---

    def apply_boundary_conditions_by_direction(self, direction, coordinate_value, disp_vals, rot_vals, tol=1e-8):
        """
        Automatically assigns displacement and rotation boundary conditions to all nodes
        whose coordinate in the specified direction is equal to coordinate_value (within a tolerance).
        
        Parameters:
            direction (str): One of 'x' or 'y' (for 2D).
            coordinate_value (float): The coordinate value to match.
            disp_vals (tuple or list of 2 floats): Prescribed displacement values (d_x, d_y).
            rot_vals (float or a one-element list/tuple): Prescribed rotation value (theta_z).
            tol (float): Tolerance for the equality check (default is 1e-8).
        """
        direction = direction.lower()
        if direction == 'x':
            coord_index = 0
        elif direction == 'y':
            coord_index = 1
        else:
            raise ValueError("Invalid direction for 2D. Choose from 'x' or 'y'.")

        # Update node count and DOFs in case the frame has been meshed.
        self.n_points = self.frame.points.shape[0]
        self.n_DoFs = self.n_points * 3

        # Loop over all nodes in the meshed frame.
        for node_idx, coord in enumerate(self.frame.points):
            if np.isclose(coord[coord_index], coordinate_value, atol=tol):
                # Assign displacement conditions (DOFs: node_idx*3, node_idx*3+1).
                self.BCs_disp_indices += [node_idx * 3, node_idx * 3 + 1]
                self.BCs_disp_values += list(disp_vals)
                # Assign rotation condition (DOF: node_idx*3+2).
                self.BCs_rot_indices.append(node_idx * 3 + 2)
                if isinstance(rot_vals, (list, tuple)):
                    self.BCs_rot_values += list(rot_vals)
                else:
                    self.BCs_rot_values.append(rot_vals)
        print(f"Boundary conditions applied for nodes with {direction} = {coordinate_value}")

    def apply_free_conditions(self):
        """
        Automatically assigns zero force and zero moment boundary conditions to every DOF
        not already constrained by prescribed displacement/rotation (Dirichlet) or force/moment (Neumann)
        boundary conditions. For 2D, each node has 3 DOFs:
         - DOFs 0 and 1 (translation): zero force is assigned.
         - DOF 2 (rotation): zero moment is assigned.
        """
        total_dofs = self.n_DoFs
        # Combine all DOFs already constrained by any boundary condition.
        constrained_dofs = set(self.BCs_disp_indices + self.BCs_rot_indices +
                               self.BCs_force_indices + self.BCs_moment_indices)
        for dof in range(total_dofs):
            if dof not in constrained_dofs:
                if dof % 3 < 2:  # Translation DOFs (0 and 1)
                    self.BCs_force_indices.append(dof)
                    self.BCs_force_values.append(0.0)
                else:           # Rotation DOF (2)
                    self.BCs_moment_indices.append(dof)
                    self.BCs_moment_values.append(0.0)
        print("Free conditions (zero force and zero moment) applied to all unconstrained DOFs.")

    def get_force_on_boundary( self, direction, coordinate_value, F):
        direction = direction.lower()
        if direction == 'x':
            coord_index = 0
        elif direction == 'y':
            coord_index = 1
        else:
            raise ValueError("Invalid direction for 2D. Choose from 'x' or 'y'.")

        points = self.frame.points
        node_indices_on_boundary, = np.where( points[:, coord_index] == coordinate_value )
        force_dof_indices = node_indices_on_boundary * 3 + 1

        total_boundary_force = F[force_dof_indices].sum()
        return total_boundary_force

