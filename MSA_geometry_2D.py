import numpy as np

class DuplicationError(Exception):
    pass

class Point2D:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.coords = np.array([self.x, self.y])

class Element2D:
    def __init__(self, p0, p1, E, A, I):
        self.p0 = p0
        self.p1 = p1
        self.E = E
        self.A = A
        self.I = I
        self.L = self.calc_connection_length()

    def calc_connection_length(self):
        length = np.sqrt((self.p0.x - self.p1.x)**2 + (self.p0.y - self.p1.y)**2)
        return length

class Frame2D:
    def __init__(self):
        self.points = None
        self.connectivities = None
        self.E_array = None
        self.A_array = None
        self.I_array = None
        self.L_array = None
        self.points_unique = {}

    def add_point(self, x, y):
        return Point2D(x, y)

    def add_element(self, p0, p1, E, A, I):
        return Element2D(p0, p1, E, A, I)

    def build_frame(self, element_lst):
        self.points_unique = {}
        point_lst = []
        connectivity_lst = []
        E_lst, A_lst, I_lst, L_lst = [], [], [], []

        for element in element_lst:
            p0_obj = element.p0
            p1_obj = element.p1
            p0 = (p0_obj.x, p0_obj.y)
            p1 = (p1_obj.x, p1_obj.y)

            if p0 in self.points_unique:
                p0_idx = self.points_unique[p0]
            else:
                p0_idx = len(self.points_unique)
                self.points_unique[p0] = p0_idx
                point_lst.append(p0)

            if p1 in self.points_unique:
                p1_idx = self.points_unique[p1]
            else:
                p1_idx = len(self.points_unique)
                self.points_unique[p1] = p1_idx
                point_lst.append(p1)

            if ([p0_idx, p1_idx] in connectivity_lst) or ([p1_idx, p0_idx] in connectivity_lst):
                raise DuplicationError('Duplication in defining the elements! Be careful!')

            connectivity_lst.append([p0_idx, p1_idx])
            E_lst.append(element.E)
            A_lst.append(element.A)
            I_lst.append(element.I)
            L_lst.append(element.L)

        self.points = np.array(point_lst)
        self.connectivities = np.sort(np.array(connectivity_lst), axis=1)
        self.E_array = np.array(E_lst)
        self.A_array = np.array(A_lst)
        self.I_array = np.array(I_lst)
        self.L_array = np.array(L_lst)

        print('Your 2D frame is ready!')

    def calc_all_connections_lengths(self):
        lengths = []
        for conn in self.connectivities:
            p0_idx, p1_idx = conn
            p0 = self.points[p0_idx]
            p1 = self.points[p1_idx]
            length = np.sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2)
            lengths.append(length)
        return np.array(lengths)

    def generate_frame_from_mesh( self, mesh, E, A, I ):
        self.points = mesh.coordinates()
        self.connectivities = np.sort(mesh.cells(), axis=1) # check this
        self.E_array = E * np.ones( ( self.connectivities.shape[0], 1))
        self.A_array = A * np.ones( ( self.connectivities.shape[0], 1))
        self.I_array = I * np.ones( ( self.connectivities.shape[0], 1))
        self.L_array = self.calc_all_connections_lengths()


    def generate_frame_directly(self, points, connectivities, E_array, A_array, I_array):
        self.points = points
        self.connectivities = np.sort(connectivities, axis=1)
        self.E_array = E_array
        self.A_array = A_array
        self.I_array = I_array
        self.L_array = self.calc_all_connections_lengths()

    def mesh_frame(self, characteristic_length):
        """
        Mesh every element by inserting intermediate nodes along each element
        based on the given characteristic_length. The meshing process creates new nodes
        on each element by interpolation in 2D and updates the nodes and connectivity.
        """
        new_points_list = []         # New node coordinates as tuples (x, y)
        new_connectivities_list = [] # New sub-element connectivity (pairs of node indices)
        new_E, new_A, new_I, new_L = [], [], [], []

        # Global dictionary to store unique nodes (using rounded coordinates to avoid floating-point issues)
        global_nodes = {}

        def get_or_add_node(coord):
            # Round coordinates for consistency; we assume a 2D point here.
            coord_rounded = tuple(np.round(coord, 8).tolist())
            if coord_rounded in global_nodes:
                return global_nodes[coord_rounded]
            else:
                idx = len(new_points_list)
                new_points_list.append(coord_rounded)
                global_nodes[coord_rounded] = idx
                return idx

        # Loop over each original element.
        for idx, conn in enumerate(self.connectivities):
            s_idx, e_idx = conn
            start = self.points[s_idx]
            end   = self.points[e_idx]

            # Compute the Euclidean distance between start and end in 2D.
            distance = np.linalg.norm(end - start)
            n_segments = int(np.ceil(distance / characteristic_length))
            if n_segments < 1:
                n_segments = 1

            seg_node_indices = []
            # Create interpolated nodes along the element.
            for j in range(n_segments + 1):
                t = j / n_segments
                p = start + t * (end - start)
                node_idx = get_or_add_node(p)
                seg_node_indices.append(node_idx)

            # Create sub-elements (connectivities) between successive nodes.
            for k in range(len(seg_node_indices) - 1):
                n1 = seg_node_indices[k]
                n2 = seg_node_indices[k + 1]
                new_connectivities_list.append([n1, n2])
                new_E.append(self.E_array[idx])
                new_A.append(self.A_array[idx])
                new_I.append(self.I_array[idx])
                p1 = np.array(new_points_list[n1])
                p2 = np.array(new_points_list[n2])
                new_L.append(np.linalg.norm(p2 - p1))

        # Update the frame with the meshed data.
        self.points = np.array(new_points_list)
        self.connectivities = np.array(new_connectivities_list)
        self.E_array = np.array(new_E)
        self.A_array = np.array(new_A)
        self.I_array = np.array(new_I)
        self.L_array = np.array(new_L)

        print("Meshing completed with characteristic length =", characteristic_length)