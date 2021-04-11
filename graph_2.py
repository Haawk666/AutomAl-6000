# By Haakon Tvedt @ NTNU
# Contributors:
"""This module contains the basic graph components of atomic graphs."""

# Internal imports:
import utils
import column_characterization
# External imports:
import numpy as np
import copy
import random
import logging
# Instantiate logger:
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Vertex:
    """A vertex is object is an atomic graph component.

    A vertex contains a lot of the local information of the parent atomic graph.

    :ivar parent_graph: A reference to the AtomicGraph object that the vertex belongs to.
    :vartype parent_graph: :code:`graph_2.AtomicGraph`

    :param index: Index of the vertex. Must be unique.
    :param im_coor_x: x-coordinate of the center of the column in the image which the vertex represents.
    :param im_coor_y: y-coordinate of the center of the column in the image which the vertex represents.
    :param r: The "approximate atomic radii" in pixels.
    :param scale: The scale of the HAADF-STEM image in pm/pixel
    :param zeta: (Optional) The initial zeta value of the vertex (0 or 1). Defaults to 0.
    :param advanced_species: (Optional) The initial "advanced species category" of the vertex. Defaults to :code:`'Un_1'`.
    :param atomic_species: (Optional) The initial "atomic species category" of the vertex. Defaults to :code:`'Un'`.
    :param void: (Optional) Whether the vertex is void or not. Defaults to :code:`False`.
    :param parent_graph: (Optional) A reference to the graph object that the vertex belongs to. Defaults to :code:`None`.
    :type parent_graph: :code:`graph_2.AtomicGraph`
    :type void: bool
    :type atomic_species: str
    :type advanced_species: str
    :type zeta: int
    :type scale: float
    :type r: int
    :type im_coor_y: float
    :type im_coor_x: float
    :type index: int

    """

    al_lattice_const = 404.95
    empty_map = {
        'district': [],
        'out_neighbourhood': set(),
        'in_neighbourhood': set(),
        'neighbourhood': set(),
        'partners': set(),
        'semi_partners': set(),
        'out_semi_partners': set(),
        'in_semi_partners': set(),
        'anti_out_neighbourhood': set(),
        'anti_in_neighbourhood': set(),
        'anti_neighbourhood': set(),
        'anti_partners': set(),
        'anti_semi_partners': set(),
        'anti_out_semi_partners': set(),
        'anti_in_semi_partners': set()
    }
    numerical_attributes = {
        'peak_gamma',
        'avg_gamma',
        'area_gamma',
        'normalized_peak_gamma',
        'normalized_avg_gamma',
        'normalized_area_gamma',
        'intensity',
        'im_coor_x',
        'im_coor_y',
        'im_coor_z',
        'spatial_coor_x',
        'spatial_coor_y',
        'spatial_coor_z',
        'alpha_min',
        'alpha_max',
        'theta_min',
        'theta_max',
        'theta_angle_variance',
        'theta_angle_mean',
        'redshift'
    }
    nominal_attributes = {
        'i',
        'zeta',
        'atomic_species',
        'advanced_species',
        'n',
        'flag_1',
        'flag_2',
        'flag_3',
        'flag_4',
        'flag_5',
        'flag_6',
        'flag_7',
        'flag_8',
        'flag_9'
    }

    def __init__(self, index, im_coor_x, im_coor_y, r, scale, zeta=0, advanced_species='Un_1', n=3, atomic_species='Un', void=False, parent_graph=None):

        # parent
        self.parent_graph = parent_graph

        # Index
        self.i = index

        # Some properties
        self.r = r
        self.zeta = zeta
        self.peak_gamma = 0
        self.avg_gamma = 0
        self.area_gamma = 0
        self.intensity = 0
        self.advanced_species = advanced_species
        self.atomic_species = atomic_species

        # Position
        self.scale = scale
        self.zeta = zeta
        self.im_coor_x = im_coor_x
        self.im_coor_y = im_coor_y
        self.im_coor_z = zeta
        self.spatial_coor_x = im_coor_x * scale
        self.spatial_coor_y = im_coor_y * scale
        self.spatial_coor_z = zeta * 0.5 * self.al_lattice_const

        # Some settings
        self.is_in_precipitate = False
        self.is_edge_column = False
        self.is_set_by_user = False
        self.show_in_overlay = True
        self.void = void
        self.flag_1 = False
        self.flag_2 = False
        self.flag_3 = False
        self.flag_4 = False
        self.flag_5 = False
        self.flag_6 = False
        self.flag_7 = False
        self.flag_8 = False
        self.flag_9 = False

        # model-analysis
        self.probability_vector = {}
        self.alpha_probability_vector = {}
        self.composite_probability_vector = {}
        self.weighted_probability_vector = {}

        self.confidence = 0
        self.alpha_confidence = 0
        self.composite_confidence = 0
        self.weighted_confidence = 0

        self.advanced_probability_vector = {}
        self.advanced_alpha_probability_vector = {}
        self.advanced_composite_probability_vector = {}
        self.advanced_weighted_probability_vector = {}

        # Model variables
        self.alpha_angles = []
        self.alpha_max = 0
        self.alpha_min = 0
        self.theta_angles = []
        self.theta_max = 0
        self.theta_min = 0
        self.theta_angle_variance = 0
        self.theta_angle_mean = 0
        self.normalized_peak_gamma = 0
        self.normalized_avg_gamma = 0
        self.normalized_area_gamma = 0
        self.redshift = 0
        self.redshift_variance = 0
        self.avg_redshift = 0
        self.avg_central_separation = 0

        # Local graph mapping
        self.district = []
        self.city = set()

        self.out_neighbourhood = set()
        self.in_neighbourhood = set()
        self.neighbourhood = set()
        self.anti_neighbourhood = set()
        self.partners = set()
        self.semi_partners = set()
        self.out_semi_partners = set()
        self.in_semi_partners = set()

        self.projected_separation_district = []
        self.zeta_district = []

        self.local_zeta_map = copy.deepcopy(Vertex.empty_map)

        # Graph parameters
        self.n = n
        self.in_degree = 0
        self.out_degree = 0
        self.degree = 0

        self.reset_probability_vector_from_atomic_species()

    def __str__(self):
        return self.report()

    def report(self):
        """Create a report summarizing the contents of the vertex instance and return it as a string

        :returns: A vertex summary
        :rtype: str

        """
        im_pos = self.im_pos()
        spatial_pos = self.spatial_pos()
        string = 'Vertex {}:\n'.format(self.i)
        string += '    General:\n'
        string += '        Image position (x, y) = ({:.3f}, {:.3f})\n'.format(im_pos[0], im_pos[1])
        string += '        Pixel position (x, y) = ({:.0f}, {:.0f})\n'.format(np.floor(im_pos[0]), np.floor(im_pos[1]))
        string += '        Spatial relative position in pm (x, y, z) = ({:.3f}, {:.3f}, {:.3f})\n'.format(spatial_pos[0], spatial_pos[1], spatial_pos[2])
        string += '        Peak gamma = {:.4f}\n'.format(self.peak_gamma)
        string += '        Average gamma = {:.4f}\n'.format(self.avg_gamma)
        string += '        Atomic species: {}\n'.format(self.atomic_species)
        string += '        Advanced species: {}\n'.format(self.advanced_species)
        string += '        Symmetry: {}\n'.format(self.n)
        string += '    Analysis:\n'
        string += '        Probability vector: {\n'
        for key, prob in self.probability_vector.items():
            string += '            {}: {:.3f}\n'.format(key, prob)
        string += '        }\n'
        string += '            Prediction: {}\n'.format(self.atomic_species)
        string += '            Confidence: {}\n'.format(self.confidence)
        string += '    Graph parameters:\n'
        string += '        In-degree: {}\n'.format(self.in_degree)
        string += '        Out-degree: {}\n'.format(self.out_degree)
        string += '        Degree: {}\n'.format(self.degree)
        string += '        Alpha angles: ['
        for alpha in self.alpha_angles:
            string += ' {:.3f}'.format(alpha)
        string += ' ]\n'
        string += '            Max: {:.3f}\n'.format(self.alpha_max)
        string += '            Min: {:.3f}\n'.format(self.alpha_min)
        string += '        Theta angles: ['
        for theta in self.theta_angles:
            string += ' {:.3f}'.format(theta)
        string += ' ]\n'
        string += '            Max: {:.3f}\n'.format(self.theta_max)
        string += '            Min: {:.3f}\n'.format(self.theta_min)
        string += '            Variance: {:.3f}\n'.format(self.theta_angle_variance)
        string += '            Mean: {:.3f}\n'.format(self.theta_angle_mean)
        string += '        Normalized peak gamma: {:.3f}\n'.format(self.normalized_peak_gamma)
        string += '        Normalized average gamma: {:.3f}\n'.format(self.normalized_avg_gamma)
        string += '        Redshift: {:.3f}\n'.format(self.redshift)
        string += '        Average central separation: {}\n'.format(self.avg_central_separation)
        string += '    Local graph mapping:\n'
        string += '        District: {}\n'.format(self.district)
        string += '        In-neighbourhood: {}\n'.format(self.in_neighbourhood)
        string += '        Out-neighbourhood: {}\n'.format(self.out_neighbourhood)
        string += '        Neighbourhood: {}\n'.format(self.neighbourhood)
        string += '        Anti-Neighbourhood: {}\n'.format(self.anti_neighbourhood)
        string += '        Partners: {}\n'.format(self.partners)
        string += '        Semi-partners: {}\n'.format(self.semi_partners)
        string += '        Out-semi-partners: {}\n'.format(self.out_semi_partners)
        string += '        In-semi-partners: {}\n'.format(self.in_semi_partners)
        string += '    Local zeta-graph mapping:\n'
        string += '        Projected separation district: {}\n'.format(self.projected_separation_district)
        string += '        District: {}\n'.format(self.local_zeta_map['district'])
        string += '        In-neighbourhood: {}\n'.format(self.local_zeta_map['in_neighbourhood'])
        string += '        Out-neighbourhood: {}\n'.format(self.local_zeta_map['out_neighbourhood'])
        string += '        Neighbourhood: {}\n'.format(self.local_zeta_map['neighbourhood'])
        string += '        Partners: {}\n'.format(self.local_zeta_map['partners'])
        string += '        Semi-partners: {}\n'.format(self.local_zeta_map['semi_partners'])
        string += '        Out-semi-partners: {}\n'.format(self.local_zeta_map['out_semi_partners'])
        string += '        In-semi-partners: {}\n'.format(self.local_zeta_map['in_semi_partners'])
        string += '    Local anti-zeta-graph mapping:\n'
        if 'anti_in_neighbourhood' in self.local_zeta_map:
            string += '        In-neighbourhood: {}\n'.format(self.local_zeta_map['anti_in_neighbourhood'])
        if 'anti_out_neighbourhood' in self.local_zeta_map:
            string += '        Out-neighbourhood: {}\n'.format(self.local_zeta_map['anti_out_neighbourhood'])
        if 'anti_neighbourhood' in self.local_zeta_map:
            string += '        Neighbourhood: {}\n'.format(self.local_zeta_map['anti_neighbourhood'])
        if 'anti_partners' in self.local_zeta_map:
            string += '        Partners: {}\n'.format(self.local_zeta_map['anti_partners'])
        if 'anti_semi_partners' in self.local_zeta_map:
            string += '        Semi-partners: {}\n'.format(self.local_zeta_map['anti_semi_partners'])
        if 'anti_out_semi_partners' in self.local_zeta_map:
            string += '        Out-semi-partners: {}\n'.format(self.local_zeta_map['anti_out_semi_partners'])
        if 'anti_in_semi_partners' in self.local_zeta_map:
            string += '        In-semi-partners: {}\n'.format(self.local_zeta_map['anti_in_semi_partners'])
        string += '    Settings:\n'
        string += '        Is in precipitate: {}\n'.format(str(self.is_in_precipitate))
        string += '        Is edge column: {}\n'.format(str(self.is_edge_column))
        string += '        Is set by user: {}\n'.format(str(self.is_set_by_user))
        string += '        Show in overlay: {}\n'.format(str(self.show_in_overlay))
        string += '        Is void: {}\n'.format(str(self.void))
        string += '        Flag 1: {}\n'.format(str(self.flag_1))
        string += '        Flag 2: {}\n'.format(str(self.flag_2))
        string += '        Flag 3: {}\n'.format(str(self.flag_3))
        string += '        Flag 4: {}\n'.format(str(self.flag_4))
        string += '        Flag 5: {}\n'.format(str(self.flag_5))
        string += '        Flag 6: {}\n'.format(str(self.flag_6))
        string += '        Flag 7: {}\n'.format(str(self.flag_7))
        string += '        Flag 8: {}\n'.format(str(self.flag_8))
        string += '        Flag 9: {}\n'.format(str(self.flag_9))
        return string

    def im_pos(self):
        """Returns a triple with the image coordinates of the vertex.

        :return: (x, y, z), where x, y and z are the image coordinates of the vertex.
        :rtype: (float, float, float)
        """
        return self.im_coor_x, self.im_coor_y, self.im_coor_z

    def spatial_pos(self):
        """Returns a triple with the spatial coordinates of the vertex.

        :return: (x, y, z), where x, y and z are the spatial coordinates of the vertex.
        :rtype: (float, float, float)
        """
        return self.spatial_coor_x, self.spatial_coor_y, self.spatial_coor_z

    def set_zeta(self, zeta):
        if zeta == 0:
            self.zeta = 0
            self.im_coor_z = 0.0
            self.spatial_coor_z = zeta * 0.5 * self.al_lattice_const
        elif zeta == 1:
            self.zeta = 1
            self.im_coor_z = 1.0
            self.spatial_coor_z = zeta * 0.5 * self.al_lattice_const
        else:
            logger.warning('Zeta can only be 0 or 1. Zeta was not set for vertex {}!'.format(self.i))

    def set_probability_from_advanced(self):
        """Set atomic species probability vector by summarizing the probabilities of the every advanced species.

        If, for instance, the advanced probability vector is:

        =================== ========================
        advanced_species    advanced probability
        =================== ========================
        Si_1                0.01
        Si_2                0.2
        Al_1                0.09
        Al_2                0.7
        =================== ========================

        Then the probability vector will be set to:

        =================== ========================
        Atomic species      probability
        =================== ========================
        Si                  0.2 + 0.01 = 0.21
        Al                  0.09 + 0.7 = 0.79
        =================== ========================

        """
        self.probability_vector = {}
        for advanced_species, advanced_prob in self.advanced_probability_vector.items():
            atomic_species = self.parent_graph.species_dict['advanced_species'][advanced_species]['atomic_species']
            if atomic_species in self.probability_vector:
                self.probability_vector[atomic_species] += advanced_prob
            else:
                self.probability_vector[atomic_species] = advanced_prob
        self.probability_vector = utils.normalize_dict(self.probability_vector)

    def set_advanced_from_probability(self):
        pass

    def reset_probability_vector(self, bias='Un_1'):
        """Reset the probability vector with a certain bias

        This method will set all the atomic species probabilities equal to each other. If an advanced species bias is
        given, then that advanced species will be slightly larger than the others. The bias defaults to 'Un_1' if no
        bias is given.


        :param bias: Default to this advanced species (optional, default='Un_1')
        :type bias: string

        """
        self.advanced_probability_vector = {}
        self.probability_vector = {}
        if self.parent_graph is not None:
            for key, value in self.parent_graph.species_dict['advanced_species'].items():
                if key == bias:
                    self.advanced_probability_vector[key] = 1.1
                    if value['atomic_species'] in self.probability_vector:
                        self.probability_vector[value['atomic_species']] += 1.1
                    else:
                        self.probability_vector[value['atomic_species']] = 1.1
                else:
                    self.advanced_probability_vector[key] = 1.0
                    if value['atomic_species'] in self.probability_vector:
                        self.probability_vector[value['atomic_species']] += 1.0
                    else:
                        self.probability_vector[value['atomic_species']] = 1.0
            self.advanced_probability_vector = utils.normalize_dict(self.advanced_probability_vector)
            self.probability_vector = utils.normalize_dict(self.probability_vector)
            self.atomic_species = self.parent_graph.species_dict['advanced_species'][bias]['atomic_species']
            self.advanced_species = bias

    def determine_species_from_probability_vector(self):
        """Set the atomic species of the vertex based on the advanced probability vector

        First determine the advanced species from the advanced probability vector, and then set the atomic species
        of the vertex based on which atomic species the advanced species maps to. Also set the symmetry number 'n'.

        """
        self.advanced_species = max(self.advanced_probability_vector, key=self.advanced_probability_vector.get)
        self.atomic_species = self.parent_graph.species_dict['advanced_species'][self.advanced_species]['atomic_species']
        self.set_probability_from_advanced()
        self.n = self.parent_graph.species_dict['advanced_species'][self.advanced_species]['n']

    def reset_probability_vector_from_atomic_species(self):
        self.reset_probability_vector(bias=self.advanced_species)

    def op_arc_pivot(self, j, k):
        if j == k or j == self.i or k == self.i:
            return False

        pos_j = -1
        pos_k = -1
        if j in self.district:
            pos_j = self.district.index(j)
        if k in self.district:
            pos_k = self.district.index(k)

        if pos_j == -1 and pos_k == -1:
            self.district[-1] = k
            return True
        elif not pos_j == -1 and not pos_k == -1:
            self.district[pos_j], self.district[pos_k] = self.district[pos_k], self.district[pos_j]
            return True
        elif pos_j == -1:
            return False
        else:
            self.district[-1] = k
            self.district[pos_j], self.district[pos_k] = self.district[pos_k], self.district[pos_j]
            return True

    def permute_zeta_j_k(self, j, k):
        if j == k:
            return False
        pos_j = -1
        pos_k = -1
        if j in self.zeta_district:
            pos_j = self.zeta_district.index(j)
        if k in self.zeta_district:
            pos_k = self.zeta_district.index(k)

        if pos_j == -1 and pos_k == -1:
            self.zeta_district[-1] = k
            return True
        elif not pos_j == -1 and not pos_k == -1:
            self.zeta_district[pos_j], self.zeta_district[pos_k] = self.zeta_district[pos_k], self.zeta_district[pos_j]
            return True
        elif pos_j == -1:
            return False
        else:
            self.zeta_district[-1] = k
            self.zeta_district[pos_j], self.zeta_district[pos_k] = self.zeta_district[pos_k], self.zeta_district[pos_j]
            return True

    def permute_pos_j_pos_k(self, pos_j, pos_k):
        self.district[pos_j], self.district[pos_k] = self.district[pos_k], self.district[pos_j]
        return True

    def shift_pos_j_pos_k(self, pos_j, pos_k):
        if pos_k == len(self.district) - 1 or pos_k == -1:
            new_district = self.district[:pos_j - 1] + self.district[pos_j + 1:pos_k] +\
                           [self.district[pos_j]]
        else:
            new_district = self.district[:pos_j - 1] + self.district[pos_j + 1:pos_k] + \
                           [self.district[pos_j]] + self.district[pos_k + 1:-1]
        self.district = new_district
        return True

    def put_j_last(self, j):
        if j in self.district:
            pos_j = self.district.index(j)
            if pos_j < self.n - 1:
                self.permute_pos_j_pos_k(pos_j, self.n - 1)
            elif pos_j == self.n - 1:
                pass
            else:
                pass
        else:
            pass

    def put_j_above_last(self, j):
        if j in self.district:
            pos_j = self.district.index(j)
            if pos_j > self.n:
                self.permute_pos_j_pos_k(pos_j, self.n)
            elif pos_j == self.n:
                pass
            else:
                pass
        else:
            pass

    def partner_query(self, j):
        if j in self.district[:self.n]:
            return True
        else:
            return False

    def anti_zeta(self):
        if self.zeta == 0:
            return 1
        else:
            return 0

    def decrement_n(self):
        for advanced_species, value in self.parent_graph.species_dict['advanced_species'].items():
            if value['n'] == self.n - 1 and not value['atomic_species'] == 'Un':
                self.reset_probability_vector(bias=advanced_species)
                self.determine_species_from_probability_vector()
                break
        else:
            return False
        return True

    def increment_n(self):
        for advanced_species, value in self.parent_graph.species_dict['advanced_species'].items():
            if value['n'] == self.n + 1:
                self.reset_probability_vector(bias=advanced_species)
                self.determine_species_from_probability_vector()
                break
        else:
            return False
        return True


class Arc:

    al_lattice_const = 404.95

    def __init__(self, j, vertex_a, vertex_b, parent_graph=None):

        self.j = j
        self.vertex_a = vertex_a
        self.vertex_b = vertex_b

        self.parent = parent_graph

        if vertex_a.i in vertex_b.out_neighbourhood:
            self.dual_arc = True
        else:
            self.dual_arc = False

        if vertex_a.zeta == vertex_b.zeta:
            self.co_planar = True
        else:
            self.co_planar = False

        self.im_separation = 0
        self.im_projected_separation = 0
        self.spatial_separation = 0
        self.spatial_projected_separation = 0
        self.hard_sphere_separation = 0
        self.redshift = 0

        self.calc_properties()

    def calc_properties(self):
        pos_i = self.vertex_a.im_pos()
        pos_j = self.vertex_b.im_pos()
        self.im_separation = np.sqrt((pos_i[0] - pos_j[0]) ** 2 + (pos_i[1] - pos_j[1]) ** 2 + (pos_i[2] - pos_j[2]) ** 2)
        self.im_projected_separation = np.sqrt((pos_i[0] - pos_j[0]) ** 2 + (pos_i[1] - pos_j[1]) ** 2)
        pos_i = self.vertex_a.spatial_pos()
        pos_j = self.vertex_b.spatial_pos()
        self.spatial_separation = np.sqrt((pos_i[0] - pos_j[0]) ** 2 + (pos_i[1] - pos_j[1]) ** 2 + (pos_i[2] - pos_j[2]) ** 2)
        self.spatial_projected_separation = np.sqrt((pos_i[0] - pos_j[0]) ** 2 + (pos_i[1] - pos_j[1]) ** 2)
        if self.parent:
            radii_1 = self.parent.species_dict[self.vertex_a.advanced_species]['atomic_radii']
            radii_2 = self.parent.species_dict[self.vertex_b.advanced_species]['atomic_radii']
        else:
            radii_1 = 100
            radii_2 = 100
        self.hard_sphere_separation = radii_1 + radii_2
        self.redshift = self.hard_sphere_separation - self.spatial_separation


class AtomicGraph:

    al_lattice_const = 404.95

    def __init__(self, scale, species_dict, active_model=None, district_size=10):

        self.species_dict = species_dict

        # Contents:
        self.vertices = []
        self.arcs = []
        self.arc_indices = []
        self.anti_arcs = []
        self.anti_arc_indices = []
        self.meshes = []
        self.mesh_indices = []

        self.particle_boarder_indices = []
        self.structures = []

        self.adjacency_matrix = None
        self.projected_separation_matrix = None

        self.scale = scale
        self.district_size = district_size

        if active_model is None:
            self.active_model = 'default_model'
        else:
            self.active_model = active_model

        # Stats
        self.chi = 0
        self.zeta_chi = 0
        self.order = 0
        self.size = 0
        self.volume = 0
        self.anti_size = 0
        self.avg_degree = 0
        self.avg_zeta_degree = 0
        self.matrix_redshift = 0
        self.particle_redshift = 0
        self.total_redshift = 0
        self.precipitate_displacement = 0

    def __str__(self):
        return self.report()

    def report(self):
        string = '    Order: {}\n'.format(self.order)
        string += '    Size: {}\n'.format(self.size)
        string += '    Chi: {:.3f}\n'.format(self.chi)
        string += '    Average degree: {:.3f}\n'.format(self.avg_degree)
        string += '    Matrix redshift: {:.3f}\n'.format(self.matrix_redshift)
        string += '    Particle redshift: {:.3f}\n'.format(self.particle_redshift)
        string += '    Total redshift: {:.3f}\n'.format(self.total_redshift)
        return string

    def vertex_report(self, i):
        string = self.vertices[i].report()
        return string

    def add_vertex(self, new_vertex):
        self.vertices.append(new_vertex)
        if not new_vertex.i == len(self.vertices) - 1:
            new_vertex.i = len(self.vertices) - 1
        i = new_vertex.i
        self.order += 1
        if self.projected_separation_matrix is not None:
            vertical_column = np.zeros([self.order - 1, 1], dtype=float)
            for j in range(0, self.order - 1):
                vertical_column[j, 0] = self.get_projected_separation(i, j)
            self.projected_separation_matrix = np.concatenate((self.projected_separation_matrix, vertical_column), axis=1)
            horizontal_column = np.zeros([1, self.order], dtype=float)
            for m, mirror in enumerate(vertical_column):
                horizontal_column[0, m] = mirror
            horizontal_column[0, i] = 0.0
            self.projected_separation_matrix = np.concatenate((self.projected_separation_matrix, horizontal_column), axis=0)
            column_characterization.determine_districts(self, [i])
            county = set()
            for citizen in self.vertices[i].district:
                county.add(citizen)
                for extended_citizen in self.vertices[citizen].district:
                    county.add(extended_citizen)
            county = list(county)
            column_characterization.determine_districts(self, county)
            # Remap the vertex maps
            self.build_local_maps(build_out=True)
            self.build_local_zeta_maps(build_out=True)
        if self.adjacency_matrix is None:
            self.adjacency_matrix = np.zeros([self.order, self.order], dtype=int)
        else:
            self.adjacency_matrix = np.concatenate((self.adjacency_matrix, np.zeros([self.order - 1, 1], dtype=int)), axis=1)
            self.adjacency_matrix = np.concatenate((self.adjacency_matrix, np.zeros([1, self.order], dtype=int)), axis=0)
        self.summarize_stats()

    def vertex_moved(self, i):
        if self.projected_separation_matrix is not None:
            for j in range(0, self.order):
                separation = self.get_projected_separation(i, j)
                self.projected_separation_matrix[i, j] = separation
                self.projected_separation_matrix[j, i] = separation
            column_characterization.determine_districts(self, [i])
            county = set()
            for citizen in self.vertices[i].district:
                county.add(citizen)
                for extended_citizen in self.vertices[citizen].district:
                    county.add(extended_citizen)
            for vertex in self.vertices:
                if i in vertex.district:
                    county.add(vertex.i)
            county = list(county)
            column_characterization.determine_districts(self, county)
            # Remap the vertex maps
            self.build_local_maps(build_out=True)
            self.build_local_zeta_maps(build_out=True)
        self.summarize_stats()

    def remove_vertex(self, i):
        logger.info('Removing vertex {}'.format(i))
        # Delete
        del self.vertices[i]
        county = set()
        for k, vertex in enumerate(self.vertices):
            vertex.i = k
            for c, citizen in enumerate(vertex.district):
                if citizen > i:
                    vertex.district[c] -= 1
                elif citizen == i:
                    county.add(vertex.i)
        self.order -= 1
        # Remove entry from distance matrix:
        if self.projected_separation_matrix is not None:
            self.projected_separation_matrix = np.delete(self.projected_separation_matrix, i, axis=0)
            self.projected_separation_matrix = np.delete(self.projected_separation_matrix, i, axis=1)
            # Repair districts:
            county = list(county)
            column_characterization.determine_districts(self, county)
            # Remap the vertex maps
            logger.info('    Rebuilding vertex maps...')
            self.build_local_maps(build_out=True)
            self.build_local_zeta_maps(build_out=True)
        if self.adjacency_matrix is not None:
            self.adjacency_matrix = np.delete(self.adjacency_matrix, i, axis=0)
            self.adjacency_matrix = np.delete(self.adjacency_matrix, i, axis=1)
        self.summarize_stats()

    def get_classes(self):
        classes = set()
        for item in self.species_dict:
            classes.add(item)
        return classes

    def get_arc(self, i, j):
        """Return the arc object adjacent from i to j, if any.

        :param i: Index of the tail vertex
        :param j: Index of the head vertex
        :type i: int
        :type j: int

        :returns Arc instance:
        :rtype graph_2.Arc:

        """
        result = None
        for arc in self.arcs:
            if arc.vertex_a.i == i and arc.vertex_b.i == j:
                result = arc
                break
        else:
            if j in self.vertices[i].out_neighbourhood:
                result = Arc(-1, self.vertices[i], self.vertices[j])
        return result

    def get_vertex_objects_from_indices(self, vertex_indices=None):
        """Get a list of Vertex objects corresponding to the indices given

        :param vertex_indices: List of indices
        :type vertex_indices: list(int)

        :returns List of vertex objects:
        :rtype list(graph_2.Vertex):

        """
        vertices = []
        if vertex_indices is not None:
            for index in vertex_indices:
                vertices.append(self.vertices[index])
        return vertices

    def get_alpha_angles(self, i, selection_type='zeta'):
        pivot = (self.vertices[i].im_coor_x, self.vertices[i].im_coor_y)
        district = self.vertices[i].district
        if len(district) == 0:
            return []
        out_neighbourhood = self.vertices[i].out_neighbourhood
        partners = self.vertices[i].partners

        if selection_type == 'partners':
            j = []
            for citizen in district:
                if citizen in partners:
                    j.append(citizen)
                if len(j) == 3:
                    break
            else:
                for citizen in district:
                    if citizen not in j:
                        j.append(citizen)
                    if len(j) == 3:
                        break
            j.append(j[0])
            for k, index in enumerate(j):
                j[k] = (self.vertices[index].im_coor_x, self.vertices[index].im_coor_y)

        elif selection_type == 'zeta':
            j = []
            for citizen in district:
                if not self.vertices[citizen].zeta == self.vertices[i].zeta:
                    j.append(citizen)
                if len(j) == 3:
                    break
            else:
                for citizen in district:
                    if citizen not in j:
                        j.append(citizen)
                    if len(j) == 3:
                        break
            j.append(j[0])
            for k, index in enumerate(j):
                j[k] = (self.vertices[index].im_coor_x, self.vertices[index].im_coor_y)

        elif selection_type == 'out':
            j = []
            for citizen in district:
                if citizen in out_neighbourhood:
                    j.append(citizen)
                if len(j) == 3:
                    break
            else:
                for citizen in district:
                    if citizen not in j:
                        j.append(citizen)
                    if len(j) == 3:
                        break
            j.append(j[0])
            for k, index in enumerate(j):
                j[k] = (self.vertices[index].im_coor_x, self.vertices[index].im_coor_y)

        else:
            j_1 = (self.vertices[district[0]].im_coor_x, self.vertices[district[0]].im_coor_y)
            j_2 = (self.vertices[district[1]].im_coor_x, self.vertices[district[1]].im_coor_y)
            j_3 = (self.vertices[district[2]].im_coor_x, self.vertices[district[2]].im_coor_y)
            j = [j_1, j_2, j_3, j_1]

        alpha = []
        for k in range(0, 3):
            alpha.append(utils.find_angle_from_points(j[k], j[k + 1], pivot))

        if sum(alpha) > 6.5:
            for x in range(0, 3):
                alpha[x] = 2 * np.pi - alpha[x]

        return alpha

    def get_zeta_alpha_angles(self, i):
        pivot = (self.vertices[i].im_coor_x, self.vertices[i].im_coor_y)
        district = self.vertices[i].local_zeta_map['district']
        if len(district) == 0:
            return []

        j = []
        for citizen in district:
            if not self.vertices[citizen].zeta == self.vertices[i].zeta:
                j.append(citizen)
            if len(j) == 3:
                break
        else:
            for citizen in district:
                if citizen not in j:
                    j.append(citizen)
                if len(j) == 3:
                    break
        j.append(j[0])
        for k, index in enumerate(j):
            j[k] = (self.vertices[index].im_coor_x, self.vertices[index].im_coor_y)

        alpha = []
        for k in range(0, 3):
            alpha.append(utils.find_angle_from_points(j[k], j[k + 1], pivot))

        if sum(alpha) > 6.5:
            for x in range(0, 3):
                alpha[x] = 2 * np.pi - alpha[x]

        return alpha

    def get_anti_zeta_alpha_angles(self, i):
        pivot = (self.vertices[i].im_coor_x, self.vertices[i].im_coor_y)
        district = self.vertices[i].local_zeta_map['district']
        if len(district) == 0:
            return []

        j = []
        for citizen in district:
            if self.vertices[citizen].zeta == self.vertices[i].zeta:
                j.append(citizen)
            if len(j) == 3:
                break
        else:
            for citizen in district:
                if citizen not in j:
                    j.append(citizen)
                if len(j) == 3:
                    break
        j.append(j[0])
        for k, index in enumerate(j):
            j[k] = (self.vertices[index].im_coor_x, self.vertices[index].im_coor_y)

        alpha = []
        for k in range(0, 3):
            alpha.append(utils.find_angle_from_points(j[k], j[k + 1], pivot))

        if sum(alpha) > 6.5:
            for x in range(0, 3):
                alpha[x] = 2 * np.pi - alpha[x]

        return alpha

    def get_theta_angles(self, i, selection_type='normal'):
        if selection_type == 'normal':
            sub_graph = self.get_column_centered_subgraph(i)
            theta = []
            for mesh in sub_graph.meshes:
                theta.append(mesh.angles[0])
        else:
            sub_graph = self.get_column_centered_subgraph(i)
            theta = []
            for mesh in sub_graph.meshes:
                if mesh.order == 4:
                    theta.append(mesh.angles[0])
        return theta

    def get_zeta_theta_angles(self, i):
        pass

    def get_anti_zeta_theta_angles(self, i):
        pass

    def get_redshift(self, i, j):
        hard_sphere_separation = self.get_hard_sphere_separation(i, j)
        actual_separation = self.get_separation(i, j)
        return hard_sphere_separation - actual_separation

    def get_image_separation(self, i, j):
        pos_i = self.vertices[i].im_pos()
        pos_j = self.vertices[j].im_pos()
        separation = np.sqrt((pos_i[0] - pos_j[0]) ** 2 + (pos_i[1] - pos_j[1]) ** 2 + (pos_i[2] - pos_j[2]) ** 2)
        return separation

    def get_projected_image_separation(self, i, j):
        pos_i = self.vertices[i].im_pos()
        pos_j = self.vertices[j].im_pos()
        projected_separation = np.sqrt((pos_i[0] - pos_j[0]) ** 2 + (pos_i[1] - pos_j[1]) ** 2)
        return projected_separation

    def get_separation(self, i, j):
        pos_i = self.vertices[i].spatial_pos()
        pos_j = self.vertices[j].spatial_pos()
        separation = np.sqrt((pos_i[0] - pos_j[0]) ** 2 + (pos_i[1] - pos_j[1]) ** 2 + (pos_i[2] - pos_j[2]) ** 2)
        return separation

    def get_projected_separation(self, i, j):
        pos_i = self.vertices[i].spatial_pos()
        pos_j = self.vertices[j].spatial_pos()
        projected_separation = np.sqrt((pos_i[0] - pos_j[0]) ** 2 + (pos_i[1] - pos_j[1]) ** 2)
        return projected_separation

    def get_hard_sphere_separation(self, i, j):
        radii_1 = self.species_dict['atomic_species'][self.vertices[i].atomic_species]['atomic_radii']
        radii_2 = self.species_dict['atomic_species'][self.vertices[j].atomic_species]['atomic_radii']
        return radii_1 + radii_2

    def get_im_separation_from_spatial_separation(self, spatial_searation):
        return spatial_searation / self.scale

    def get_spatial_separation_from_im_separation(self, im_separation):
        return im_separation * self.scale

    def get_adjacency_matrix(self):
        self.summarize_stats()
        if self.adjacency_matrix is None:
            self.adjacency_matrix = np.zeros([self.order, self.order], dtype=int)
        for vertex in self.vertices:
            for out_neighbour in vertex.out_neighbourhood:
                self.adjacency_matrix[vertex.i, out_neighbour] = 1
        return self.adjacency_matrix

    def get_anti_graph(self):
        logger.info('Building anti-graph..')
        anti_graph = AntiGraph(self)
        logger.info('Anti-graph built.')
        return anti_graph

    def get_mesh(self, i, j, mesh_index=0):
        indices, angles, vectors = self.get_mesh_numerics(i, j)
        mesh = Mesh(mesh_index, self.get_vertex_objects_from_indices(indices))
        mesh.angles = angles
        mesh.angle_vectors = vectors
        return mesh

    def get_induced_subgraph(self, vertex_indices):
        pass

    def get_mesh_centered_subgraph(self, i, j, order=1):
        sub_graph = SubGraph()
        mesh = self.get_mesh(i, j)
        sub_graph.add_mesh(mesh)
        sub_graph.finalize_init()
        return sub_graph

    def get_arc_centered_subgraph(self, i, j, order=1):
        mesh_0 = self.get_mesh(i, j, 0)
        mesh_1 = self.get_mesh(j, i, 1)
        sub_graph = SubGraph()
        sub_graph.add_mesh(mesh_0)
        sub_graph.add_mesh(mesh_1)
        sub_graph.finalize_init()
        return sub_graph

    def get_column_centered_subgraph(self, i, order=1):
        sub_graph = SubGraph()
        for mesh_index, neighbour in enumerate(self.vertices[i].neighbourhood):
            mesh = self.get_mesh(i, neighbour, mesh_index)
            sub_graph.add_mesh(mesh)
        sub_graph.finalize_init()
        return sub_graph

    def get_mesh_numerics(self, i, j, order=0):
        corners = [i, j]
        counter = 0
        backup_counter = 0
        stop = False
        while not stop:
            angle, next_index = self.angle_sort(i, j)
            if next_index == corners[0] or counter > 14:
                _, nextnext = self.angle_sort(j, next_index)
                if not nextnext == corners[1]:
                    corners, i, j = self.rebase(corners, nextnext, next_index, append=False)
                stop = True
            elif next_index in corners:
                corners, i, j = self.rebase(corners, next_index, j)
                counter = len(corners) - 2
            else:
                corners.append(next_index)
                counter += 1
                i, j = j, next_index
            backup_counter += 1
            if backup_counter > 25 or len(corners) == 1:
                # logger.warning('Emergency stop!')
                stop = True
        angles = []
        vectors = []
        for m, corner in enumerate(corners):
            pivot = self.vertices[corner].im_pos()
            if m == 0:
                p1 = self.vertices[corners[len(corners) - 1]].im_pos()
                p2 = self.vertices[corners[m + 1]].im_pos()
            elif m == len(corners) - 1:
                p1 = self.vertices[corners[m - 1]].im_pos()
                p2 = self.vertices[corners[0]].im_pos()
            else:
                p1 = self.vertices[corners[m - 1]].im_pos()
                p2 = self.vertices[corners[m + 1]].im_pos()
            angle = utils.find_angle_from_points(p1[:2], p2[:2], pivot[:2])
            theta = angle / 2
            vector = (p1[0] - pivot[0], p1[1] - pivot[1])
            length = utils.vector_magnitude(vector)
            vector = (vector[0] / length, vector[1] / length)
            vector = (vector[0] * np.cos(theta) + vector[1] * np.sin(theta),
                      -vector[0] * np.sin(theta) + vector[1] * np.cos(theta))
            angles.append(angle)
            vectors.append(vector)
        return corners, angles, vectors

    def get_advanced_species_list(self, species_dict=None):
        advanced_species_list = []
        if species_dict is None:
            for key in self.species_dict['advanced_species']:
                advanced_species_list.append(key)
        else:
            for key in species_dict['advanced_species']:
                advanced_species_list.append(key)
        return advanced_species_list

    def get_atomic_species_list(self, species_dict=None):
        atomic_species_list = []
        if species_dict is None:
            for key in self.species_dict['atomic_species']:
                atomic_species_list.append(key)
        else:
            for key in species_dict['atomic_species']:
                atomic_species_list.append(key)
        return atomic_species_list

    def set_zeta(self, i, zeta):

        """Set the zeta value of the vertex with index i.

        :param i: vertex index
        :param zeta: Zeta
        :type i: int
        :type zeta: bool

        """

        if zeta is True:
            zeta = 1
        elif zeta is False:
            zeta = 0
        self.vertices[i].set_zeta(zeta)
        self.build_local_zeta_map([i] + self.vertices[i].district)

    def set_n(self, i, n):
        for key, item in self.species_dict['advanced_species'].items():
            if item['n'] == n and not item['atomic_species'] == 'Un':
                self.set_advanced_species(i, key)
                break
        else:
            logger.info('Could not set the symmetry of vertex {}, n = {} is not present in the species dictionary!'.format(i, n))

    def set_advanced_species(self, i, advanced_species):
        if not self.vertices[i].advanced_species == advanced_species:
            if advanced_species in self.species_dict['advanced_species']:
                self.vertices[i].reset_probability_vector(bias=advanced_species)
                self.vertices[i].determine_species_from_probability_vector()
                self.build_local_map([i] + self.vertices[i].district)
                self.build_local_zeta_map([i] + self.vertices[i].district)
            else:
                self.vertices[i].reset_probability_vector()
                self.vertices[i].determine_species_from_probability_vector()
                self.build_local_map([i] + self.vertices[i].district)
                self.build_local_zeta_map([i] + self.vertices[i].district)

    def set_atomic_species(self, i, atomic_species):
        if not self.vertices[i].atomic_species == atomic_species:
            if atomic_species in self.species_dict['atomic_species']:
                for advanced_species, details in self.species_dict['advanced_species'].items():
                    if details['atomic_species'] == atomic_species:
                        self.vertices[i].reset_probability_vector(bias=advanced_species)
                        self.vertices[i].determine_species_from_probability_vector()
                        self.build_local_map([i] + self.vertices[i].district)
                        self.build_local_zeta_map([i] + self.vertices[i].district)
                        break
                else:
                    self.vertices[i].reset_probability_vector()
                    self.vertices[i].determine_species_from_probability_vector()
                    self.build_local_map([i] + self.vertices[i].district)
                    self.build_local_zeta_map([i] + self.vertices[i].district)
            else:
                self.vertices[i].reset_probability_vector()
                self.vertices[i].determine_species_from_probability_vector()
                self.build_local_map([i] + self.vertices[i].district)
                self.build_local_zeta_map([i] + self.vertices[i].district)

    @staticmethod
    def rebase(corners, next_, j, append=True):
        for k, corner in enumerate(corners):
            if corner == next_:
                del corners[k + 1:]
                if append:
                    corners.append(j)
                break
        return corners, next_, j

    def angle_sort(self, i, j):
        min_angle = 1000
        next_index = -1
        p1 = self.vertices[i].im_pos()
        pivot = self.vertices[j].im_pos()
        search_list = self.vertices[j].neighbourhood
        for k in search_list:
            if not k == i:
                p2 = self.vertices[k].im_pos()
                alpha = utils.find_angle_from_points(p1[:2], p2[:2], pivot[:2])
                if alpha < min_angle:
                    min_angle = alpha
                    next_index = k
        logger.debug('Found next: {}'.format(next_index))
        return min_angle, next_index

    def invert_levels(self):
        for vertex in self.vertices:
            if vertex.level == 0:
                vertex.level = 1
            else:
                vertex.level = 0

    def reset_all_flags(self):
        for vertex in self.vertices:
            vertex.flag_1 = False
            vertex.flag_2 = False
            vertex.flag_3 = False
            vertex.flag_4 = False
            vertex.flag_5 = False
            vertex.flag_6 = False
            vertex.flag_7 = False
            vertex.flag_8 = False
            vertex.flag_9 = False
        logger.info('All flags reset!')

    def determine_zeta_districts(self):
        for vertex in self.vertices:
            vertex.zeta_district = []
            for citizen in vertex.district:
                if self.vertices[citizen].zeta == vertex.anti_zeta():
                    vertex.zeta_district.append(citizen)
            if 'district' not in vertex.local_zeta_map:
                vertex.local_zeta_map['district'] = vertex.projected_separation_district

    def build_local_maps(self, build_out=True):
        self.build_local_map([vertex.i for vertex in self.vertices], build_out=build_out)

    def build_local_map(self, indices, build_out=True):
        # Determine out_neighbourhoods:
        if build_out:
            for vertex in self.get_vertex_objects_from_indices(indices):
                vertex.out_neighbourhood = set()
                counter = 0
                for citizen in vertex.district:
                    vertex.out_neighbourhood.add(citizen)
                    counter += 1
                    if counter == vertex.n:
                        break

        # Determine in_neighbourhoods:
        for vertex in self.get_vertex_objects_from_indices(indices):
            vertex.in_neighbourhood = set()
            if not vertex.void:
                for candidate in self.vertices:
                    if vertex.i in candidate.out_neighbourhood:
                        vertex.in_neighbourhood.add(candidate.i)
        # Determine neighbourhood:
        for vertex in self.get_vertex_objects_from_indices(indices):
            vertex.neighbourhood = set()
            if not vertex.void:
                vertex.neighbourhood = vertex.out_neighbourhood.union(vertex.in_neighbourhood)
        # Determine anti_neighbourhood:
        for vertex in self.get_vertex_objects_from_indices(indices):
            vertex.anti_neighbourhood = set()
            if not vertex.void:
                for citizen in vertex.district:
                    if citizen not in vertex.neighbourhood:
                        vertex.anti_neighbourhood.add(citizen)
        # Determine partners and semi-partners:
        for vertex in self.get_vertex_objects_from_indices(indices):
            vertex.partners = set()
            vertex.semi_partners = set()
            vertex.in_semi_partners = set()
            vertex.out_semi_partners = set()
            if not vertex.void:
                vertex.partners = vertex.out_neighbourhood.intersection(vertex.in_neighbourhood)
                vertex.semi_partners = vertex.neighbourhood - vertex.partners
                vertex.out_semi_partners = vertex.semi_partners.intersection(vertex.out_neighbourhood)
                vertex.in_semi_partners = vertex.semi_partners.intersection(vertex.in_neighbourhood)

                vertex.in_degree = len(vertex.in_neighbourhood)
                vertex.out_degree = len(vertex.out_neighbourhood)
                vertex.degree = len(vertex.neighbourhood)

    def build_local_zeta_maps(self, build_out=True):
        self.build_local_zeta_map([vertex.i for vertex in self.vertices], build_out=build_out)

    def build_local_zeta_map(self, indices, build_out=True):
        # Rebuild zeta districts:
        self.determine_zeta_districts()
        # Build out-neighbourhoods
        if build_out:
            for vertex in self.get_vertex_objects_from_indices(indices):
                vertex.local_zeta_map['out_neighbourhood'] = set()
                counter = 0
                for citizen in vertex.local_zeta_map['district']:
                    if self.vertices[citizen].zeta == vertex.anti_zeta():
                        vertex.local_zeta_map['out_neighbourhood'].add(citizen)
                        counter += 1
                    if counter == vertex.n:
                        break
                else:
                    for citizen in vertex.local_zeta_map['district']:
                        if citizen not in vertex.local_zeta_map['out_neighbourhood']:
                            vertex.local_zeta_map['out_neighbourhood'].add(citizen)
                            counter += 1
                        if counter == vertex.n:
                            break
        # Determine in_neighbourhoods:
        for vertex in self.get_vertex_objects_from_indices(indices):
            vertex.local_zeta_map['in_neighbourhood'] = set()
            for candidate in self.vertices:
                if vertex.i in candidate.local_zeta_map['out_neighbourhood']:
                    vertex.local_zeta_map['in_neighbourhood'].add(candidate.i)
        # Determine neighbourhood:
        for vertex in self.get_vertex_objects_from_indices(indices):
            vertex.local_zeta_map['neighbourhood'] = set()
            vertex.local_zeta_map['neighbourhood'] = vertex.local_zeta_map['out_neighbourhood'].union(vertex.local_zeta_map['in_neighbourhood'])
        # Determine partners and semi-partners:
        for vertex in self.get_vertex_objects_from_indices(indices):
            vertex.local_zeta_map['partners'] = set()
            vertex.local_zeta_map['semi_partners'] = set()
            vertex.local_zeta_map['in_semi_partners'] = set()
            vertex.local_zeta_map['out_semi_partners'] = set()
            vertex.local_zeta_map['partners'] = vertex.local_zeta_map['out_neighbourhood'].intersection(vertex.local_zeta_map['in_neighbourhood'])
            vertex.local_zeta_map['semi_partners'] = vertex.local_zeta_map['neighbourhood'] - vertex.local_zeta_map['partners']
            vertex.local_zeta_map['out_semi_partners'] = vertex.local_zeta_map['semi_partners'].intersection(vertex.local_zeta_map['out_neighbourhood'])
            vertex.local_zeta_map['in_semi_partners'] = vertex.local_zeta_map['semi_partners'].intersection(vertex.local_zeta_map['in_neighbourhood'])

    def build_cities(self):
        for vertex in self.vertices:
            vertex.city = set()
            vertex.city.add(vertex.i)
            for citizen in vertex.district:
                vertex.city.add(citizen)
                for extended_citizen in self.vertices[citizen].district:
                    vertex.city.add(extended_citizen)

    def match_zeta_graph(self):
        for vertex in self.vertices:
            reordered_district = []
            for zeta_citizen in vertex.local_zeta_map['district']:
                if vertex.zeta == self.vertices[zeta_citizen].anti_zeta():
                    reordered_district.append(zeta_citizen)
            for zeta_citizen in vertex.local_zeta_map['district']:
                if zeta_citizen not in reordered_district:
                    reordered_district.append(zeta_citizen)
            vertex.district = reordered_district
        self.build_local_maps()

    def op_arc_pivot(self, i, j, k):
        if j in self.vertices[i].out_neighbourhood:
            if k in self.vertices[i].out_neighbourhood:
                return False
            else:
                if self.vertices[i].op_arc_pivot(j, k):
                    self.vertices[i].out_neighbourhood.discard(j)
                    self.vertices[i].out_neighbourhood.add(k)
                    self.build_local_map([i] + self.vertices[i].district, build_out=False)
                    return True
                else:
                    return False
        else:
            return False

    def zeta_op_arc_pivot(self, i, j, k):
        if self.vertices[i].zeta_op_arc_pivot(j, k):
            self.vertices[i].local_zeta_map['out_neighbourhood'].discard(j)
            self.vertices[i].local_zeta_map['out_neighbourhood'].add(k)
            self.build_local_zeta_map([i] + self.vertices[i].zeta_district, build_out=False)

    def op_weak_arc_termination(self, i, j, aggressive=False):

        config = self.get_column_centered_subgraph(i)
        options = []

        for mesh in config.meshes:
            for m, corner in enumerate(mesh.vertex_indices):
                if m not in [0, 1, len(mesh.vertex_indices) - 1]:
                    options.append(corner)
                if m == 1 and corner not in self.vertices[i].out_neighbourhood:
                    options.append(corner)

        for option in options:

            mesh_1 = self.get_mesh(i, option)
            mesh_2 = self.get_mesh(option, i)
            if mesh_1.order == 4 and mesh_2.order == 4:
                k = option
                break

        else:

            if aggressive:
                for option in options:
                    mesh_1 = self.get_mesh(i, option)
                    mesh_2 = self.get_mesh(option, i)
                    if mesh_1.order == 4 or mesh_2.order == 4:
                        k = option
                        break

                else:
                    return -1

            else:
                return -1

        return k

    def op_weak_arc_preservation(self, i, j):

        config = self.get_column_centered_subgraph(j)
        options = []

        for mesh in config.meshes:
            for m, corner in enumerate(mesh.vertex_indices):
                if (m == 1 or m == len(mesh.vertex_indices) - 1) and corner in self.vertices[j].out_neighbourhood and not corner == i:
                    options.append(corner)

        k = -1
        for option in options:
            mesh_1 = self.get_mesh(j, option)
            mesh_2 = self.get_mesh(option, j)
            if mesh_1.order == 3 and mesh_2.order == 3:
                k = option
                break

        return k

    def op_strong_arc_termination(self, i, j):
        self.vertices[i].put_j_last(j)
        if self.vertices[i].decrement_n():
            self.build_local_map([i] + self.vertices[i].district, build_out=True)
            return True

    def op_strong_arc_preservation(self, i, j):
        self.vertices[j].put_j_above_last(i)
        if self.vertices[j].increment_n():
            self.build_local_map([j] + self.vertices[j].district, build_out=True)
            return True

    def find_intersections(self):

        intersecting_segments = []

        for a in self.vertices:
            a_coor = a.im_pos()
            a_coor = (a_coor[0], a_coor[1])
            for b in [self.vertices[index] for index in a.out_neighbourhood]:
                if not a.is_edge_column and not b.is_edge_column:
                    b_coor = b.im_pos()
                    b_coor = (b_coor[0], b_coor[1])
                    for c in [self.vertices[index] for index in a.out_neighbourhood]:
                        if not c.i == b.i:
                            c_coor = c.im_pos()
                            c_coor = (c_coor[0], c_coor[1])
                            for d in [self.vertices[index] for index in c.out_neighbourhood]:
                                d_coor = d.im_pos()
                                d_coor = (d_coor[0], d_coor[1])
                                intersects = utils.closed_segment_intersect(a_coor, b_coor, c_coor, d_coor)
                                if intersects and (a.i, b.i, c.i, d.i) not in intersecting_segments and \
                                        (c.i, d.i, a.i, b.i) not in intersecting_segments:
                                    intersecting_segments.append((a.i, b.i, c.i, d.i))
                    for c in [self.vertices[index] for index in a.out_neighbourhood]:
                        c_coor = c.im_pos()
                        c_coor = (c_coor[0], c_coor[1])
                        for d in [self.vertices[index] for index in c.out_neighbourhood]:
                            d_coor = d.im_pos()
                            d_coor = (d_coor[0], d_coor[1])
                            for e in [self.vertices[index] for index in d.out_neighbourhood]:
                                e_coor = e.im_pos()
                                e_coor = (e_coor[0], e_coor[1])
                                intersects = utils.closed_segment_intersect(a_coor, b_coor, d_coor, e_coor)
                                if intersects and (a.i, b.i, d.i, e.i) not in intersecting_segments and \
                                        (d.i, e.i, a.i, b.i) not in intersecting_segments:
                                    intersecting_segments.append((a.i, b.i, d.i, e.i))

        return intersecting_segments

    def find_zeta_intersections(self):

        intersecting_segments = []

        for a in self.vertices:
            a_coor = a.im_pos()
            a_coor = (a_coor[0], a_coor[1])
            for b in [self.vertices[index] for index in a.local_zeta_map['out_neighbourhood']]:
                if not a.is_edge_column and not b.is_edge_column:
                    b_coor = b.im_pos()
                    b_coor = (b_coor[0], b_coor[1])
                    for c in [self.vertices[index] for index in a.local_zeta_map['out_neighbourhood']]:
                        if not c.i == b.i:
                            c_coor = c.im_pos()
                            c_coor = (c_coor[0], c_coor[1])
                            for d in [self.vertices[index] for index in c.local_zeta_map['out_neighbourhood']]:
                                d_coor = d.im_pos()
                                d_coor = (d_coor[0], d_coor[1])
                                intersects = utils.closed_segment_intersect(a_coor, b_coor, c_coor, d_coor)
                                if intersects:
                                    if (a.i, b.i, c.i, d.i) not in intersecting_segments and \
                                            (a.i, b.i, d.i, c.i) not in intersecting_segments and \
                                            (b.i, a.i, c.i, d.i) not in intersecting_segments and \
                                            (b.i, a.i, d.i, c.i) not in intersecting_segments and \
                                            (c.i, d.i, a.i, b.i) not in intersecting_segments and \
                                            (c.i, d.i, b.i, a.i) not in intersecting_segments and \
                                            (d.i, c.i, a.i, b.i) not in intersecting_segments and \
                                            (d.i, c.i, b.i, a.i) not in intersecting_segments:
                                        intersecting_segments.append((a.i, b.i, c.i, d.i))
                    for c in [self.vertices[index] for index in a.local_zeta_map['out_neighbourhood']]:
                        for d in [self.vertices[index] for index in c.local_zeta_map['out_neighbourhood']]:
                            d_coor = d.im_pos()
                            d_coor = (d_coor[0], d_coor[1])
                            for e in [self.vertices[index] for index in d.local_zeta_map['out_neighbourhood']]:
                                e_coor = e.im_pos()
                                e_coor = (e_coor[0], e_coor[1])
                                intersects = utils.closed_segment_intersect(a_coor, b_coor, d_coor, e_coor)
                                if intersects:
                                    if (a.i, b.i, d.i, e.i) not in intersecting_segments and \
                                            (a.i, b.i, e.i, d.i) not in intersecting_segments and \
                                            (b.i, a.i, d.i, e.i) not in intersecting_segments and \
                                            (b.i, a.i, e.i, d.i) not in intersecting_segments and \
                                            (d.i, e.i, a.i, b.i) not in intersecting_segments and \
                                            (d.i, e.i, b.i, a.i) not in intersecting_segments and \
                                            (e.i, d.i, a.i, b.i) not in intersecting_segments and \
                                            (e.i, d.i, b.i, a.i) not in intersecting_segments:
                                        intersecting_segments.append((a.i, b.i, d.i, e.i))

        return intersecting_segments

    def find_anti_zeta_intersections(self):
        intersecting_segments = []

        for a in self.vertices:
            a_coor = a.im_pos()
            a_coor = (a_coor[0], a_coor[1])
            for b in [self.vertices[index] for index in a.local_zeta_map['anti_out_neighbourhood']]:
                if not a.is_edge_column and not b.is_edge_column:
                    b_coor = b.im_pos()
                    b_coor = (b_coor[0], b_coor[1])
                    for c in [self.vertices[index] for index in a.local_zeta_map['out_neighbourhood']]:
                        if not c.i == b.i:
                            c_coor = c.im_pos()
                            c_coor = (c_coor[0], c_coor[1])
                            for d in [self.vertices[index] for index in c.local_zeta_map['out_neighbourhood']]:
                                d_coor = d.im_pos()
                                d_coor = (d_coor[0], d_coor[1])
                                intersects = utils.closed_segment_intersect(a_coor, b_coor, c_coor, d_coor)
                                if intersects:
                                    if (a.i, b.i, c.i, d.i) not in intersecting_segments and \
                                            (a.i, b.i, d.i, c.i) not in intersecting_segments and \
                                            (b.i, a.i, c.i, d.i) not in intersecting_segments and \
                                            (b.i, a.i, d.i, c.i) not in intersecting_segments and \
                                            (c.i, d.i, a.i, b.i) not in intersecting_segments and \
                                            (c.i, d.i, b.i, a.i) not in intersecting_segments and \
                                            (d.i, c.i, a.i, b.i) not in intersecting_segments and \
                                            (d.i, c.i, b.i, a.i) not in intersecting_segments:
                                        intersecting_segments.append((a.i, b.i, c.i, d.i))
                    for c in [self.vertices[index] for index in a.local_zeta_map['anti_out_neighbourhood']]:
                        for d in [self.vertices[index] for index in c.local_zeta_map['out_neighbourhood']]:
                            d_coor = d.im_pos()
                            d_coor = (d_coor[0], d_coor[1])
                            for e in [self.vertices[index] for index in d.local_zeta_map['out_neighbourhood']]:
                                e_coor = e.im_pos()
                                e_coor = (e_coor[0], e_coor[1])
                                intersects = utils.closed_segment_intersect(a_coor, b_coor, d_coor, e_coor)
                                if intersects:
                                    if (a.i, b.i, d.i, e.i) not in intersecting_segments and \
                                            (a.i, b.i, e.i, d.i) not in intersecting_segments and \
                                            (b.i, a.i, d.i, e.i) not in intersecting_segments and \
                                            (b.i, a.i, e.i, d.i) not in intersecting_segments and \
                                            (d.i, e.i, a.i, b.i) not in intersecting_segments and \
                                            (d.i, e.i, b.i, a.i) not in intersecting_segments and \
                                            (e.i, d.i, a.i, b.i) not in intersecting_segments and \
                                            (e.i, d.i, b.i, a.i) not in intersecting_segments:
                                        intersecting_segments.append((a.i, b.i, d.i, e.i))

        return intersecting_segments

    def terminate_arc(self, i, j):
        if self.vertices[i].op_arc_pivot(j, self.vertices[i].district[self.vertices[i].n - 1]):
            if self.vertices[i].decrement_n():
                self.build_local_map([i] + self.vertices[i].district)
                return True
            else:
                return False
        else:
            return True

    def calc_vertex_parameters(self, i):
        vertex = self.vertices[i]
        vertex.alpha_angles = self.get_alpha_angles(i, selection_type='zeta')
        if vertex.alpha_angles is not None and not len(vertex.alpha_angles) == 0:
            vertex.alpha_max = max(vertex.alpha_angles)
            vertex.alpha_min = min(vertex.alpha_angles)
        else:
            vertex.alpha_max = 0
            vertex.alpha_min = 0
        vertex.theta_angles = self.get_theta_angles(i, selection_type='selective')
        if vertex.theta_angles is not None and not len(vertex.theta_angles) == 0:
            vertex.theta_max = max(vertex.theta_angles)
            vertex.theta_min = min(vertex.theta_angles)
            vertex.theta_angle_variance = utils.variance(vertex.theta_angles)
            vertex.theta_angle_mean = utils.mean_val(vertex.theta_angles)
        else:
            vertex.theta_max = 0
            vertex.theta_min = 0
            vertex.theta_angle_variance = 0
            vertex.theta_angle_mean = 0

    def calc_normalized_gamma(self):
        peak_gammas = []
        avg_gammas = []
        for vertex in self.vertices:
            if not vertex.is_in_precipitate and vertex.atomic_species == 'Al':
                peak_gammas.append(vertex.peak_gamma)
                avg_gammas.append(vertex.avg_gamma)
        peak_mean = utils.mean_val(peak_gammas)
        avg_mean = utils.mean_val(avg_gammas)
        peak_mean_diff = peak_mean - 0.3
        avg_mean_diff = avg_mean - 0.3
        for vertex in self.vertices:
            vertex.normalized_peak_gamma = vertex.peak_gamma - peak_mean_diff
            vertex.normalized_avg_gamma = vertex.avg_gamma - avg_mean_diff

    def calc_redshifts(self):
        self.total_redshift = 0
        self.matrix_redshift = 0
        self.particle_redshift = 0
        anti_graph = AntiGraph(self)
        for vertex in self.vertices:
            if not vertex.is_edge_column and not vertex.void:
                vertex.redshift = 0
                vertex.avg_central_separation = 0
                counter = 0
                for partner in vertex.partners:
                    vertex.redshift += self.get_redshift(vertex.i, partner)
                    vertex.avg_central_separation += self.get_separation(vertex.i, partner)
                    counter += 1
                for partner in anti_graph.vertices[vertex.i].partners:
                    vertex.redshift += self.get_redshift(vertex.i, partner)
                    vertex.avg_central_separation += self.get_separation(vertex.i, partner)
                    counter += 1
                if not counter == 0:
                    vertex.avg_central_separation /= counter
                    vertex.avg_redshift = vertex.redshift / counter
                self.total_redshift += vertex.redshift
                if vertex.is_in_precipitate:
                    self.particle_redshift += vertex.redshift
                else:
                    self.matrix_redshift += vertex.redshift

    def calc_all_parameters(self):
        for vertex in self.vertices:
            self.calc_vertex_parameters(vertex.i)
        self.calc_normalized_gamma()

    def check_species(self):
        for vertex in self.vertices:
            if vertex.advanced_species not in self.species_dict['advanced_species']:
                for key, item in self.species_dict['advanced_species'].items():
                    if item['atomic_species'] == vertex.atomic_species:
                        self.set_advanced_species(vertex.i, key)
                        break
                else:
                    self.set_advanced_species(vertex.i, 'Un_1')
            new_apv = {}
            for advanced_prob in self.species_dict['advanced_species']:
                if advanced_prob in vertex.advanced_probability_vector:
                    new_apv[advanced_prob] = vertex.advanced_probability_vector[advanced_prob]
                else:
                    new_apv[advanced_prob] = 0.0
            vertex.advanced_probability_vector = utils.normalize_dict(new_apv, 1)
            vertex.set_probability_from_advanced()
            vertex.determine_species_from_probability_vector()

    def check_integrity(self):
        for vertex in self.vertices:
            if vertex.i in vertex.district:
                logger.info('    {}'.format(vertex.i))

    def refresh_graph(self):
        self.build_local_maps()
        self.calc_all_parameters()
        self.evaluate_sub_categories()
        self.map_arcs()
        self.summarize_stats()

    def evaluate_sub_categories(self):
        for vertex in self.vertices:
            if vertex.atomic_species == 'Si':
                sub_graph = self.get_column_centered_subgraph(vertex.i)
                for mesh in sub_graph.meshes:
                    if len(mesh.vertices) == 4:
                        if mesh.vertices[1].atomic_species == 'Mg' and mesh.vertices[2].atomic_species == 'Cu' and mesh.vertices[3].atomic_species == 'Mg':
                            self.set_advanced_species(vertex.i, 'Si_2')
                            break
                else:
                    self.set_advanced_species(vertex.i, 'Si_1')
            elif vertex.atomic_species == 'Al':
                if vertex.is_in_precipitate:
                    self.set_advanced_species(vertex.i, 'Al_2')
                else:
                    self.set_advanced_species(vertex.i, 'Al_1')
            elif vertex.atomic_species == 'Mg':
                if vertex.alpha_max < 3.175:
                    self.set_advanced_species(vertex.i, 'Mg_1')
                else:
                    self. set_advanced_species(vertex.i, 'Mg_2')

    def calc_condensed_property_data(self, keys, filter_=None, recalc=False, include_un=False):

        if filter_ is None:
            filter_ = {
                'exclude_edge_columns': True,
                'exclude_matrix_columns': False,
                'exclude_particle_columns': False,
                'exclude_hidden_columns': False,
                'exclude_flag_1_columns': False,
                'exclude_flag_2_columns': False,
                'exclude_flag_3_columns': False,
                'exclude_flag_4_columns': False
            }

        if recalc:
            self.refresh_graph()

        data = []
        for vertex in self.vertices:
            if not vertex.void:
                edge = False
                if vertex.is_edge_column:
                    edge = True
                for neighbour in self.get_vertex_objects_from_indices(vertex.neighbourhood):
                    if neighbour.is_edge_column:
                        edge = True
                        break
                if not (filter_['exclude_edge_columns'] and edge):
                    if not (filter_['exclude_matrix_columns'] and not vertex.is_in_precipitate):
                        if not (filter_['exclude_particle_columns'] and vertex.is_in_precipitate):
                            if not (filter_['exclude_hidden_columns'] and not vertex.show_in_overlay):
                                if not (filter_['exclude_flag_1_columns'] and vertex.flag_1):
                                    if not (filter_['exclude_flag_2_columns'] and vertex.flag_2):
                                        if not (filter_['exclude_flag_3_columns'] and vertex.flag_3):
                                            if not (filter_['exclude_flag_4_columns'] and vertex.flag_4):

                                                data_item = {}
                                                for key in keys:
                                                    data_item[key] = getattr(vertex, key)
                                                    if key == 'theta_angle_mean':
                                                        data_item[key] += random.gauss(0, 0.1)
                                                if not (not include_un and vertex.atomic_species == 'Un'):
                                                    if not vertex.advanced_species == 'Cu_2':
                                                        if not vertex.theta_angle_mean < 0.7:
                                                            data.append(data_item)

        return data

    def map_arcs(self):
        self.arcs = []
        self.size = 0
        for vertex in self.vertices:
            for out_neighbour in self.get_vertex_objects_from_indices(vertex.out_neighbourhood):
                arc = Arc(len(self.arcs), vertex, out_neighbour)
                self.arcs.append(arc)
                self.size += 1

    def map_meshes(self):
        """Automatically generate a connected relational map of all meshes in graph.

        The index of a mesh is temporarily indexed during the mapping by the following algorithm: Take the indices of
        its corners, and circularly permute them such that the lowest index comes first. After the mapping is complete,
        these indices are replaced by the integers 0 to the number of meshes.

        """
        logger.info('mapping meshes')
        self.map_arcs()
        self.meshes = []
        for arc in self.arcs:
            mesh_1 = self.get_mesh(arc.vertex_a.i, arc.vertex_b.i, 0)
            mesh_2 = self.get_mesh(arc.vertex_b.i, arc.vertex_a.i, 0)
            if mesh_1 not in self.meshes:
                self.meshes.append(mesh_1)
            if mesh_2 not in self.meshes:
                self.meshes.append(mesh_2)
        self.volume = len(self.meshes)
        logger.info('Meshes mapped.')

    @staticmethod
    def determine_temp_index(mesh):
        return utils.make_int_from_list(utils.cyclic_sort(mesh.vertex_indices))

    def calc_chi(self):
        # Calc size
        self.size = 0
        for vertex in self.vertices:
            self.size += len(vertex.out_neighbourhood)

        # Calc chi (# weak arcs / # num strong arcs)
        num_weak_arcs = 0
        for vertex in self.vertices:
            if not vertex.is_edge_column:
                num_weak_arcs += len(vertex.out_semi_partners)
        if not self.size == 0:
            self.chi = num_weak_arcs / self.size
        else:
            self.chi = 1

    def summarize_stats(self):
        # Calc order
        self.order = len(self.vertices)

        # Calc size and chi
        self.calc_chi()

        # Calc average degree
        counted_columns = 0
        degrees = 0
        for vertex in self.vertices:
            if not vertex.is_edge_column:
                degrees += vertex.degree
                counted_columns += 1
        if not counted_columns == 0:
            self.avg_degree = degrees / counted_columns
        else:
            self.avg_degree = 0


class Mesh:

    def __init__(self, mesh_index, vertices):

        self.mesh_index = mesh_index
        self.vertices = vertices
        self.vertex_indices = []
        self.arcs = []
        self.angles = []
        self.angle_vectors = []
        self.surrounding_meshes = []

        self.is_enclosed = True
        self.is_consistent = True
        self.order = 0
        self.size = 0
        self.cm = (0, 0)

        for vertex in self.vertices:
            self.vertex_indices.append(vertex.i)
            self.order += 1
        for vertex in self.vertices:
            for out_neighbour in vertex.out_neighbourhood:
                if out_neighbour in self.vertex_indices:
                    index = self.vertex_indices.index(out_neighbour)
                    self.arcs.append(Arc(self.size, vertex, self.vertices[index]))
                    self.size += 1
        self.calc_cm()

    def __str__(self):

        string = ''

        for k, index in enumerate(self.vertex_indices):

            if self.vertices[utils.circularize_next_index(k + 1, len(self.vertices) - 1)].partner_query(k):
                end_left = '<'
            else:
                end_left = ''

            if self.vertices[k].partner_query(utils.circularize_next_index(k + 1, len(self.vertices) - 1)):
                end_right = '>'
            else:
                end_right = ''

            string += '{} {}-{} '.format(index, end_left, end_right)

        return string

    def __eq__(self, other):
        return utils.is_circularly_identical(self.vertex_indices, other.vertex_indices)

    def calc_cm(self):
        x_fit = 0
        y_fit = 0
        mass = len(self.vertices)
        for corner in self.vertices:
            x_fit += corner.im_coor_x
            y_fit += corner.im_coor_y
        x_fit = x_fit / mass
        y_fit = y_fit / mass
        self.cm = (x_fit, y_fit)

    def test_sidedness(self):
        if not self.order == 4:
            return False
        else:
            return True

    def get_ind_from_mother(self, i):

        for index, mother_index in enumerate(self.vertex_indices):
            if mother_index == i:
                sub_index = index
                break
        else:
            sub_index = -1
        return sub_index

    @staticmethod
    def neighbour_test(mesh_1, mesh_2):

        found = 0
        edge = []

        for corner in mesh_1.vertex_indices:

            if corner in mesh_2.vertex_indices:
                found += 1
                edge.append(corner)

            if found == 2:
                break

        else:

            return False

        assure_edge = True

        for k, corner in enumerate(mesh_1.vertex_indices):
            k_max = len(mesh_1.vertex_indices) - 1

            if corner == edge[0]:

                if mesh_1.vertex_indices[utils.circularize_next_index(k - 1, k_max)] == edge[1] or \
                        mesh_1.vertex_indices[utils.circularize_next_index(k + 1, k_max) == edge[1]]:
                    pass
                else:
                    assure_edge = False

            elif corner == edge[1]:

                if mesh_1.vertex_indices[utils.circularize_next_index(k - 1, k_max)] == edge[0] or \
                        mesh_1.vertex_indices[utils.circularize_next_index(k + 1, k_max) == edge[0]]:
                    pass
                else:
                    assure_edge = False

        for k, corner in enumerate(mesh_2.vertex_indices):
            k_max = len(mesh_2.vertex_indices) - 1

            if corner == edge[0]:

                if mesh_2.vertex_indices[utils.circularize_next_index(k - 1, k_max)] == edge[1] or \
                        mesh_2.vertex_indices[utils.circularize_next_index(k + 1, k_max) == edge[1]]:
                    pass
                else:
                    assure_edge = False

            elif corner == edge[1]:

                if mesh_2.vertex_indices[utils.circularize_next_index(k - 1, k_max)] == edge[0] or \
                        mesh_2.vertex_indices[utils.circularize_next_index(k + 1, k_max) == edge[0]]:
                    pass
                else:
                    assure_edge = False

        if assure_edge:
            return True, edge
        else:
            return False, edge


class MeshCenteredSubGraph:

    def __init__(self, super_graph, order, i, j):

        self.super_graph = super_graph
        self.order = order
        self.vertices = []
        self.arcs = []
        self.i = i
        self.j = j
        self.exists = True

        self.build()

    def build(self):
        logger.info('Building subgraph..')
        self.super_graph.map_arcs()
        if self.i not in self.super_graph.vertices[self.j].out_neighbourhood:
            if self.j not in self.super_graph.vertices[self.i].out_neighbourhood:
                self.exists = False
                logger.warning('Trying to build an undefined subgraph!')


class SubGraph:

    def __init__(self):
        self.vertices = []
        self.vertex_indices = []
        self.arcs = []
        self.meshes = []

        self.num_vertices = 0
        self.num_arcs = 0
        self.num_meshes = 0

        self.class_ = None
        self.configuration = None

    def finalize_init(self):
        self.redraw_edges()
        self.sort_meshes()
        self.summarize_stats()

    def summarize_stats(self):
        self.num_vertices = len(self.vertices)
        self.num_arcs = len(self.arcs)
        self.num_meshes = len(self.meshes)

    def sort_meshes(self):
        new_list = []
        for mesh in self.meshes:
            if mesh.vertex_indices[1] == self.vertices[0].district[0]:
                new_list.append(mesh)
                break
        else:
            new_list.append(self.meshes[0])
        closed = False
        if len(self.meshes) == 0:
            closed = True
        backup_counter = 0
        while not closed:
            for mesh in self.meshes:
                if mesh.vertex_indices[1] ==\
                        new_list[-1].vertex_indices[-1]:
                    new_list.append(mesh)

                    if new_list[-1].vertex_indices[-1] ==\
                            new_list[0].vertex_indices[1]:
                        closed = True
                    break
            backup_counter += 1
            if backup_counter > 26:
                break
        if len(new_list) == len(self.meshes) and closed:
            self.meshes = new_list

    def add_vertex(self, vertex):
        self.vertices.append(vertex)
        self.vertex_indices.append(vertex.i)
        self.num_vertices += 1

    def add_mesh(self, mesh):
        self.meshes.append(mesh)
        self.num_meshes += 1
        for vertex in self.meshes[-1].vertices:
            if vertex.i not in self.vertex_indices:
                self.add_vertex(vertex)

    def get_ind_from_mother(self, i):
        for index, mother_index in enumerate(self.vertex_indices):
            if mother_index == i:
                sub_index = index
                break
        else:
            sub_index = -1
        return sub_index

    def remove_vertex(self, vertex_index):
        raise NotImplemented

    def increase_h(self, i):
        i = self.get_ind_from_mother(i)
        if not i == -1:
            changed = self.vertices[i].increment_species_index()
        else:
            changed = False
        return changed

    def decrease_h(self, i):
        i = self.get_ind_from_mother(i)
        if not i == -1:
            changed = self.vertices[i].decrement_species_index()
        else:
            changed = False
        return changed

    def add_arc(self, j, vertex_a, vertex_b):
        arc = Arc(j, vertex_a, vertex_b)
        self.arcs.append(arc)
        self.num_arcs += 1

    def remove_arcs(self, arc_index):
        raise NotImplemented

    def redraw_edges(self):
        self.arcs = []
        self.num_arcs = 0
        for vertex in self.vertices:
            for out_neighbour in vertex.out_neighbourhood:
                if out_neighbour in self.vertex_indices:
                    self.add_arc(self.num_arcs, vertex, self.vertices[self.get_ind_from_mother(out_neighbour)])


class AntiGraph:

    def __init__(self, graph):

        self.graph = graph
        self.vertices = copy.deepcopy(graph.vertices)
        self.arcs = []
        self.size = 0

        self.build()

    def build(self):
        for i, vertex in enumerate(self.graph.vertices):
            if not vertex.is_edge_column and not vertex.void:
                sub_graph = self.graph.get_column_centered_subgraph(vertex.i)
                for mesh in sub_graph.meshes:
                    self.vertices[i].op_arc_pivot(mesh.vertex_indices[1], mesh.vertex_indices[2])
        self.graph = AtomicGraph(self.graph.scale, self.graph.species_dict)
        for vertex in self.vertices:
            self.graph.add_vertex(vertex)
        self.graph.build_local_maps()
        self.graph.summarize_stats()

    def map_arcs(self):
        self.arcs = []
        self.size = 0
        for vertex in self.vertices:
            if not vertex.void:
                for out_neighbour in self.graph.get_vertex_objects_from_indices(vertex.out_neighbourhood):
                    if not out_neighbour.void:
                        arc = Arc(len(self.arcs), vertex, out_neighbour)
                        self.arcs.append(arc)
                        self.size += 1



