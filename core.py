# By Haakon Tvedt @ NTNU
# Contributors:
"""Module for the 'Project' class that handles a *project instance*."""

# Program imports:
import utils
import graph_2
import compatibility
import column_characterization
import conversion_tools
# External imports:
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import scipy.ndimage
import dm3_lib as dm3
import sys
import copy
import time
import pickle
import configparser
import logging
# Instantiate logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Project:
    """The main API through which to build and access the data extracted from HAADF-STEM images.

    :param filename_full: The full path and/or relative path and filename of the .dm3 image to import. A project can
        be instantiated with filename_full='empty', but this is only meant to be used as a placeholder.
        Use Project.import(filename, file_type) to instantiate projects from other files than dm3.
    :param debug_obj: (Optional, default=None) An instance of an AutomAl 6000 GUI.MainUI().
    :type filename_full: string
    :type debug_obj: GUI.MainUI()

    """

    # Version
    version = [0, 2, 1]

    default_species_dict = {
        'advanced_species': {
            'Si_1': {'n': 3, 'atomic_species': 'Si', 'color': (255, 20, 20), 'description': 'Beta-pprime Si'},
            'Si_2': {'n': 3, 'atomic_species': 'Si', 'color': (160, 40, 40), 'description': 'Q-prime Si'},
            'Si_3': {'n': 3, 'atomic_species': 'Si', 'color': (120, 0, 0), 'description': 'interstitial beta-prime Si'},
            'Cu_1': {'n': 3, 'atomic_species': 'Cu', 'color': (180, 180, 0), 'description': '3-fold symmetry Cu'},
            'Cu_2': {'n': 4, 'atomic_species': 'Cu', 'color': (220, 220, 20), 'description': '4-fold symmetry Cu'},
            'Al_1': {'n': 4, 'atomic_species': 'Al', 'color': (20, 255, 20), 'description': 'FCC matrix Al'},
            'Al_2': {'n': 4, 'atomic_species': 'Al', 'color': (40, 140, 40), 'description': 'Precipitate Al'},
            'Mg_1': {'n': 5, 'atomic_species': 'Mg', 'color': (138, 43, 226), 'description': 'Mg with alpha_max < 3.12'},
            'Mg_2': {'n': 5, 'atomic_species': 'Mg', 'color': (76, 24, 151), 'description': 'Mg with alpha_max > 3.12'},
            'Un_1': {'n': 3, 'atomic_species': 'Un', 'color': (20, 20, 255), 'description': 'Unknown species'},
            'Un_2': {'n': 4, 'atomic_species': 'Un', 'color': (0, 0, 0), 'description': 'Vacant/void column'}
        },
        'atomic_species': {
            'Si': {'atomic_radii': 117.50, 'color': (255, 20, 20)},
            'Cu': {'atomic_radii': 127.81, 'color': (180, 180, 0)},
            'Al': {'atomic_radii': 143.00, 'color': (20, 255, 20)},
            'Mg': {'atomic_radii': 160.00, 'color': (138, 43, 226)},
            'Un': {'atomic_radii': 100.00, 'color': (20, 20, 255)}
        }
    }

    # District size
    district_size = 10
    al_lattice_const = 404.95

    def __init__(self, filename_full, debug_obj=None, species_dict=None):

        self.filename_full = filename_full
        self.im_mat = None
        self.fft_im_mat = None
        self.scale = 1
        self.im_height = 0
        self.im_width = 0
        self.version_saved = None
        self.starting_index = None

        # In AutomAl 6000, each column is modelled by a single atomic species, more advanced categorization can be
        # applied however. Each advanced category must map to one of the simple categories and also to a symmetry.
        parser = configparser.ConfigParser()
        parser.read('config.ini')
        if species_dict is None:
            self.species_dict = Project.default_species_dict
        else:
            self.species_dict = species_dict

        # For communicating with the interface, if any:
        self.debug_obj = debug_obj
        self.debug_mode = False

        self.im_meta_data = {}
        self.has_been_up_sampled = False
        self.maps = []
        if not (filename_full == 'Empty' or filename_full == 'empty'):
            dm3f = dm3.DM3(self.filename_full)
            self.im_mat = dm3f.imagedata
            (self.scale, _) = dm3f.pxsize
            self.scale = 1000 * self.scale  # Scale is now in nm/pixel
            if self.scale > 7.0:
                logger.info('Up-sampling image to prevent over-granulation in column detection!')
                self.has_been_up_sampled = True
                self.im_mat = scipy.ndimage.zoom(self.im_mat, 2, order=1)
                self.scale = self.scale / 2
            (self.im_height, self.im_width) = self.im_mat.shape
            self.im_mat = utils.normalize_static(self.im_mat)
            self.fft_im_mat = utils.gen_fft(self.im_mat)

        # Data matrices:
        self.search_mat = copy.deepcopy(self.im_mat)

        # Counting and statistical variables.
        self.num_columns = 0
        self.pixel_average = 0
        if not (filename_full == 'Empty' or filename_full == 'empty'):
            self.calc_avg_pixel_value()

        # These are hyper-parameters of the algorithms. See the documentation.
        self.threshold = 0.2586
        self.search_size = 10
        self.r = int(100 / self.scale)
        self.overhead = int(7 * (self.r / 10))

        # Initialize an empty graph
        self.graph = graph_2.AtomicGraph(self.scale, active_model=None, species_dict=self.species_dict, district_size=self.district_size)

        logger.info('Generated instance from {}'.format(filename_full))

    def __str__(self):
        return self.report()

    def report(self, supress_log=True):
        """Build a string representation of the current instance.

        """

        string = 'Project summary:\n'
        string += '    General:\n'
        string += '        Number of columns: {}\n'.format(self.num_columns)
        string += '        Alloy: {}\n'.format(self.get_alloy_string())
        string += '        Precipitate column composition:\n'
        for key, value in self.get_precipitate_composition().items():
            string += '            {}:\t{:.4f}\n'.format(key, value)
        string += '        Image column composition:\n'
        for key, value in self.get_image_composition().items():
            string += '            {}:\t{:.4f}\n'.format(key, value)
        string += '    Image:\n'
        string += '        Dimension (height, width): ({}, {})\n'.format(self.im_height, self.im_width)
        string += '        Average pixel intensity: {:.4f}\n'.format(self.pixel_average)
        string += '        Scale (pm/pixel): {:.7f}\n'.format(self.scale)
        string += '    Atomic graph:\n'
        for line in self.graph.report().splitlines(keepends=True):
            string += '    ' + line

        if supress_log:
            return string
        else:
            logger.info(string)
            return None

    def vertex_report(self, i, supress_log=False):
        string = self.graph.vertices[i].report()
        if supress_log:
            return string
        else:
            logger.info(string)
            return None

    def get_alloy_string(self):
        """Generate a string representation of the alloy based on the alloys present in the species dictionary.

        :returns A species string on the form Al-Mg-Si.
        :rtype string:

        """
        species_string = ''
        species = set()
        for sp in self.species_dict['atomic_species']:
            species.add(sp)
        species.discard('Un')
        if 'Al' in species:
            species_string += 'Al-'
            species.discard('Al')
        if 'Mg' in species:
            species_string += 'Mg-'
            species.discard('Mg')
        if 'Si' in species:
            species_string += 'Si-'
            species.discard('Si')
        if 'Cu' in species:
            species_string += 'Cu-'
            species.discard('Cu')
        for remaining in species:
            species_string += '{}-'.format(remaining)
        return species_string[:-1]

    def get_image_composition(self):

        image_composition = {}
        columns = 0

        for vertex in self.graph.vertices:
            if vertex.atomic_species in image_composition:
                image_composition[vertex.atomic_species] += 1
            else:
                image_composition[vertex.atomic_species] = 1
            columns += 1

        for key, value in image_composition.items():
            image_composition[key] = image_composition[key] / columns

        return image_composition

    def get_precipitate_composition(self):

        precipitate_composition = {}
        precipitate_columns = 0

        for vertex in self.graph.vertices:
            if vertex.is_in_precipitate:
                if vertex.atomic_species in precipitate_composition:
                    precipitate_composition[vertex.atomic_species] += 1
                else:
                    precipitate_composition[vertex.atomic_species] = 1
                precipitate_columns += 1

        for key, value in precipitate_composition.items():
            precipitate_composition[key] = precipitate_composition[key] / precipitate_columns

        return precipitate_composition

    def get_precipitate_packing_fraction(self):

        pass

    def get_precipitate_displacement(self):

        pass

    def save(self, filename_full):
        """Save the current project as a pickle file.

        :param filename_full: Path and name of save-file. The project will be pickled as *filename_full* without any
            file-extension identifier.
        :type filename_full: string

        """

        with open(filename_full, 'wb') as f:
            self.debug_obj = None
            self.version_saved = self.version
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
        logger.info('Saved {}'.format(filename_full))

    def experimental_save(self, filename_full):
        pass

    @staticmethod
    def load(filename_full):
        """Load an instance from a pickle-file.

        :param filename_full: Path-name of the file to be loaded.
        :type filename_full: string

        :return: project instance.
        :rtype: core.Project

        """
        with open(filename_full, 'rb') as f:
            try:
                obj = pickle.load(f)
            except:
                obj = None
                logger.error('Failed to load save-file!')
            else:
                if not obj.version_saved == Project.version:
                    logger.info('Attempted to load un-compatible save-file. Running conversion script...')
                    obj = compatibility.convert(obj, obj.version_saved)
                    if obj is None:
                        logger.error('Conversion unsuccessful, failed to load save-file!')
                    else:
                        logger.info('Conversion successful, loaded {}'.format(filename_full))
                else:
                    logger.info('Loaded {}'.format(filename_full))
        return obj

    @staticmethod
    def experimental_load(filename_full):
        pass

    def column_detection(self, search_type='s', plot=False):
        """Column detection algorithm.

        The search type indicates the termination condition of the column detection (see table).

        ----------------------------------- ------------------------------------------------------
        search_type                         explanation
        ----------------------------------- ------------------------------------------------------
        s                                   Search until self.search_size columns are found.
        t                                   Search until self.search_mat.max() < self.threshold.
        o                                   Undefined.
        ----------------------------------- ------------------------------------------------------

        This algorithm will check the column state and continue appropriately if there are already columns
        present. It will also "roll back" columns if the threshold value is higher than the peak gamma of the
        already existing columns. This means that if the initial threshold is guessed too low, then it will
        remove columns.

        :param search_type: (Optional, default='s') The algorithm termination condition.
        :type search_type: str
        :param plot: (Optional, default=False) Plotting on/off.
        :type plot: bool

        """
        time_1 = time.time()
        search_peak_values = []
        actual_peak_values = []
        avg_values = []
        if len(self.graph.vertices) == 0:
            logger.info('Starting column detection. Search mode is \'{}\''.format(search_type))
            self.search_mat = copy.deepcopy(self.im_mat)
            cont = True
        else:
            logger.info('Continuing column detection. Search mode is \'{}\''.format(search_type))
            min_val = 2
            for vertex in self.graph.vertices:
                if vertex.peak_gamma < min_val:
                    min_val = vertex.peak_gamma
            if min_val < self.threshold:
                logger.info('Columns overdetected. Rolling back..')
                new_graph = graph_2.AtomicGraph(
                    self.scale,
                    active_model=self.graph.active_model,
                    species_dict=self.graph.species_dict
                )
                self.num_columns = 0
                self.search_mat = copy.deepcopy(self.im_mat)
                for vertex in self.graph.vertices:
                    if vertex.peak_gamma > self.threshold:
                        self.num_columns += 1
                        search_peak_values.append(self.search_mat.max())
                        actual_peak_values.append(vertex.peak_gamma)
                        avg_values.append(vertex.avg_gamma)
                        new_graph.add_vertex(vertex)
                        self.search_mat = utils.delete_pixels(
                            self.search_mat,
                            int(vertex.im_coor_x),
                            int(vertex.im_coor_y),
                            self.r + self.overhead
                        )
                self.graph = new_graph
                cont = False
            else:
                self.redraw_search_mat()
                for vertex in self.graph.vertices:
                    actual_peak_values.append(vertex.peak_gamma)
                    avg_values.append(vertex.avg_gamma)
                cont = True

        counter = self.num_columns
        original_counter = copy.deepcopy(counter)

        while cont:

            pos = np.unravel_index(self.search_mat.argmax(), (self.im_height, self.im_width))
            max_val = self.search_mat[pos]
            search_peak_values.append(max_val)

            x_fit, y_fit = utils.cm_fit(self.im_mat, pos[1], pos[0], self.r)
            x_fit_int = int(x_fit)
            y_fit_int = int(y_fit)

            self.search_mat = utils.delete_pixels(self.search_mat, x_fit_int, y_fit_int, self.r + self.overhead)

            vertex = graph_2.Vertex(counter, x_fit, y_fit, self.r, self.scale, parent_graph=self.graph)
            vertex.avg_gamma, vertex.peak_gamma = utils.circular_average(self.im_mat, x_fit_int, y_fit_int, self.r)
            actual_peak_values.append(vertex.peak_gamma)
            avg_values.append(vertex.avg_gamma)

            self.graph.add_vertex(vertex)

            self.num_columns += 1
            counter += 1

            if search_type == 's':
                if counter >= self.search_size:
                    cont = False
            elif search_type == 't':
                if np.max(self.search_mat) < self.threshold:
                    cont = False
            else:
                logger.error('Invalid search_type')

        self.find_edge_columns()

        time_2 = time.time()

        logger.info('Column detection complete! Found {} columns in {} seconds.'.format(
            self.num_columns - original_counter,
            time_2 - time_1
        ))

        if plot:
            fig = plt.figure(constrained_layout=True)
            gs = GridSpec(1, 1, figure=fig)
            ax_values = fig.add_subplot(gs[0, 0])

            ax_values.plot(
                range(0, len(search_peak_values)),
                search_peak_values,
                c='k',
                label='Search mat peak intensity'
            )
            ax_values.plot(
                range(0, len(search_peak_values)),
                actual_peak_values,
                c='r',
                label='Vertex peak intensity'
            )
            ax_values.plot(
                range(0, len(search_peak_values)),
                avg_values,
                c='b',
                label='Vertex average intensity'
            )
            ax_values.plot(
                range(0, len(search_peak_values)),
                [self.pixel_average] * len(search_peak_values),
                c='g',
                label='Average image pixel intensitiy'
            )
            ax_values.set_title('Column detection summary')
            ax_values.set_xlabel('# Column')
            ax_values.set_ylabel('Pixel intensity')
            ax_values.legend()

            plt.show()

    def column_characterization(self, starting_index, search_type=0, ui_obj=None):
        """Column characterization algorithm.

        Assumes a *starting_index*, that is taken to be an Al column in the matrix. The *search_type* enables access to
        the different sub-proccesses of the algorithm. These are:

        ===================     ========================================================================
        :code:`search_type`     Process
        ===================     ========================================================================
        0                       The full column characterization
        1                       The most important parts of the algorithm
        2                       Finalize the parameter state
        3                       Spatial mapping
        4                       Identify edge columns
        5                       Not in use
        6                       Zeta analysis
        7                       Active model predictions from alpha angles
        8                       Particle detection
        9                       Calculate normalized gamma attributes
        10                      Evaluate advanced species
        11                      Active model predictions from all attributes
        12                      Reset all probability vectors
        13                      Reset user-set columns
        14                      Look for intersections
        15                      Not in use
        16                      Build the vertex maps of the graph (including out-neighbourhoods)
        17                      Build the vertex maps of the graph (excluding out-neighbourhoods)
        18                      Run untangling (including strong)
        19                      Run untangling (excluding strong)
        20                      Experimental zeta analysis
        21                      Calculate area gamma
        22                      Not in use
        23                      Not in use
        ===================     ========================================================================

        Columns that have been set by the user will not be have their species changed by the algorithm. Apply
        :code:`search_type=13` to reset columns if needed.

        :param starting_index: Index of a vertex that is a part of the matrix.
        :param search_type: (optional, default=0) Which sub-process to access.
        :type starting_index: int
        :type search_type: int

        """

        sys.setrecursionlimit(10000)

        if search_type == 0:
            self.column_characterization(starting_index, search_type=1, ui_obj=ui_obj)
            self.column_characterization(starting_index, search_type=2, ui_obj=ui_obj)

        elif search_type == 1:
            logger.info('Doing the basics...')
            # Spatial map:
            self.column_characterization(starting_index, search_type=3, ui_obj=ui_obj)
            # Detect edges:
            self.column_characterization(starting_index, search_type=4, ui_obj=ui_obj)
            # Map connectivity:
            self.column_characterization(starting_index, search_type=16, ui_obj=ui_obj)
            # Alpha model:
            self.column_characterization(starting_index, search_type=7, ui_obj=ui_obj)
            # Map connectivity
            self.column_characterization(starting_index, search_type=16, ui_obj=ui_obj)
            # Find particle:
            self.column_characterization(starting_index, search_type=8, ui_obj=ui_obj)
            # Calc gamma:
            self.column_characterization(starting_index, search_type=9, ui_obj=ui_obj)
            # Advanced zeta analysis:
            self.column_characterization(starting_index, search_type=6, ui_obj=ui_obj)
            # Untangle
            self.column_characterization(starting_index, search_type=18, ui_obj=ui_obj)
            # Composite model:
            self.column_characterization(starting_index, search_type=11, ui_obj=ui_obj)
            # Find particle:
            self.column_characterization(starting_index, search_type=8, ui_obj=ui_obj)
            # Calc area gamma:
            self.column_characterization(starting_index, search_type=21, ui_obj=ui_obj)
            # Calc gamma:
            self.column_characterization(starting_index, search_type=9, ui_obj=ui_obj)
            # Untangle
            self.column_characterization(starting_index, search_type=18, ui_obj=ui_obj)
            # Advanced zeta analysis:
            self.column_characterization(starting_index, search_type=6, ui_obj=ui_obj)
            logger.info('Basics done')

        elif search_type == 2:
            logger.info('Analyzing results...')
            # Find particle:
            self.column_characterization(starting_index, search_type=8, ui_obj=ui_obj)
            # Calc gamma:
            self.column_characterization(starting_index, search_type=9, ui_obj=ui_obj)
            # Refresh graph parameters
            self.graph.refresh_graph()

        elif search_type == 3:
            # Run spatial mapping
            logger.info('Mapping spatial locality...')
            column_characterization.determine_districts(self.graph)
            logger.info('Spatial mapping complete.')

        elif search_type == 4:
            # Identify edge columns
            logger.info('Finding edge columns....')
            column_characterization.find_edge_columns(self.graph, self.im_width, self.im_height)
            for vertex in self.graph.vertices:
                if vertex.is_edge_column:
                    self.graph.set_advanced_species(vertex.i, 'Al_1')
            logger.info('Edge columns found.')

        elif search_type == 5:
            # not in use
            pass

        elif search_type == 6:
            # Advanced zetas
            logger.info('Running zeta analysis...')
            column_characterization.zeta_analysis(self.graph)
            logger.info('zeta\'s set.')

        elif search_type == 7:
            # Applying alpha model
            logger.info('Calculating probabilities from alpha attributes...')
            column_characterization.apply_alpha_model(self.graph)
            logger.info('Calculated probabilities from alpha attributes.')

        elif search_type == 8:
            # Particle detection
            logger.info('Finding particle...')
            column_characterization.particle_detection(self.graph)
            logger.info('Found particle.')

        elif search_type == 9:
            # Determine normalized intensities
            logger.info('Finding normalized intensities...')
            self.normalize_gamma()
            logger.info('Found intensities.')

        elif search_type == 10:
            # Evaluate sub-species:
            logger.info('Evaluating advanced species...')
            self.graph.evaluate_sub_categories()
            logger.info('Advanced species set.')

        elif search_type == 11:
            # Applying composite model
            logger.info('Calculating probabilities from all attributes...')
            column_characterization.apply_composite_model(self.graph)
            logger.info('Calculated probabilities from all attributes.')

        elif search_type == 12:
            # reset probs
            logger.info('Resetting probability vectors with zero bias...')
            for vertex in self.graph.vertices:
                if not vertex.void and not vertex.is_set_by_user and not vertex.is_edge_column:
                    vertex.reset_probability_vector()
            logger.info('Probability vectors reset.')

        elif search_type == 13:
            # Reset user-input
            logger.info('Resetting user-set columns...')
            for i in range(0, self.num_columns):
                if self.graph.vertices[i].is_set_by_user:
                    self.graph.vertices[i].is_set_by_user = False
            logger.info('User-set columns was re-set.')

        elif search_type == 14:
            # Locate and remove edge intersections
            logger.info('Looking for intersections')
            intersections = self.graph.find_intersections()
            num_intersections = len(intersections)
            column_characterization.arc_intersection_denial(self.graph)
            intersections = self.graph.find_intersections()
            if ui_obj is not None:
                ui_obj.update_overlay()
                ui_obj.update_graph()
            logger.info('Found {} intersections'.format(num_intersections))
            logger.info('{} literal intersections still remain'.format(len(intersections)))

        elif search_type == 15:
            pass

        elif search_type == 16:
            # Run local graph mapping
            logger.info('Mapping vertex connectivity...')
            self.graph.build_local_maps(build_out=True)
            self.graph.build_local_zeta_maps(build_out=True)
            logger.info('Vertices mapped.')

        elif search_type == 17:
            # Run local graph mapping
            logger.info('Mapping vertex connectivity...')
            self.graph.build_local_maps(build_out=False)
            self.graph.build_local_zeta_maps(build_out=False)
            logger.info('Vertices mapped.')

        elif search_type == 18:
            # Untangle (strong)
            logger.info('Running untangling algorithm')
            column_characterization.untangle(self.graph, ui_obj=ui_obj)
            logger.info('Untangling complete')

        elif search_type == 19:
            # Untangle (weak)
            logger.info('Running untangling algorithm')
            column_characterization.untangle(self.graph, ui_obj=ui_obj, strong=False)
            logger.info('Untangling complete')

        elif search_type == 20:
            # Binary zeta analysis
            logger.info('Running binary zeta analysis')
            column_characterization.zeta_analysis(
                self.graph,
                starting_index,
                self.graph.vertices[starting_index].zeta,
                method='binary'
            )
            logger.info('Zeta analysis complete.')

        elif search_type == 21:
            # Calc area gamma
            logger.info('Calculating area gamma')
            self.calc_area_gamma()
            logger.info('Area gamma calculated')

        elif search_type == 22:
            pass

        elif search_type == 23:
            pass

        else:
            logger.error('No such search type!')

    def find_edge_columns(self):
        """Locate vertices that are close to the edge of the image.

        These vertices get special treatment throughout the program because information about their surroundings will
        be incomplete. This method will find all vertices that are within a distance 6 * self.r from the edge, and set
        the field self.graph.vertices[i].is_edge_column = True.

        """

        for vertex in self.graph.vertices:

            x_coor = vertex.im_coor_x
            y_coor = vertex.im_coor_y
            margin = 4 * self.r

            if x_coor < margin or x_coor > self.im_width - margin - 1 or y_coor < margin or y_coor > self.im_height - margin - 1:
                self.graph.vertices[vertex.i].is_edge_column = True
            else:
                self.graph.vertices[vertex.i].is_edge_column = False

    def normalize_gamma(self):
        """Find the mean of the intensity of the Al-matrix, and scale all intensities.

        Scale all intensities such that the mean of the Al-matrix is a fixed point. Store the result in each vertex *i*
        :code:`self.graph.vertices[i].normalized_peak_gamma` and
        :code:`self.graph.vertices[i].normalized_avg_gamma` fields.

        """
        self.graph.calc_normalized_gamma()

    def calc_avg_gamma(self):
        """Calculate average intensity for every vertex based on image information.

        """
        if self.graph.order > 0:
            for vertex in self.graph.vertices:
                vertex.avg_gamma, vertex.peak_gamma = utils.circular_average(
                    self.im_mat,
                    int(vertex.im_coor_x),
                    int(vertex.im_coor_y),
                    int(self.r)
                )

    def calc_area_gamma(self):
        """Calculate the area averaged pixel intensity"""
        if self.graph.order > 0:
            for vertex in self.graph.vertices:
                vertex.area_gamma, _ = utils.circular_average(
                    self.im_mat,
                    int(vertex.im_coor_x),
                    int(vertex.im_coor_y),
                    int(0.5 * self.species_dict['atomic_species'][vertex.atomic_species]['atomic_radii'] / self.scale)
                )

    def calc_avg_pixel_value(self):
        """Calculate the average pixel intensity of the image."""
        pixel_avg = 0
        for i in range(0, self.im_width):
            for j in range(0, self.im_height):
                pixel_avg += self.im_mat[j, i]
        pixel_avg /= self.im_width * self.im_height
        self.pixel_average = pixel_avg

    def redraw_search_mat(self):
        """Redraw the search matrix.

        """

        self.search_mat = copy.deepcopy(self.im_mat)
        if self.num_columns > 0:
            for vertex in self.graph.vertices:
                self.search_mat = utils.delete_pixels(
                    self.search_mat,
                    int(vertex.im_coor_x),
                    int(vertex.im_coor_y),
                    self.r + self.overhead)

    def reset_graph(self):
        """Reset the graph (delete all vertices)"""
        self.num_columns = 0
        self.graph = graph_2.AtomicGraph(self.scale, active_model=None, species_dict=self.species_dict, district_size=self.district_size)
        self.search_mat = copy.deepcopy(self.im_mat)
        self.starting_index = None

    def get_im_length_from_spatial(self, spatial_length):
        """Returns the a spatial length in image pixel length.

        :param spatial_length: A spatial length in pm
        :type spatial_length: float

        :returns An image length in px:
        :rtype float:

        """
        return self.scale * spatial_length

    def get_spatial_length_from_im(self, im_length):
        """Returns the an image length in spatial length.

        :param spatial_length: An image length in px
        :type spatial_length: float

        :returns A spatial length in pm:
        :rtype float:

        """
        return im_length / self.scale

    def make_heat_map(self, attribute, kernel_size=1, kernel_step_size=1, measure_type='variance', kernel_type='square', title='heat map'):
        """

        :param attribute:
        :type attribute: str
        :param kernel_size:
        :param kernel_step_size:
        :param kernel_shape:

        :returns A heat map of the image:
        :rtype numpy.array:

        """
        logger.info('Generating heat map (this may take a long time).')
        time_1 = time.time()
        self.graph.build_cities()
        heat_map = np.zeros((self.im_height, self.im_width), dtype=np.float)
        kernel_pixel_size = int(kernel_size * self.al_lattice_const / self.scale)
        x_steps = int(self.im_width/kernel_step_size)
        y_steps = int(self.im_height/kernel_step_size)
        warning_triggered = False
        for x_step in range(0, x_steps):
            x = x_step * kernel_step_size
            print('{:.2f}%'.format(100 * x / (x_steps - 1)))
            for y_step in range(0, y_steps):
                y = y_step * kernel_step_size
                attribute_list = []
                for vertex in self.graph.vertices:
                    if x - kernel_pixel_size < vertex.im_coor_x < x + kernel_pixel_size and \
                            y - kernel_pixel_size < vertex.im_coor_y < y + kernel_pixel_size:
                        for citizen in self.graph.get_vertex_objects_from_indices(list(vertex.city)):
                            if x - kernel_pixel_size < citizen.im_coor_x < x + kernel_pixel_size and \
                                    y - kernel_pixel_size < citizen.im_coor_y < y + kernel_pixel_size:
                                attribute_list.append(getattr(citizen, attribute))
                        break
                if measure_type.lower() == 'variance':
                    measure = np.var(attribute_list)
                elif measure_type.lower() == 'mean':
                    measure = np.mean(attribute_list)
                elif measure_type.lower() == 'sum':
                    measure = np.sum(attribute_list)
                else:
                    if not warning_triggered:
                        logger.warning('Project.make_heat_map(): Unrecognized measure type!')
                        warning_triggered = True
                    measure = np.var(attribute_list)

                heat_map[y, x] = measure
        time_2 = time.time()
        total_time = time_2 - time_1
        logger.info(
            'Heat map summary:\n    Kernel size: {} * a (pm)\n    Step size: {} (px)\n'
            '    Attribute: {}\n    Measure type: {}\n    Total time: {:.4f} (s)\n'.format(
                kernel_size,
                kernel_step_size,
                attribute,
                measure_type,
                total_time,
            )
        )
        details = {
            'title': title,
            'heat_mat': heat_map,
            'kernel': {'size': kernel_size, 'step_size': kernel_step_size},
            'attribute': attribute,
            'measure': measure_type,
            'min': heat_map.min(),
            'max': heat_map.max()
        }
        for i, map_ in enumerate(self.maps):
            if map_['title'] == title:
                logger.warning('There is already a heat map with this title. Overwriting...')
                self.maps[i] = details
                break
        else:
            self.maps.append(details)

        return heat_map

    def make_heat_map_2(self, attribute, kernel_size=1, kernel_step_size=1, measure_type='variance', kernel_type='square', title='heat map'):
        logger.info('Generating heat map (this may take a long time).')
        time_1 = time.time()
        columns = [[[]] * self.im_height] * self.im_width
        heat_map = np.zeros((self.im_height, self.im_width), dtype=np.float)
        kernel_pixel_size = int(kernel_size * self.al_lattice_const / self.scale)
        for vertex in self.graph.vertices:
            value = getattr(self.graph.vertices[vertex.i], attribute)
            x_min = max(0, int(vertex.im_coor_x) - kernel_pixel_size)
            x_max = min(self.im_width, int(vertex.im_coor_x) + kernel_pixel_size)
            y_min = max(0, int(vertex.im_coor_y) - kernel_pixel_size)
            y_max = min(self.im_height, int(vertex.im_coor_y) + kernel_pixel_size)
            for x in range(x_min, x_max):
                for y in range(y_min, y_max):
                    columns[y][x].append(value)
        if measure_type.lower() == 'variance':
            for x in range(0, self.im_width):
                print('{:.2f}%'.format(100 * x / (self.im_width - 1)))
                for y in range(0, self.im_height):
                    heat_map[y, x] = np.var(columns[y][x])
        elif measure_type.lower() == 'mean':
            for x in range(0, self.im_width):
                for y in range(0, self.im_height):
                    heat_map[y, x] = np.mean(columns[y][x])
        else:
            for x in range(0, self.im_width):
                for y in range(0, self.im_height):
                    heat_map[y, x] = np.sum(columns[y][x])
        time_2 = time.time()
        total_time = time_2 - time_1
        logger.info(
            'Heat map summary:\n    Kernel size: {} * a (pm)\n    Step size: {} (px)\n'
            '    Attribute: {}\n    Measure type: {}\n    Total time: {:.4f} (s)\n'.format(
                kernel_size,
                kernel_step_size,
                attribute,
                measure_type,
                total_time,
            )
        )
        details = {
            'title': title,
            'heat_mat': heat_map,
            'kernel': {'size': kernel_size, 'step_size': kernel_step_size},
            'attribute': attribute,
            'measure': measure_type,
            'min': heat_map.min(),
            'max': heat_map.max()
        }
        for i, map_ in enumerate(self.maps):
            if map_['title'] == title:
                logger.warning('There is already a heat map with this title. Overwriting...')
                self.maps[i] = details
                break
        else:
            self.maps.append(details)

        return heat_map

    def export_svg(self, savefile, filter_, species_type='atomic_species', graph=None, image=True, radii_factor=0.5):
        svg_string = conversion_tools.make_svg(
            self,
            filter_,
            species_type=species_type,
            graph=graph,
            image=image,
            radii_factor=radii_factor
        )
        with open(savefile[0], mode='w', newline='') as f:
            for line in svg_string.splitlines(keepends=True):
                f.write(line)

    @ staticmethod
    def import_from_file(filename, file_type, gui=None):
        if file_type == 'dm3':
            try:
                project = Project(filename, debug_obj=gui, species_dict=Project.default_species_dict)
            except:
                logger.info('Could not open .dm3 file.')
                project = None
            return project
        elif file_type == 'AtoMap':
            project = conversion_tools.import_from_atomap(filename, debug_obj=gui)
            return project
        else:
            logger.info('Unknown import format')
            return None


