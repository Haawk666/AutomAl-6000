# By Haakon Tvedt @ NTNU
# Contributors:

# Internal imports:
import data_module
import utils
import untangling_2
# External imports:
import numpy as np
import time
import copy
import logging
# Instantiate logger:
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def find_edge_columns(graph_obj, im_width, im_height):
    logger.info('Detecting edge columns')
    time_1 = time.time()
    for vertex in graph_obj.vertices:
        x_coor = vertex.im_coor_x
        y_coor = vertex.im_coor_y
        margin = 3 * vertex.r
        if x_coor < margin or x_coor > im_width - margin - 1 or y_coor < margin or y_coor > im_height - margin - 1:
            graph_obj.vertices[vertex.i].is_edge_column = True
        else:
            graph_obj.vertices[vertex.i].is_edge_column = False
    time_2 = time.time()
    logger.info('Found edge columns in {} seconds'.format(time_2 - time_1))


def calc_separation_matrix(graph_obj):
    logger.info('Starting separation calculation.\n    This could take some time, especially if there are many columns')
    time_1 = time.time()
    projected_separation_matrix = np.zeros([graph_obj.order, graph_obj.order], dtype=float)
    for i in range(0, len(graph_obj.vertices) - 1):
        for j in range(i + 1, len(graph_obj.vertices)):
            dist = graph_obj.get_projected_separation(i, j)
            projected_separation_matrix[j, i] = dist
            projected_separation_matrix[i, j] = dist
    time_2 = time.time()
    logger.info('Distance calculations took {} seconds.\n'.format(time_2 - time_1))
    return projected_separation_matrix


def determine_districts(graph_obj, search_list=None):
    logger.info('Determining districts based on projected separations..')
    if graph_obj.projected_separation_matrix is None:
        logger.info('    The projected separation matrix was None, calculating matrix..')
        graph_obj.projected_separation_matrix = calc_separation_matrix(graph_obj)

    if not graph_obj.projected_separation_matrix.shape == (len(graph_obj.vertices), len(graph_obj.vertices)):
        logger.warning('    The projected separation matrix dimensions does not match the graph order.')

    time_1 = time.time()
    if search_list is None:
        search_list = graph_obj.vertices
    else:
        search_list = graph_obj.get_vertex_objects_from_indices(search_list)
    for vertex in search_list:
        vertex.projected_separation_district = np.argsort(graph_obj.projected_separation_matrix[vertex.i, :])[1:graph_obj.district_size + 1].tolist()
        vertex.district = copy.deepcopy(vertex.projected_separation_district)
    time_2 = time.time()
    logger.info('Sorting took {} seconds.'.format(time_2 - time_1))


def particle_detection(graph_obj):
    for vertex in graph_obj.vertices:
        if vertex.void:
            vertex.is_in_precipitate = False
        else:
            if vertex.atomic_species == 'Al':
                num_foreign_species = 0
                for neighbour in vertex.neighbourhood:
                    if not graph_obj.vertices[neighbour].atomic_species == 'Al':
                        num_foreign_species += 1
                    if num_foreign_species == 2:
                        vertex.is_in_precipitate = True
                        break
                else:
                    vertex.is_in_precipitate = False
            else:
                vertex.is_in_precipitate = True


def find_precipitate_edge_sequence(graph_obj):
    """Find a cycle that represents the boarder of the precipitate.

    This method assumes a correct precipitate identification and that the precipitate is fully contained by the image



    """
    # Locate a vertex that is on the edge of the precipitate. Definition:
    #     A precipitate vertex that has a non-precipitate vertex in its partner set.
    predecessor = -1
    sequence = set()
    for vertex in graph_obj.vertices:
        if vertex.is_in_precipitate:
            for partner in vertex.partners:
                if not partner.is_in_precipitate:
                    sequence.add(vertex.i)
                    predecessor = vertex.i
                    break
    # Order the sequence into a path:
    if not predecessor == -1:
        ordered_sequence = [predecessor]
        sequence.remove(predecessor)
        successor = -1
        for partner in graph_obj.vertices[predecessor].partners:
            if partner in sequence:
                successor = partner
                break
        if not successor == -1:
            ordered_sequence.append(successor)
            sequence.remove(successor)
            while len(sequence) > 0:
                for partner in graph_obj.vertices[successor].partners:
                    if partner in sequence:
                        successor = partner
                        break
                else:
                    logger.info('Could not find precipitate boarder. Current attempt resulted in {}. Exiting method.'.format(ordered_sequence))
                    return ordered_sequence
            ordered_sequence.append(ordered_sequence[0])
            return ordered_sequence
        else:
            return [predecessor]
    else:
        return []


def calculate_precipitate_packing(graph_obj):
    """Calculate the precipitate packing fraction, assuming the precipitate is correctly identified.

    This method is only well defined for images that fully enclose the precipitate.

    """
    # Get boarder:
    boarder = find_precipitate_edge_sequence(graph_obj)

    # Calculate precipitate area:


def calculate_precipitate_displacement(graph_obj):
    """Calculate the ration between the number of columns in the precipitate and the number of columns in the FCC it displaces.


    """
    # Calculate the number of columns in the precipitate:
    num_precipitate_columns = 0
    for vertex in graph_obj.vertices:
        if vertex.is_in_precipitate:
            num_precipitate_columns += 1
    # Find a precipitate bounding rect:
    left_x = -1
    right_x = -1
    upper_y = -1
    lower_y = -1
    for vertex in graph_obj.vertices:
        if vertex.is_in_precipitate:
            if vertex.spatial_coor_x > right_x or right_x == -1:
                right_x = vertex.spatial_coor_x
            if vertex.spatial_coor_x < left_x or left_x == -1:
                left_x = vertex.spatial_coor_x
            if vertex.spatial_coor_y > upper_y or upper_y == -1:
                upper_y = vertex.spatial_coor_y
            if vertex.spatial_coor_y < lower_y or lower_y == -1:
                lower_y = vertex.spatial_coor_y


def zeta_analysis(graph_obj, starting_index, starting_zeta=0, method='district', use_n=False):
    logger.info('Starting zeta analysis')
    if method == 'weighted':
        logger.info('    Mapping cities..')
        for vertex in graph_obj.vertices:
            vertex.city = set()
            vertex.city.add(vertex.i)
            for citizen in vertex.district:
                vertex.city.add(citizen)
                for extended_citizen in graph_obj.vertices[citizen].district:
                    vertex.city.add(extended_citizen)
        graph_obj.build_local_maps(build_out=True)
        logger.info('    Cities mapped.')
    time_1 = time.time()
    votes = [0] * graph_obj.order
    if starting_zeta == 0:
        votes[starting_index] = 1
    else:
        votes[starting_index] = -1
    counter = 0
    cont = True
    while cont:
        for vertex in graph_obj.vertices:
            if use_n:
                n = vertex.n
            else:
                n = 3
            if method == 'district':
                candidates = vertex.district[0:n]
            elif method == 'separation':
                candidates = vertex.projected_separation_district[0:n]
            elif method == 'partners':
                candidates = list(vertex.partners)
            elif method == 'out_neighbours':
                candidates = list(vertex.out_neighbourhood)
            elif method == 'weighted':
                candidates = vertex.city
            elif method == 'binary':
                candidates = vertex.district
                if vertex.is_edge_column:
                    candidates = vertex.district[0:n]
            else:
                candidates = vertex.district[0:n]
            if method == 'binary':
                for citizen in candidates:
                    if citizen in vertex.partners:
                        votes[citizen] -= 0.5 * votes[vertex.i]
                    elif citizen in vertex.anti_neighbourhood:
                        votes[citizen] += 0.2 * votes[vertex.i]
                    if votes[citizen] > 100:
                        votes[citizen] = 100
                    elif votes[citizen] < -100:
                        votes[citizen] = -100
            else:
                if vertex.is_edge_column:
                    for citizen in candidates:
                        votes[citizen] -= 0.01 * votes[vertex.i]
                        if votes[citizen] > 100:
                            votes[citizen] = 100
                        elif votes[citizen] < -100:
                            votes[citizen] = -100
                else:
                    for citizen in candidates:
                        r = graph_obj.projected_separation_matrix[vertex.i, citizen]
                        factor = 1 - r / 404.95
                        votes[citizen] -= factor * votes[vertex.i]
                        if votes[citizen] > 100:
                            votes[citizen] = 100
                        elif votes[citizen] < -100:
                            votes[citizen] = -100
        counter += 1
        if counter > 1000:
            cont = False

    for vertex in graph_obj.vertices:
        if votes[vertex.i] > 0:
            vertex.set_zeta(0)
        else:
            vertex.set_zeta(1)

    graph_obj.build_local_zeta_maps()
    # graph_obj.build_local_maps()

    time_2 = time.time()
    logger.info('Zeta analysis completed in {} seconds'.format(time_2 - time_1))


def arc_intersection_denial(graph_obj):
    time_1 = time.time()
    intersections = graph_obj.find_intersections()
    for intersection in intersections:
        if not graph_obj.vertices[intersection[0]].partner_query(intersection[1]):
            if not graph_obj.terminate_arc(intersection[0], intersection[1]):
                if not graph_obj.terminate_arc(intersection[2], intersection[3]):
                    logger.info('Could not remove intersection {}'.format(intersection))
        elif not graph_obj.vertices[intersection[2]].partner_query(intersection[3]):
            if not graph_obj.terminate_arc(intersection[2], intersection[3]):
                if not graph_obj.terminate_arc(intersection[0], intersection[1]):
                    logger.info('Could not remove intersection {}'.format(intersection))
        else:
            if not graph_obj.terminate_arc(intersection[2], intersection[3]):
                if not graph_obj.terminate_arc(intersection[0], intersection[1]):
                    logger.info('Could not remove intersection {}'.format(intersection))
    graph_obj.build_local_maps()
    intersections = graph_obj.find_intersections()
    for intersection in intersections:
        if not graph_obj.vertices[intersection[0]].partner_query(intersection[1]):
            if not graph_obj.terminate_arc(intersection[0], intersection[1]):
                if not graph_obj.terminate_arc(intersection[2], intersection[3]):
                    logger.info('Could not remove intersection {}'.format(intersection))
        elif not graph_obj.vertices[intersection[2]].partner_query(intersection[3]):
            if not graph_obj.terminate_arc(intersection[2], intersection[3]):
                if not graph_obj.terminate_arc(intersection[0], intersection[1]):
                    logger.info('Could not remove intersection {}'.format(intersection))
        else:
            if not graph_obj.terminate_arc(intersection[2], intersection[3]):
                if not graph_obj.terminate_arc(intersection[0], intersection[1]):
                    logger.info('Could not remove intersection {}'.format(intersection))
    graph_obj.build_local_maps()
    time_3 = time.time()
    logger.info('Performed arc intersection denial in {} seconds'.format(time_3 - time_1))


def apply_alpha_model(graph_obj, model=None, alpha_selection_type='zeta'):
    if model is None:
        this_model = data_module.VertexDataManager.load(graph_obj.active_model)
    else:
        this_model = model
    for vertex in graph_obj.vertices:
        vertex.alpha_angles = graph_obj.get_alpha_angles(vertex.i, selection_type=alpha_selection_type)
        if vertex.alpha_angles is not None and not len(vertex.alpha_angles) == 0:
            vertex.alpha_max = max(vertex.alpha_angles)
            vertex.alpha_min = min(vertex.alpha_angles)
        else:
            vertex.alpha_max = 0
            vertex.alpha_min = 0
    for vertex in graph_obj.vertices:
        if not vertex.is_edge_column and not vertex.void and not vertex.is_set_by_user:
            vertex.advanced_probability_vector = this_model.calc_prediction(
                {
                    'alpha_max': vertex.alpha_max,
                    'alpha_min': vertex.alpha_min
                },
                graph_obj.get_advanced_species_list()
            )
            vertex.advanced_probability_vector['Un_1'] = 0.0
            vertex.determine_species_from_probability_vector()
    graph_obj.build_local_maps()


def apply_composite_model(graph_obj, model=None, alpha_selection_type='zeta'):
    if model is None:
        this_model = data_module.VertexDataManager.load(graph_obj.active_model)
    else:
        this_model = model
    for vertex in graph_obj.vertices:
        vertex.alpha_angles = graph_obj.get_alpha_angles(vertex.i, selection_type=alpha_selection_type)
        if vertex.alpha_angles is not None and not len(vertex.alpha_angles) == 0:
            vertex.alpha_max = max(vertex.alpha_angles)
            vertex.alpha_min = min(vertex.alpha_angles)
        else:
            vertex.alpha_max = 0
            vertex.alpha_min = 0
        vertex.theta_angles = graph_obj.get_theta_angles(vertex.i, selection_type='selective')
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
    for vertex in graph_obj.vertices:
        if not vertex.is_edge_column and not vertex.void and not vertex.is_set_by_user:
            dict_ = {
                'alpha_max': vertex.alpha_max,
                'alpha_min': vertex.alpha_min,
                'theta_angle_mean': vertex.theta_angle_mean,
                'normalized_peak_gamma': vertex.normalized_peak_gamma,
                'normalized_avg_gamma': vertex.normalized_peak_gamma
            }
            keys_to_pop = []
            for key, value in dict_.items():
                if value == 0:
                    keys_to_pop.append(key)
            for key in keys_to_pop:
                dict_.pop(key)
            vertex.advanced_probability_vector = this_model.calc_prediction(
                dict_,
                graph_obj.get_advanced_species_list()
            )
            vertex.advanced_probability_vector['Un_1'] = 0.0
            vertex.determine_species_from_probability_vector()
    graph_obj.build_local_maps()
    graph_obj.build_local_zeta_maps()


def untangle(graph_obj, ui_obj=None, strong=True):

    # Untangling:
    for vertex in graph_obj.vertices:
        if not vertex.is_edge_column and not vertex.void:
            if len(vertex.out_semi_partners) > 0:
                district_copy = copy.deepcopy(vertex.district)
                for citizen in district_copy:
                    if citizen in vertex.out_semi_partners:
                        sub_graph = graph_obj.get_arc_centered_subgraph(vertex.i, citizen)
                        untangling_2.determine_sub_graph_class([sub_graph])
                        untangling_2.determine_sub_graph_configuration(graph_obj, [sub_graph])
                        untangling_2.weak_resolve(graph_obj, [sub_graph], ui_obj=ui_obj)
    if strong:
        for vertex in graph_obj.vertices:
            if not vertex.is_edge_column and not vertex.void:
                if len(vertex.out_semi_partners) > 0:
                    district_copy = copy.deepcopy(vertex.district)
                    for citizen in district_copy:
                        if citizen in vertex.out_semi_partners:
                            sub_graph = graph_obj.get_arc_centered_subgraph(vertex.i, citizen)
                            untangling_2.determine_sub_graph_class([sub_graph])
                            untangling_2.determine_sub_graph_configuration(graph_obj, [sub_graph])
                            untangling_2.strong_resolve(graph_obj, [sub_graph], ui_obj=ui_obj)





