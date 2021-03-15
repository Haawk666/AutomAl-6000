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
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPM
import logging
# Instantiate logger:
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def find_edge_columns(graph_obj, im_width, im_height):
    for vertex in graph_obj.vertices:
        x_coor = vertex.im_coor_x
        y_coor = vertex.im_coor_y
        margin = 3 * vertex.r
        if x_coor < margin or x_coor > im_width - margin - 1 or y_coor < margin or y_coor > im_height - margin - 1:
            graph_obj.vertices[vertex.i].is_edge_column = True
        else:
            graph_obj.vertices[vertex.i].is_edge_column = False


def calc_separation_matrix(graph_obj):
    projected_separation_matrix = np.zeros([graph_obj.order, graph_obj.order], dtype=float)
    for i in range(0, len(graph_obj.vertices) - 1):
        for j in range(i + 1, len(graph_obj.vertices)):
            dist = graph_obj.get_projected_separation(i, j)
            projected_separation_matrix[j, i] = dist
            projected_separation_matrix[i, j] = dist
    return projected_separation_matrix


def determine_districts(graph_obj, search_list=None):
    if graph_obj.projected_separation_matrix is None:
        graph_obj.projected_separation_matrix = calc_separation_matrix(graph_obj)

    if not graph_obj.projected_separation_matrix.shape == (len(graph_obj.vertices), len(graph_obj.vertices)):
        logger.warning('The projected separation matrix dimensions does not match the graph order.')

    if search_list is None:
        search_list = graph_obj.vertices
    else:
        search_list = graph_obj.get_vertex_objects_from_indices(search_list)
    for vertex in search_list:
        vertex.projected_separation_district = np.argsort(graph_obj.projected_separation_matrix[vertex.i, :])[1:graph_obj.district_size + 1].tolist()
        vertex.district = copy.deepcopy(vertex.projected_separation_district)
        vertex.local_zeta_map['district'] = copy.deepcopy(vertex.projected_separation_district)


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


def zeta_analysis(graph_obj, starting_index, print_states=False):
    # First determine the matrix zeta:
    votes = [0.0] * graph_obj.order
    votes[starting_index] = 1.0
    cont = True
    counter = 0
    while cont:
        new_votes = copy.deepcopy(votes)
        for vertex in graph_obj.vertices:
            if not vertex.is_in_precipitate:
                _sum = 0
                for partner in vertex.partners:
                    if graph_obj.vertices[partner].is_edge_column:
                        _sum += 0.1 * votes[partner]
                    else:
                        _sum += 0.5 * votes[partner]
                for in_semi_partner in vertex.in_semi_partners:
                    if graph_obj.vertices[in_semi_partner].is_edge_column:
                        _sum += 0.1 * votes[in_semi_partner]
                    else:
                        _sum += 0.2 * votes[in_semi_partner]
                new_votes[vertex.i] = min(100, max(-100, votes[vertex.i] - _sum))
        if print_states:
            print_state(graph_obj, new_votes, True, counter, starting_index)
        counter += 1
        votes = copy.deepcopy(new_votes)
        if counter > int(0.05 * graph_obj.order):
            cont = False
    for vertex in graph_obj.vertices:
        if not vertex.is_in_precipitate:
            if votes[vertex.i] > 0:
                votes[vertex.i] = 100
            else:
                votes[vertex.i] = -100
        else:
            votes[vertex.i] = 0.0

    # Then determine precipitate zeta:
    cont = True
    counter = 1
    while cont:
        new_votes = copy.deepcopy(votes)
        for vertex in graph_obj.vertices:
            if vertex.is_in_precipitate:
                _sum = 0
                for partner in vertex.partners:
                    if graph_obj.vertices[partner].is_edge_column:
                        _sum += 0.1 * votes[partner]
                    else:
                        _sum += 0.5 * votes[partner]
                for in_semi_partner in vertex.in_semi_partners:
                    if graph_obj.vertices[in_semi_partner].is_edge_column:
                        _sum += 0.1 * votes[in_semi_partner]
                    else:
                        _sum += 0.2 * votes[in_semi_partner]
                new_votes[vertex.i] = min(100, max(-100, votes[vertex.i] - _sum))
        if print_states:
            print_state(graph_obj, new_votes, False, counter, starting_index)
        counter += 1
        votes = copy.deepcopy(new_votes)
        if counter > int(0.05 * graph_obj.order):
            cont = False
    for vertex in graph_obj.vertices:
        if votes[vertex.i] > 0:
            votes[vertex.i] = 100
            vertex.set_zeta(0)
        else:
            votes[vertex.i] = -100
            vertex.set_zeta(1)

    graph_obj.build_local_zeta_maps()
    graph_obj.build_local_maps()

    if print_states:
        print_state(graph_obj, votes, False, counter, starting_index, complete=True)


def print_state(graph_obj, votes, matrix, cycle, seed_index, complete=False):

    max_x = 0
    min_x = 0
    max_y = 0
    min_y = 0

    for vertex in graph_obj.vertices:
        if vertex.im_coor_x > max_x:
            max_x = vertex.im_coor_x
        if vertex.im_coor_x < min_x:
            min_x = vertex.im_coor_x
        if vertex.im_coor_y > max_y:
            max_y = vertex.im_coor_y
        if vertex.im_coor_y < min_y:
            min_y = vertex.im_coor_y

    xml_string = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
    xml_string += '<!-- Created with AutomAl 6000 (http://www.automal.org/) -->\n\n'

    xml_string += '<svg\n   xmlns:xlink="http://www.w3.org/1999/xlink"\n'
    xml_string += '   height="{}"\n'.format(max_y - min_y)
    xml_string += '   width="{}">\n'.format(max_x - min_x)

    xml_string += '  <g\n'
    xml_string += '     inkscape:groupmode="layer"\n'
    xml_string += '     inkscape:label="Atomic graph [arcs]"\n'
    xml_string += '     id="Atomic graph [arcs]"\n'
    xml_string += '     style="display:inline" >\n'

    graph_obj.map_arcs()

    for arc in graph_obj.arcs:

        if arc.dual_arc:

            xml_string += '    <path\n'
            xml_string += '       id="arc{}"\n'.format(arc.j)
            xml_string += '       d="M {},{} {},{}"\n'.format(arc.vertex_a.im_coor_x, arc.vertex_a.im_coor_y,
                                                              arc.vertex_b.im_coor_x, arc.vertex_b.im_coor_y)

            if arc.co_planar:

                xml_string += '       stroke="{}"\n'.format(utils.rgb_to_hex((20, 20, 220)))
                xml_string += '       stroke-width="4" />\n'

            else:

                xml_string += '       stroke="{}"\n'.format(utils.rgb_to_hex((0, 0, 0)))
                xml_string += '       stroke-width="1" />\n'

        else:

            xml_string += '    <path\n'
            xml_string += '       id="arc{}"\n'.format(arc.j)
            xml_string += '       d="M {},{} {},{}"\n'.format(arc.vertex_a.im_coor_x, arc.vertex_a.im_coor_y,
                                                              arc.vertex_b.im_coor_x, arc.vertex_b.im_coor_y)
            xml_string += '       stroke="{}"\n'.format(utils.rgb_to_hex((220, 20, 20)))
            xml_string += '       stroke-width="4" />\n'

    xml_string += '  </g>\n'

    xml_string += '  <g\n'
    xml_string += '     inkscape:groupmode="layer"\n'
    xml_string += '     inkscape:label="Atomic graph [vertices]"\n'
    xml_string += '     id="Atomic graph [vertices]"\n'
    xml_string += '     style="display:inline" >\n'

    for vertex in graph_obj.vertices:

        grey_code = int((255 / 200) * (votes[vertex.i] + 100))
        if vertex.i == seed_index:
            boarder_colour = 220
        else:
            boarder_colour = 0

        xml_string += '    <circle\n'
        xml_string += '       id="graph_vertex_{}"\n'.format(vertex.i)
        xml_string += '       cx="{}"\n'.format(vertex.im_coor_x)
        xml_string += '       cy="{}"\n'.format(vertex.im_coor_y)
        xml_string += '       r="{}"\n'.format(vertex.r / 2 - 8 / graph_obj.scale)
        xml_string += '       style="fill:{};fill-opacity:1;stroke:{};stroke-width:{};stroke-opacity:1"\n'.format(
            utils.rgb_to_hex((grey_code, grey_code, grey_code)),
            utils.rgb_to_hex((0, boarder_colour, 0)),
            15 / graph_obj.scale
        )

    xml_string += '  </g>\n'

    xml_string += '</svg>\n'

    with open('temp.svg'.format(cycle), mode='w', newline='') as f:
        for line in xml_string.splitlines(keepends=True):
            f.write(line)
    drawing = svg2rlg('temp.svg')

    if complete:
        renderPM.drawToFile(drawing, 'snapshot_result.png', fmt='PNG')
    else:
        if matrix:
            renderPM.drawToFile(drawing, 'snapshot_matrix_{}.png'.format(cycle), fmt='PNG')
        else:
            renderPM.drawToFile(drawing, 'snapshot_precipitate_{}.png'.format(cycle), fmt='PNG')


def arc_intersection_denial(graph_obj):
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
    graph_obj.build_local_zeta_maps()


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
    if not strong:
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
    else:
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


def zeta_untangle(graph_obj):
    for k in [0, 1]:
        zeta_intersections = graph_obj.find_zeta_intersections()
        for intersection in zeta_intersections:
            vertex_a = graph_obj.vertices[intersection[0]]
            vertex_b = graph_obj.vertices[intersection[1]]
            vertex_c = graph_obj.vertices[intersection[2]]
            vertex_d = graph_obj.vertices[intersection[3]]
            if vertex_b.i in vertex_a.local_zeta_map['out_neighbourhood'] and vertex_a.i not in vertex_b.local_zeta_map['out_neighbourhood']:
                # type 1
                if vertex_d.i in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i not in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant a
                    if k > 0:
                        if graph_obj.get_projected_separation(vertex_a.i, vertex_b.i) < graph_obj.get_projected_separation(vertex_c.i, vertex_d.i):
                            if not vertex_c.is_edge_column:
                                vertex_c.decrement_n()
                        else:
                            if not vertex_a.is_edge_column:
                                vertex_a.decrement_n()
                elif vertex_d.i not in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant b
                    if k > 0:
                        if graph_obj.get_projected_separation(vertex_a.i, vertex_b.i) < graph_obj.get_projected_separation(vertex_c.i, vertex_d.i):
                            if not vertex_d.is_edge_column:
                                vertex_d.decrement_n()
                        else:
                            if not vertex_a.is_edge_column:
                                vertex_a.decrement_n()
                elif vertex_d.i in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant c
                    if not vertex_a.is_edge_column:
                        vertex_a.decrement_n()
            elif vertex_b.i not in vertex_a.local_zeta_map['out_neighbourhood'] and vertex_a.i in vertex_b.local_zeta_map['out_neighbourhood']:
                # type 2
                if vertex_d.i in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i not in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant a
                    if k > 0:
                        if graph_obj.get_projected_separation(vertex_a.i, vertex_b.i) < graph_obj.get_projected_separation(vertex_c.i, vertex_d.i):
                            if not vertex_c.is_edge_column:
                                vertex_c.decrement_n()
                        else:
                            if not vertex_b.is_edge_column:
                                vertex_b.decrement_n()
                elif vertex_d.i not in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant b
                    if k > 0:
                        if graph_obj.get_projected_separation(vertex_a.i, vertex_b.i) < graph_obj.get_projected_separation(vertex_c.i, vertex_d.i):
                            if not vertex_d.is_edge_column:
                                vertex_d.decrement_n()
                        else:
                            if not vertex_b.is_edge_column:
                                vertex_b.decrement_n()
                elif vertex_d.i in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant c
                    if not vertex_b.is_edge_column:
                        vertex_b.decrement_n()
            elif vertex_b.i in vertex_a.local_zeta_map['out_neighbourhood'] and vertex_a.i in vertex_b.local_zeta_map['out_neighbourhood']:
                # type 3
                if vertex_d.i in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i not in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant a
                    if not vertex_c.is_edge_column:
                        vertex_c.decrement_n()
                elif vertex_d.i not in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant b
                    if not vertex_d.is_edge_column:
                        vertex_d.decrement_n()
                elif vertex_d.i in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant c
                    if graph_obj.get_projected_separation(vertex_a.i, vertex_b.i) < graph_obj.get_projected_separation(vertex_c.i, vertex_d.i):
                        if not vertex_c.is_edge_column:
                            vertex_c.decrement_n()
                        if not vertex_d.is_edge_column:
                            vertex_d.decrement_n()
                    else:
                        if not vertex_a.is_edge_column:
                            vertex_a.decrement_n()
                        if not vertex_b.is_edge_column:
                            vertex_b.decrement_n()
        graph_obj.build_local_zeta_maps()
    for k in [0, 1]:
        zeta_intersections = graph_obj.find_anti_zeta_intersections()
        for intersection in zeta_intersections:
            print(intersection)
            vertex_a = graph_obj.vertices[intersection[0]]
            vertex_b = graph_obj.vertices[intersection[1]]
            vertex_c = graph_obj.vertices[intersection[2]]
            vertex_d = graph_obj.vertices[intersection[3]]
            if vertex_b.i in vertex_a.local_zeta_map['anti_out_neighbourhood'] and vertex_a.i not in vertex_b.local_zeta_map['anti_out_neighbourhood']:
                # type 1
                if vertex_d.i in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i not in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant a
                    pass
                elif vertex_d.i not in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant b
                    pass
                elif vertex_d.i in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant c
                    for index in vertex_a.local_zeta_map['district']:
                        if index in vertex_a.local_zeta_map['anti_in_semi_partners']:
                            j = vertex_b.i
                            k = index
                            if not j == k:
                                pos_j = -1
                                pos_k = -1
                                if j in vertex_a.local_zeta_map['district']:
                                    pos_j = vertex_a.local_zeta_map['district'].index(j)
                                if k in vertex_a.local_zeta_map['district']:
                                    pos_k = vertex_a.local_zeta_map['district'].index(k)

                                if pos_j == -1 and pos_k == -1:
                                    vertex_a.local_zeta_map['district'][-1] = k
                                    break
                                elif not pos_j == -1 and not pos_k == -1:
                                    vertex_a.local_zeta_map['district'][pos_j], vertex_a.local_zeta_map['district'][pos_k] = vertex_a.local_zeta_map['district'][pos_k], vertex_a.local_zeta_map['district'][pos_j]
                                    break
                                elif pos_j == -1:
                                    break
                                else:
                                    vertex_a.local_zeta_map['district'][-1] = k
                                    vertex_a.local_zeta_map['district'][pos_j], vertex_a.local_zeta_map['district'][pos_k] = vertex_a.local_zeta_map['district'][pos_k], vertex_a.local_zeta_map['district'][pos_j]
                                    break

            elif vertex_b.i not in vertex_a.local_zeta_map['anti_out_neighbourhood'] and vertex_a.i in vertex_b.local_zeta_map['anti_out_neighbourhood']:
                # type 2
                if vertex_d.i in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i not in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant a
                    pass
                elif vertex_d.i not in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant b
                    pass
                elif vertex_d.i in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant c
                    for index in vertex_b.local_zeta_map['district']:
                        if index in vertex_b.local_zeta_map['anti_in_semi_partners']:
                            j = vertex_a.i
                            k = index
                            if not j == k:
                                pos_j = -1
                                pos_k = -1
                                if j in vertex_a.local_zeta_map['district']:
                                    pos_j = vertex_a.local_zeta_map['district'].index(j)
                                if k in vertex_a.local_zeta_map['district']:
                                    pos_k = vertex_a.local_zeta_map['district'].index(k)

                                if pos_j == -1 and pos_k == -1:
                                    vertex_a.local_zeta_map['district'][-1] = k
                                    break
                                elif not pos_j == -1 and not pos_k == -1:
                                    vertex_a.local_zeta_map['district'][pos_j], vertex_a.local_zeta_map['district'][pos_k] = vertex_a.local_zeta_map['district'][pos_k], vertex_a.local_zeta_map['district'][pos_j]
                                    break
                                elif pos_j == -1:
                                    break
                                else:
                                    vertex_a.local_zeta_map['district'][-1] = k
                                    vertex_a.local_zeta_map['district'][pos_j], vertex_a.local_zeta_map['district'][pos_k] = vertex_a.local_zeta_map['district'][pos_k], vertex_a.local_zeta_map['district'][pos_j]
                                    break
            elif vertex_b.i in vertex_a.local_zeta_map['anti_out_neighbourhood'] and vertex_a.i in vertex_b.local_zeta_map['anti_out_neighbourhood']:
                # type 3
                if vertex_d.i in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i not in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant a
                    pass
                elif vertex_d.i not in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant b
                    pass
                elif vertex_d.i in vertex_c.local_zeta_map['out_neighbourhood'] and vertex_c.i in vertex_d.local_zeta_map['out_neighbourhood']:
                    # variant c
                    pass
        graph_obj.build_local_zeta_maps()
    graph_obj.match_zeta_graph()
    graph_obj.build_local_maps()



