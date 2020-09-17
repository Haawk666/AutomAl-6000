from copy import deepcopy
import numpy as np
import utils


def base_angle_score(graph_obj, i, apply=True):

    alpha = graph_obj.get_alpha_angles(i)

    if apply:

        cu_min_mean = 1.92
        si_1_min_mean = 1.94
        si_2_min_mean = 1.56
        al_min_mean = 1.56
        mg_1_min_mean = 1.30
        mg_2_min_mean = 1.26

        cu_min_std = 0.19
        si_1_min_std = 0.14
        si_2_min_std = 0.08
        al_min_std = 0.05
        mg_1_min_std = 0.09
        mg_2_min_std = 0.05

        cu_max_mean = 2.28
        si_1_max_mean = 2.25
        si_2_max_mean = 2.48
        al_max_mean = 3.11
        mg_1_max_mean = 2.57
        mg_2_max_mean = 3.69

        cu_max_std = 0.26
        si_1_max_std = 0.12
        si_2_max_std = 0.15
        al_max_std = 0.07
        mg_1_max_std = 0.10
        mg_2_max_std = 0.09

        cf_cu_min = utils.normal_dist(min(alpha), cu_min_mean, cu_min_std)
        cf_si_1_min = utils.normal_dist(min(alpha), si_1_min_mean, si_1_min_std)
        cf_si_2_min = utils.normal_dist(min(alpha), si_2_min_mean, si_2_min_std)
        cf_al_min = utils.normal_dist(min(alpha), al_min_mean, al_min_std)
        cf_mg_1_min = utils.normal_dist(min(alpha), mg_1_min_mean, mg_1_min_std)
        cf_mg_2_min = utils.normal_dist(min(alpha), mg_2_min_mean, mg_2_min_std)

        cf_cu_max = utils.normal_dist(max(alpha), cu_max_mean, cu_max_std)
        cf_si_1_max = utils.normal_dist(max(alpha), si_1_max_mean, si_1_max_std)
        cf_si_2_max = utils.normal_dist(max(alpha), si_2_max_mean, si_2_max_std)
        cf_al_max = utils.normal_dist(max(alpha), al_max_mean, al_max_std)
        cf_mg_1_max = utils.normal_dist(max(alpha), mg_1_max_mean, mg_1_max_std)
        cf_mg_2_max = utils.normal_dist(max(alpha), mg_2_max_mean, mg_2_max_std)

        cf_min = [cf_cu_min, cf_si_1_min, cf_si_2_min, cf_al_min, cf_mg_1_min, cf_mg_2_min]
        cf_max = [cf_cu_max, cf_si_1_max, cf_si_2_max, cf_al_max, cf_mg_1_max, cf_mg_2_max]

        cf = [a * b for a, b in zip(cf_min, cf_max)]
        probs = utils.normalize_list(cf)
        sum_probs = [probs[1] + probs[2], probs[0], probs[3], probs[4] + probs[5], 0]
        advanced_probability_vector = {
            'Si_1': cf_si_1_max * cf_si_1_min,
            'Si_2': cf_si_2_max * cf_si_2_min,
            'Cu_1': cf_cu_max * cf_cu_min,
            'Al_1': cf_al_max * cf_al_min / 2,
            'Al_2': cf_al_max * cf_al_min / 2,
            'Mg_1': cf_mg_1_max * cf_mg_1_min,
            'Mg_2': cf_mg_2_max * cf_mg_2_min,
            'Un_1': 0
        }
        advanced_probability_vector = utils.normalize_dict(advanced_probability_vector, norm_sum=1)

        return advanced_probability_vector

    else:

        return max(alpha), min(alpha)


def precipitate_controller(graph, i):

    graph.particle_boarder_indices = []

    graph.reset_all_flags()

    precipitate_finder(graph, i)

    counter = 0

    for x in range(0, graph.order):

        if graph.vertices[x].flag_1 or graph.vertices[x].atomic_species == 'Un':
            graph.vertices[x].is_in_precipitate = False
        else:
            graph.vertices[x].is_in_precipitate = True

        if graph.vertices[x].flag_2:
            graph.particle_boarder_indices.append(x)
            counter = counter + 1

    graph.reset_all_flags()
    sort_boarder(graph)


def sort_boarder(graph):

    temp_boarder = deepcopy(graph.particle_boarder_indices)
    if len(temp_boarder) > 1:
        selected = []
        for y in range(0, len(graph.particle_boarder_indices)):
            selected.append(False)
        next_index = 0
        index = 0
        cont_var = True
        selected[0] = True

        while cont_var:

            distance = 1000000

            for x in range(0, len(graph.particle_boarder_indices)):

                current_distance = graph.get_projected_image_separation(graph.particle_boarder_indices[x], temp_boarder[index])

                if current_distance < distance and not temp_boarder[index] == graph.particle_boarder_indices[x] and not selected[x]:
                    distance = current_distance
                    next_index = x

            selected[next_index] = True
            index = index + 1

            temp_boarder[index] = graph.particle_boarder_indices[next_index]

            if index == len(graph.particle_boarder_indices) - 1:
                cont_var = False

            graph.particle_boarder_indices = deepcopy(temp_boarder)


def precipitate_finder(graph, i):

    graph.vertices[i].flag_1 = True

    for partner in graph.vertices[i].partners:

        if not graph.vertices[partner].atomic_species == 'Al':

            if not graph.vertices[partner].atomic_species == 'Un':
                graph.vertices[i].flag_2 = True

        else:

            if not graph.vertices[partner].flag_1:
                precipitate_finder(graph, partner)


def define_levels(graph, i, level=0):
    graph.reset_all_flags()

    mesh_levels(graph, i, level)

    complete = False
    emer_abort = False
    overcounter = 0
    neighbour_level = 0

    while not complete and not emer_abort:

        found = False
        counter = 0

        while counter < graph.order:

            if graph.vertices[counter].is_in_precipitate and not graph.vertices[counter].flag_1:

                x = 0

                while x <= neighbour_level:

                    if graph.vertices[graph.vertices[counter].district[x]].is_in_precipitate and \
                            graph.vertices[graph.vertices[counter].district[x]].flag_1:

                        neighbour = graph.vertices[counter].district[x]
                        if graph.vertices[neighbour].zeta == 0:
                            graph.vertices[counter].zeta = 1
                        else:
                            graph.vertices[counter].zeta = 0
                        graph.vertices[counter].flag_1 = True
                        found = True

                    x = x + 1

            counter = counter + 1

        complete = True

        for y in range(0, graph.order):

            if graph.vertices[y].is_in_precipitate and not graph.vertices[y].flag_1:

                complete = False

        if found and neighbour_level > 0:

            neighbour_level = neighbour_level - 1

        if not found and neighbour_level < 2:

            neighbour_level = neighbour_level + 1

        overcounter = overcounter + 1
        if overcounter > 100:

            emer_abort = True

    graph.reset_all_flags()


def mesh_levels(graph, i, level):

    if graph.vertices[i].is_in_precipitate:

        graph.vertices[i].flag_1 = True
        set_level(graph, i, level)

    else:

        graph.vertices[i].flag_1 = True

        next_level = 0
        if level == 0:
            next_level = 1
        elif level == 1:
            next_level = 0

        set_level(graph, i, level)

        indices = graph.vertices[i].district

        for x in range(0, graph.vertices[i].n):
            reciprocal = graph.vertices[indices[x]].partner_query(i)

            if not graph.vertices[indices[x]].flag_1 and not graph.vertices[i].is_edge_column and reciprocal:

                mesh_levels(graph, indices[x], next_level)


def precipitate_levels(graph, i, level):

    if not graph.vertices[i].is_in_precipitate:

        graph.vertices[i].flag_1 = True

    else:

        graph.vertices[i].flag_1 = True

        next_level = 0
        if level == 0:
            next_level = 1
        elif level == 1:
            next_level = 0

        set_level(graph, i, level)

        indices = graph.vertices[i].neighbour_indices

        complete = False
        counter_1 = 0
        counter_2 = 0

        while not complete:

            if not graph.vertices[indices[counter_1]].flag_1:

                if graph.test_reciprocality(i, indices[counter_1]):

                    precipitate_levels(indices[counter_1], next_level)
                    counter_1 = counter_1 + 1
                    counter_2 = counter_2 + 1

                else:

                    counter_1 = counter_1 + 1

            else:

                counter_1 = counter_1 + 1

            if counter_2 == graph.vertices[i].n() - 2 or counter_1 == graph.vertices[i].n() - 2:

                complete = True


def set_level(graph, i, level):

    previous_level = graph.vertices[i].zeta
    graph.vertices[i].zeta = level

    if level == previous_level:
        return False
    else:
        return True


def classify_pair(graph, i, j):

    neighbour_type = 0
    partner_type = 0
    intersects = False

    i_neighbour_to_j = False
    j_neighbour_to_i = False
    i_partner_to_j = False
    j_partner_to_i = False

    for x in range(0, 8):

        if graph.vertices[i].neighbour_indices[x] == j:
            j_neighbour_to_i = True
            if x < graph.vertices[i].n():
                j_partner_to_i = True

        if graph.vertices[j].neighbour_indices[x] == i:
            i_neighbour_to_j = True
            if x < graph.vertices[j].n():
                i_partner_to_j = True

    if not i_neighbour_to_j and not j_neighbour_to_i:
        neighbour_type = 0
    elif not i_neighbour_to_j and j_neighbour_to_i:
        neighbour_type = 1
    elif i_neighbour_to_j and not j_neighbour_to_i:
        neighbour_type = 2
    elif i_neighbour_to_j and j_neighbour_to_i:
        neighbour_type = 3

    if not i_partner_to_j and not j_partner_to_i:
        partner_type = 0
    elif not i_partner_to_j and j_partner_to_i:
        partner_type = 1
    elif i_partner_to_j and not j_partner_to_i:
        partner_type = 2
    elif i_partner_to_j and j_partner_to_i:
        partner_type = 3

    if graph.vertices[i].level == graph.vertices[j].level:
        level_type = 0
    else:
        level_type = 1

    if partner_type == 0:

        geometry_type_clockwise = -1
        geometry_type_anticlockwise = -1
        geo_type_symmetry = -1

    else:

        indices, num_edges_clockwise_right = find_shape(graph, i, j, clockwise=True)
        indices, num_edges_clockwise_left = find_shape(graph, j, i, clockwise=True)
        indices, num_edges_anticlockwise_right = find_shape(graph, j, i, clockwise=False)
        indices, num_edges_anticlockwise_left = find_shape(graph, i, j, clockwise=False)

        if num_edges_clockwise_right == 3 and num_edges_clockwise_left == 3:
            geometry_type_clockwise = 1
        elif num_edges_clockwise_right == 5 and num_edges_clockwise_left == 3:
            geometry_type_clockwise = 2
        elif num_edges_clockwise_right == 3 and num_edges_clockwise_left == 5:
            geometry_type_clockwise = 3
        elif num_edges_clockwise_right == 4 and num_edges_clockwise_left == 3:
            geometry_type_clockwise = 4
        elif num_edges_clockwise_right == 3 and num_edges_clockwise_left == 4:
            geometry_type_clockwise = 5
        elif num_edges_clockwise_right == 4 and num_edges_clockwise_left == 4:
            geometry_type_clockwise = 6
        elif num_edges_clockwise_right == 5 and num_edges_clockwise_left == 5:
            geometry_type_clockwise = 7
        else:
            geometry_type_clockwise = 0

        if num_edges_anticlockwise_right == 3 and num_edges_anticlockwise_left == 3:
            geometry_type_anticlockwise = 1
        elif num_edges_anticlockwise_right == 5 and num_edges_anticlockwise_left == 3:
            geometry_type_anticlockwise = 2
        elif num_edges_anticlockwise_right == 3 and num_edges_anticlockwise_left == 5:
            geometry_type_anticlockwise = 3
        elif num_edges_anticlockwise_right == 4 and num_edges_anticlockwise_left == 3:
            geometry_type_anticlockwise = 4
        elif num_edges_anticlockwise_right == 3 and num_edges_anticlockwise_left == 4:
            geometry_type_anticlockwise = 5
        elif num_edges_anticlockwise_right == 4 and num_edges_anticlockwise_left == 4:
            geometry_type_anticlockwise = 6
        elif num_edges_anticlockwise_right == 5 and num_edges_anticlockwise_left == 5:
            geometry_type_anticlockwise = 7
        else:
            geometry_type_anticlockwise = 0

        if geometry_type_clockwise == geometry_type_anticlockwise:
            geo_type_symmetry = 0
        else:
            geo_type_symmetry = 1

    # Implement method to find intersections

    return neighbour_type, partner_type, level_type, geometry_type_clockwise, geometry_type_anticlockwise,\
        geo_type_symmetry, intersects


def resolve_edge_inconsistency(graph, i, j, clockwise=True):

    neighbour_type, partner_type, level_type, geometry_type_clockwise, geometry_type_anticlockwise, \
        geo_type_symmetry, intersects = classify_pair(graph, i, j)

    if geo_type_symmetry == 0:
        geometry_type = geometry_type_clockwise
    else:
        geometry_type = 0

    if neighbour_type == 0:
        i_neighbour_to_j = False
        j_neighbour_to_i = False
    elif neighbour_type == 1:
        i_neighbour_to_j = False
        j_neighbour_to_i = True
    elif neighbour_type == 2:
        i_neighbour_to_j = True
        j_neighbour_to_i = False
    elif neighbour_type == 3:
        i_neighbour_to_j = True
        j_neighbour_to_i = True
    else:
        i_neighbour_to_j = False
        j_neighbour_to_i = False

    if partner_type == 0:
        i_partner_to_j = False
        j_partner_to_i = False
    elif partner_type == 1:
        i_partner_to_j = False
        j_partner_to_i = True
    elif partner_type == 2:
        i_partner_to_j = True
        j_partner_to_i = False
    elif partner_type == 3:
        i_partner_to_j = True
        j_partner_to_i = True
    else:
        i_partner_to_j = False
        j_partner_to_i = False

    i_index_in_j = -1
    j_index_in_i = -1

    for x in range(0, 8):
        if graph.vertices[i].neighbour_indices[x] == j:
            j_index_in_i = x

    if j_index_in_i == -1:
        graph.vertices[i].neighbour_indices[7] = j
        j_index_in_i = 7

    for x in range(0, 8):
        if graph.vertices[j].neighbour_indices[x] == i:
            i_index_in_j = x

    if i_index_in_j == -1:
        graph.vertices[j].neighbour_indices[7] = i
        i_index_in_j = 7

    if i_partner_to_j:
        # Perturb neighbours of j such that i is last element in k^j
        perturbator(graph, j, i_index_in_j, graph.vertices[j].n() - 1)
    else:
        # Perturb neighbours of j such that i is last element in k^j
        perturbator(graph, j, i_index_in_j, graph.vertices[j].n())

    if j_partner_to_i:
        # Perturb neighbours of i such that j is last k
        perturbator(graph, i, j_index_in_i, graph.vertices[i].n() - 1)
    else:
        # Perturb neighbours of i such that j is last k
        perturbator(graph, i, j_index_in_i, graph.vertices[i].n())

    if clockwise:

        shape_1_indices, num_edge_1 = find_shape(graph, i, j, clockwise=clockwise)
        shape_2_indices, num_edge_2 = find_shape(graph, j, i, clockwise=clockwise)

    else:

        shape_1_indices, num_edge_1 = find_shape(graph, j, i, clockwise=clockwise)
        shape_2_indices, num_edge_2 = find_shape(graph, i, j, clockwise=clockwise)

    if geometry_type == 1:
        # This means we want to break the connection!

        if partner_type == 1:
            if not try_connect(graph, i, j_index_in_i):
                if not graph.decrease_h_value(i):
                    print('Could not reconnect!')

        elif partner_type == 2:
            if not try_connect(graph, j, i_index_in_j):
                if not graph.decrease_h_value(j):
                    print('Could not reconnect!')

        elif partner_type == 3:
            if not try_connect(graph, j, i_index_in_j):
                if not graph.decrease_h_value(j):
                    print('Could not reconnect!')
            if not try_connect(graph, i, j_index_in_i):
                if not graph.decrease_h_value(i):
                    print('Could not reconnect!')

    if geometry_type == 2 or geometry_type == 3:
        # This means we want to switch connections to make geometry type 6

        loser_index_in_stayer = -1
        loser_connected_to_stayer = False
        stayer_connected_to_loser = False
        new_index_in_stayer = -1

        if geometry_type == 2:

            ind_1 = shape_1_indices[2]
            ind_2 = shape_1_indices[4]

        else:

            ind_1 = shape_2_indices[4]
            ind_2 = shape_2_indices[2]

        distance_1 = np.sqrt((graph.vertices[ind_1].x - graph.vertices[i].x) ** 2 + (
                graph.vertices[ind_1].y - graph.vertices[i].y) ** 2)
        distance_2 = np.sqrt((graph.vertices[j].x - graph.vertices[ind_2].x) ** 2 + (
                graph.vertices[j].y - graph.vertices[ind_2].y) ** 2)

        if distance_1 < distance_2:
            index_stayer = i
            index_loser = j
            stayer_index_in_loser = i_index_in_j
            index_new = ind_1
            if i_partner_to_j:
                loser_connected_to_stayer = True
            if j_partner_to_i:
                stayer_connected_to_loser = True
        else:
            index_stayer = j
            index_loser = i
            stayer_index_in_loser = j_index_in_i
            index_new = ind_2
            if j_partner_to_i:
                loser_connected_to_stayer = True
            if i_partner_to_j:
                stayer_connected_to_loser = True

        for x in range(graph.vertices[index_stayer].n(), 8):
            if graph.vertices[index_stayer].neighbour_indices[x] == index_new:
                new_index_in_stayer = x

        if new_index_in_stayer == -1:
            graph.vertices[index_stayer].neighbour_indices[7] = index_new
            new_index_in_stayer = 7

        perturbator(graph, index_stayer, graph.vertices[index_stayer].n(), new_index_in_stayer)

        if loser_connected_to_stayer:
            if not try_connect(graph, index_loser, stayer_index_in_loser):
                if not graph.decrease_h_value(index_loser):
                    print('Could not reconnect!')

        if stayer_connected_to_loser:
            perturbator(graph, index_stayer, loser_index_in_stayer, new_index_in_stayer)
        else:
            if not graph.increase_h_value(index_stayer):
                print('Could not reconnect!')

    if geometry_type == 4 or geometry_type == 5:

        pass

    if geometry_type == 6:
        # This means we want to keep the connection

        if partner_type == 1:
            if not graph.increase_h_value(j):
                print('Could not reconnect!')

        elif partner_type == 2:
            if not graph.increase_h_value(i):
                print('Could not reconnect!')

    if geometry_type == 0:

        print(str(num_edge_1) + ', ' + str(num_edge_2))
        print(shape_1_indices)
        print(shape_2_indices)


def try_connect(graph, i, j_index_in_i):

    changed = False
    better_friend = False
    friend_index_in_i = -1

    for x in range(graph.vertices[i].n(), 8):

        if graph.test_reciprocality(i, graph.vertices[i].neighbour_indices[x]):

            if not graph.vertices[i].level == graph.vertices[graph.vertices[i].neighbour_indices[x]].level:

                better_friend = True
                friend_index_in_i = x

            else:

                print('Maybe should have?')

    if better_friend:

        perturbator(graph, i, j_index_in_i, friend_index_in_i)
        changed = True

    return changed


def perturbator(graph, i, a, b):

    val_a = graph.vertices[i].neighbour_indices[a]
    val_b = graph.vertices[i].neighbour_indices[b]

    graph.vertices[i].neighbour_indices[a] = val_b
    graph.vertices[i].neighbour_indices[b] = val_a


def find_consistent_perturbations_simple(graph, y, sub=False):

    if not graph.vertices[y].is_edge_column:

        n = graph.vertices[y].n()

        if sub:

            n = 3

        indices = graph.vertices[y].neighbour_indices
        new_indices = np.zeros([len(indices)], dtype=int)
        found = 0

        for x in range(0, len(indices)):

            n2 = 3

            if graph.vertices[indices[x]].h_index == 0 or graph.vertices[indices[x]].h_index == 1:
                n2 = 3
            elif graph.vertices[indices[x]].h_index == 3:
                n2 = 4
            elif graph.vertices[indices[x]].h_index == 5:
                n2 = 5
            else:
                print('Problem in find_consistent_perturbations_simple!')

            neighbour_indices = graph.vertices[indices[x]].neighbour_indices

            for z in range(0, n2):

                if neighbour_indices[z] == y:
                    new_indices[found] = indices[x]
                    found = found + 1

        if found == n:

            index_positions = np.zeros([found], dtype=int)

            for k in range(0, found):

                for z in range(0, len(indices)):

                    if indices[z] == new_indices[k]:
                        index_positions[k] = z

            counter = found - 1

            for i in range(0, len(indices)):

                are_used = False

                for z in range(0, found):

                    if i == index_positions[z]:

                        are_used = True

                if not are_used:

                    counter = counter + 1
                    new_indices[counter] = indices[i]

            graph.vertices[y].neighbour_indices = new_indices
            graph.vertices[y].is_popular = False
            graph.vertices[y].is_unpopular = False

        elif found > n:

            # Here

            index_positions = np.zeros([found], dtype=int)

            for k in range(0, found):

                for z in range(0, len(indices)):

                    if int(indices[z]) == int(new_indices[k]):
                        index_positions[k] = z

            counter = found - 1

            print(index_positions)
            print(indices)
            print(new_indices)

            for i in range(0, len(indices)):

                are_used = False

                for z in range(0, found):

                    if i == index_positions[z]:
                        are_used = True

                if not are_used:
                    counter = counter + 1
                    new_indices[counter] = indices[i]

            graph.vertices[y].neighbour_indices = new_indices
            graph.vertices[y].is_unpopular = False
            graph.vertices[y].is_popular = True

        else:

            index_positions = np.zeros([found], dtype=int)

            for k in range(0, found):

                for z in range(0, len(indices)):

                    if indices[z] == new_indices[k]:
                        index_positions[k] = z

            counter = found - 1

            for i in range(0, len(indices)):

                are_used = False

                for z in range(0, found):

                    if i == index_positions[z]:
                        are_used = True

                if not are_used:
                    counter = counter + 1
                    new_indices[counter] = indices[i]

            graph.vertices[y].neighbour_indices = new_indices
            graph.vertices[y].is_unpopular = True
            graph.vertices[y].is_popular = False


def sort_neighbours_by_level(graph, y):

    n = graph.vertices[y].n()

    num_wrong_flags = 0

    for x in range(0, n):

        if graph.vertices[graph.vertices[y].neighbour_indices[x]].level == graph.vertices[y].level:
            num_wrong_flags = num_wrong_flags + 1

    if num_wrong_flags >= n - 1:

        if graph.vertices[y].level == 0:
            graph.vertices[y].level = 1
        else:
            graph.vertices[y].level = 0

        num_wrong_flags = 0

        for x in range(0, n):

            if graph.vertices[graph.vertices[y].neighbour_indices[x]].level == graph.vertices[y].level:
                num_wrong_flags = num_wrong_flags + 1

    finished = False
    debug_counter = 0

    while not finished:

        print(debug_counter)

        num_perturbations = 0

        for x in range(0, n - 1):

            if not graph.vertices[graph.vertices[y].neighbour_indices[x]].level == graph.vertices[y].level:
                pass
            else:
                perturbator(graph, y, x, x + 1)
                num_perturbations = num_perturbations + 1

            if x == n - num_wrong_flags - 2 and num_perturbations == 0:

                finished = True

        debug_counter = debug_counter + 1


def find_consistent_perturbations_advanced(graph, y, experimental=False):

    if not graph.vertices[y].is_edge_column:

        n = graph.vertices[y].n()

        indices = deepcopy(graph.vertices[y].neighbour_indices)
        new_indices = np.zeros([indices.shape[0]], dtype=int)
        index_of_unpopular_neighbours = np.zeros([indices.shape[0]], dtype=int)
        found = 0

        for x in range(0, len(indices)):

            if graph.vertices[indices[x]].is_unpopular:
                index_of_unpopular_neighbours[x] = indices[x]
            else:
                index_of_unpopular_neighbours[x] = -1

            n2 = 3

            if graph.vertices[indices[x]].h_index == 0 or graph.vertices[indices[x]].h_index == 1:
                n2 = 3
            elif graph.vertices[indices[x]].h_index == 3:
                n2 = 4
            elif graph.vertices[indices[x]].h_index == 5:
                n2 = 5
            else:
                print('Problem in find_consistent_perturbations_simple!')

            neighbour_indices = graph.vertices[indices[x]].neighbour_indices

            for z in range(0, n2):

                if neighbour_indices[z] == y:
                    new_indices[found] = indices[x]
                    found = found + 1
                    index_of_unpopular_neighbours[x] = -1

        if found == n:

            index_positions = np.zeros([found], dtype=int)

            for k in range(0, found):

                for z in range(0, len(indices)):

                    if indices[z] == new_indices[k]:
                        index_positions[k] = z

            counter = found - 1

            for i in range(0, len(indices)):

                are_used = False

                for z in range(0, found):

                    if i == index_positions[z]:

                        are_used = True

                if not are_used:

                    counter = counter + 1
                    new_indices[counter] = indices[i]

            graph.vertices[y].neighbour_indices = new_indices
            graph.vertices[y].is_unpopular = False
            graph.vertices[y].is_popular = False

            if experimental:
                sort_neighbours_by_level(graph, y)

        elif found > n:

            index_positions = np.zeros([found], dtype=int)

            for k in range(0, found):

                for z in range(0, len(indices)):

                    if indices[z] == new_indices[k]:
                        index_positions[k] = z

            counter = found - 1

            for i in range(0, len(indices)):

                are_used = False

                for z in range(0, found):

                    if i == index_positions[z]:
                        are_used = True

                if not are_used:
                    counter = counter + 1
                    new_indices[counter] = indices[i]

            graph.vertices[y].neighbour_indices = new_indices
            graph.vertices[y].is_popular = True
            graph.vertices[y].is_unpopular = False

            if experimental:
                sort_neighbours_by_level(graph, y)

        else:

            print(index_of_unpopular_neighbours)

            index_positions = np.zeros([found], dtype=int)

            for k in range(0, found):

                for z in range(0, len(indices)):

                    if indices[z] == new_indices[k]:
                        index_positions[k] = z

            counter = found - 1

            for i in range(0, len(indices)):

                are_used = False

                for z in range(0, found):

                    if i == index_positions[z]:
                        are_used = True

                if not are_used:
                    counter = counter + 1
                    new_indices[counter] = indices[i]

            graph.vertices[y].neighbour_indices = new_indices
            graph.vertices[y].is_unpopular = True
            graph.vertices[y].is_popular = False

            friend_index = -1
            distance = 100000

            for x in range(0, len(indices)):

                if not index_of_unpopular_neighbours[x] == -1:

                    temp_distance = np.sqrt((graph.vertices[y].real_coor_x -
                                             graph.vertices[index_of_unpopular_neighbours[x]].real_coor_x)**2 +
                                            (graph.vertices[y].real_coor_y -
                                             graph.vertices[index_of_unpopular_neighbours[x]].real_coor_y)**2)

                    if temp_distance < distance:
                        distance = temp_distance
                        friend_index = index_of_unpopular_neighbours[x]

            if not friend_index == -1:

                i_1 = -1

                for j in range(0, len(indices)):

                    if new_indices[j] == friend_index:
                        i_1 = j

                print('y: ' + str(y) + ', found: ' + str(found) + ', i_1: ' + str(i_1) + ', friend_index: ' + str(friend_index))

                graph.vertices[y].neighbour_indices[i_1] = graph.vertices[y].neighbour_indices[found]
                graph.vertices[y].neighbour_indices[found] = friend_index

                find_consistent_perturbations_simple(graph, friend_index)
                find_consistent_perturbations_simple(graph, y)

            else:

                distance = 10000
                friend_index = -1

                for x in range(found, len(indices)):

                    if not graph.vertices[graph.vertices[y].neighbour_indices[x]].level == graph.vertices[y].level:

                        temp_distance = np.sqrt((graph.vertices[y].real_coor_x -
                                                 graph.vertices[graph.vertices[y].neighbour_indices[x]].real_coor_x) ** 2 +
                                                (graph.vertices[y].real_coor_y -
                                                 graph.vertices[graph.vertices[y].neighbour_indices[x]].real_coor_y) ** 2)

                        if temp_distance < distance:
                            distance = temp_distance
                            friend_index = graph.vertices[y].neighbour_indices[x]

                if not friend_index == -1:

                    i_1 = -1

                    for j in range(0, len(indices)):

                        if new_indices[j] == friend_index:
                            i_1 = j

                    graph.vertices[y].neighbour_indices[i_1] = graph.vertices[y].neighbour_indices[found]
                    graph.vertices[y].neighbour_indices[found] = friend_index

                    find_consistent_perturbations_simple(graph, friend_index)
                    find_consistent_perturbations_simple(graph, y)

            if experimental:
                sort_neighbours_by_level(graph, y)


def connection_shift_on_level(graph, i, experimental=False):

    n = graph.vertices[i].n()
    indices = graph.vertices[i].neighbour_indices

    bad_index = -1
    good_index = -1

    for x in range(0, n):

        if graph.vertices[indices[x]].level == graph.vertices[i].level:

            bad_index = x

    if experimental:

        high = n + 1

    else:

        high = 8

    for x in range(n, high):

        if not graph.vertices[indices[n + high - 1 - x]].level == graph.vertices[i].level:

            good_index = n + high - 1 - x

    if not bad_index == -1 and not good_index == -1:
        perturbator(graph, i, bad_index, good_index)


def reset_popularity_flags(graph):

    for x in range(0, graph.num_vertices):
        graph.vertices[x].is_popular = False
        graph.vertices[x].is_unpopular = False

    graph.num_unpopular = 0
    graph.num_popular = 0
    graph.num_inconsistencies = 0


def find_shape(graph, i, j, clockwise=True):

    closed = False
    start_index = i
    shape_indices = np.ndarray([2], dtype=int)
    shape_indices[0] = i
    shape_indices[1] = j

    while not closed:

        i = shape_indices[shape_indices.shape[0] - 2]
        j = shape_indices[shape_indices.shape[0] - 1]

        if j == start_index or shape_indices.shape[0] > 7:

            closed = True

        else:

            next_index = -1

            if not graph.test_reciprocality(i, j):

                if clockwise:
                    sorted_indices, alpha = clockwise_neighbour_sort(graph, j, j=i)
                else:
                    sorted_indices, alpha = anticlockwise_neighbour_sort(graph, j, j=i)

                next_index = graph.vertices[j].n()

            else:

                if clockwise:
                    sorted_indices, alpha = clockwise_neighbour_sort(graph, j)
                else:
                    sorted_indices, alpha = anticlockwise_neighbour_sort(graph, j)

                for x in range(0, graph.vertices[j].n()):

                    if sorted_indices[x] == i:

                        if x == 0:
                            next_index = graph.vertices[j].n() - 1
                        else:
                            next_index = x - 1

            next_index = sorted_indices[next_index]

            shape_indices = np.append(shape_indices, next_index)

    return shape_indices, shape_indices.shape[0] - 1


def clockwise_neighbour_sort(graph, i, j=-1):

    n = graph.vertices[i].n()

    if not j == -1:
        n = graph.vertices[i].n() + 1

    a = np.ndarray([n], dtype=np.int)
    b = np.ndarray([n], dtype=np.int)
    alpha = np.ndarray([n - 1], dtype=np.float64)
    indices = np.ndarray([n - 1], dtype=int)
    sorted_indices = np.ndarray([n], dtype=int)

    if not j == -1:

        sorted_indices[0] = j

        a[0] = graph.vertices[j].real_coor_x - graph.vertices[i].real_coor_x
        b[0] = graph.vertices[j].real_coor_y - graph.vertices[i].real_coor_y

        for x in range(1, n):
            a[x] = graph.vertices[graph.vertices[i].neighbour_indices[x - 1]].real_coor_x - graph.vertices[i].real_coor_x
            b[x] = graph.vertices[graph.vertices[i].neighbour_indices[x - 1]].real_coor_y - graph.vertices[i].real_coor_y

        for x in range(0, n - 1):
            indices[x] = graph.vertices[i].neighbour_indices[x]

        for x in range(1, n):

            alpha[x - 1] = utils.find_angle(a[0], a[x], b[0], b[x])

            if utils.vector_cross_product_magnitude(a[0], a[x], b[0], b[x]) < 0:
                alpha[x - 1] = 2 * np.pi - alpha[x - 1]

        alpha, indices = utils.dual_sort(alpha, indices)

        for x in range(0, n - 1):
            sorted_indices[x + 1] = indices[x]

        alpha = np.append(alpha, 2 * np.pi)

    else:

        sorted_indices[0] = graph.vertices[i].neighbour_indices[0]

        for x in range(1, n):
            indices[x - 1] = graph.vertices[i].neighbour_indices[x]

        for x in range(0, n):
            a[x] = graph.vertices[graph.vertices[i].neighbour_indices[x]].real_coor_x - graph.vertices[i].real_coor_x
            b[x] = graph.vertices[graph.vertices[i].neighbour_indices[x]].real_coor_y - graph.vertices[i].real_coor_y

        for x in range(1, n):

            indices[x - 1] = graph.vertices[i].neighbour_indices[x]

            alpha[x - 1] = utils.find_angle(a[0], a[x], b[0], b[x])

            if utils.vector_cross_product_magnitude(a[0], a[x], b[0], b[x]) < 0:

                alpha[x - 1] = 2 * np.pi - alpha[x - 1]

        alpha, indices = utils.dual_sort(alpha, indices)

        for x in range(0, n - 1):

            sorted_indices[x + 1] = indices[x]

        alpha = np.append(alpha, 2 * np.pi)

    return sorted_indices, alpha


def anticlockwise_neighbour_sort(graph, i, j=-1):

    n = graph.vertices[i].n()

    if not j == -1:
        n = graph.vertices[i].n() + 1

    a = np.ndarray([n], dtype=np.int)
    b = np.ndarray([n], dtype=np.int)
    alpha = np.ndarray([n - 1], dtype=np.float64)
    indices = np.ndarray([n - 1], dtype=int)
    sorted_indices = np.ndarray([n], dtype=int)

    if not j == -1:

        sorted_indices[0] = j

        a[0] = graph.vertices[j].real_coor_x - graph.vertices[i].real_coor_x
        b[0] = graph.vertices[j].real_coor_y - graph.vertices[i].real_coor_y

        for x in range(1, n):
            a[x] = graph.vertices[graph.vertices[i].neighbour_indices[x - 1]].real_coor_x - graph.vertices[i].real_coor_x
            b[x] = graph.vertices[graph.vertices[i].neighbour_indices[x - 1]].real_coor_y - graph.vertices[i].real_coor_y

        for x in range(0, n - 1):
            indices[x] = graph.vertices[i].neighbour_indices[x]

        for x in range(1, n):

            alpha[x - 1] = utils.find_angle(a[0], a[x], b[0], b[x])

            if utils.vector_cross_product_magnitude(a[0], a[x], b[0], b[x]) > 0:
                alpha[x - 1] = 2 * np.pi - alpha[x - 1]

        alpha, indices = utils.dual_sort(alpha, indices)

        for x in range(0, n - 1):
            sorted_indices[x + 1] = indices[x]

        alpha = np.append(alpha, 2 * np.pi)

    else:

        sorted_indices[0] = graph.vertices[i].neighbour_indices[0]

        for x in range(1, n):
            indices[x - 1] = graph.vertices[i].neighbour_indices[x]

        for x in range(0, n):
            a[x] = graph.vertices[graph.vertices[i].neighbour_indices[x]].real_coor_x - graph.vertices[i].real_coor_x
            b[x] = graph.vertices[graph.vertices[i].neighbour_indices[x]].real_coor_y - graph.vertices[i].real_coor_y

        for x in range(1, n):

            indices[x - 1] = graph.vertices[i].neighbour_indices[x]

            alpha[x - 1] = utils.find_angle(a[0], a[x], b[0], b[x])

            if utils.vector_cross_product_magnitude(a[0], a[x], b[0], b[x]) > 0:

                alpha[x - 1] = 2 * np.pi - alpha[x - 1]

        alpha, indices = utils.dual_sort(alpha, indices)

        for x in range(0, n - 1):

            sorted_indices[x + 1] = indices[x]

        alpha = np.append(alpha, 2 * np.pi)

    return sorted_indices, alpha


