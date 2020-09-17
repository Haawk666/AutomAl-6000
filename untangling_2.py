# By Haakon Tvedt @ NTNU
# Contributors:

import logging
import graph_2
# Instantiate logger:
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def determine_sub_graph_class(sub_graphs):
    for sub_graph in sub_graphs:
        if sub_graph.meshes[0].order == 3 and sub_graph.meshes[1].order == 3:
            sub_graph.class_ = 1
        elif sub_graph.meshes[0].order == 4 and sub_graph.meshes[1].order == 3:
            sub_graph.class_ = 2
        elif sub_graph.meshes[0].order == 3 and sub_graph.meshes[1].order == 4:
            sub_graph.class_ = 3
        elif sub_graph.meshes[0].order == 3 and sub_graph.meshes[1].order == 5:
            sub_graph.class_ = 4
        elif sub_graph.meshes[0].order == 5 and sub_graph.meshes[1].order == 3:
            sub_graph.class_ = 5
        elif sub_graph.meshes[0].order == 4 and sub_graph.meshes[1].order == 4:
            sub_graph.class_ = 6


def determine_sub_graph_configuration(atomic_graph, sub_graphs):
    for sub_graph in sub_graphs:
        if sub_graph.class_ == 1:
            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            a = sub_graph.meshes[0].vertex_indices[2]
            b = sub_graph.meshes[1].vertex_indices[2]
            s = [(j, a), (a, i), (i, b), (b, j)]
            s_types = []

            for f in range(0, 4):
                if atomic_graph.vertices[s[f][0]].partner_query(s[f][1]):
                    if atomic_graph.vertices[s[f][1]].partner_query(s[f][0]):
                        s_types.append(0)
                    else:
                        s_types.append(1)
                else:
                    s_types.append(2)

            if s_types == [0, 0, 0, 0]:
                sub_graph.configuration = 'A_1'
            elif s_types == [0, 2, 0, 2]:
                sub_graph.configuration = 'B_1'
            elif s_types == [1, 0, 1, 0]:
                sub_graph.configuration = 'B_2'
            elif s_types == [0, 1, 0, 0]:
                sub_graph.configuration = 'C_1'
            elif s_types == [0, 0, 2, 0]:
                sub_graph.configuration = 'C_2'
            elif s_types == [0, 1, 1, 0]:
                sub_graph.configuration = 'D_1'
            elif s_types == [0, 2, 2, 0]:
                sub_graph.configuration = 'D_2'
            elif s_types == [0, 2, 0, 1]:
                sub_graph.configuration = 'E_1'
            elif s_types == [2, 0, 1, 0]:
                sub_graph.configuration = 'E_2'
            elif s_types == [0, 2, 0, 0]:
                sub_graph.configuration = 'F_1'
            elif s_types == [0, 0, 1, 0]:
                sub_graph.configuration = 'F_2'
            elif s_types == [0, 0, 0, 2]:
                sub_graph.configuration = 'G_1'
            elif s_types == [1, 0, 0, 0]:
                sub_graph.configuration = 'G_2'
            elif s_types == [2, 1, 0, 0]:
                sub_graph.configuration = 'H_1'
            elif s_types == [0, 0, 2, 1]:
                sub_graph.configuration = 'H_2'
            elif s_types == [0, 1, 2, 0]:
                sub_graph.configuration = 'I_1'

        elif sub_graph.class_ == 2:
            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            a = sub_graph.meshes[0].vertex_indices[2]
            b = sub_graph.meshes[0].vertex_indices[3]
            c = sub_graph.meshes[1].vertex_indices[2]
            s = [(j, a), (a, b), (b, i), (i, c), (c, j)]
            s_types = []

            for f in range(0, 5):
                if atomic_graph.vertices[s[f][0]].partner_query(s[f][1]):
                    if atomic_graph.vertices[s[f][1]].partner_query(s[f][0]):
                        s_types.append(0)
                    else:
                        s_types.append(1)
                else:
                    s_types.append(2)

            if s_types == [0, 0, 0, 0, 0]:
                sub_graph.configuration = 'A_1'
            elif s_types == [0, 1, 0, 0, 0]:
                sub_graph.configuration = 'B_1'
            elif s_types == [0, 1, 0, 1, 0]:
                sub_graph.configuration = 'C_1'
            elif s_types == [0, 0, 0, 0, 2]:
                sub_graph.configuration = 'D_1'
            elif s_types == [0, 0, 1, 0, 0]:
                sub_graph.configuration = 'E_1'

        elif sub_graph.class_ == 3:
            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            a = sub_graph.meshes[0].vertex_indices[2]
            b = sub_graph.meshes[1].vertex_indices[2]
            c = sub_graph.meshes[1].vertex_indices[3]
            s = [(j, a), (a, i), (i, b), (b, c), (c, j)]
            s_types = []

            for f in range(0, 5):
                if atomic_graph.vertices[s[f][0]].partner_query(s[f][1]):
                    if atomic_graph.vertices[s[f][1]].partner_query(s[f][0]):
                        s_types.append(0)
                    else:
                        s_types.append(1)
                else:
                    s_types.append(2)

            if s_types == [0, 0, 0, 0, 0]:
                sub_graph.configuration = 'A_1'
            elif s_types == [0, 0, 0, 2, 0]:
                sub_graph.configuration = 'B_1'
            elif s_types == [0, 2, 0, 2, 0]:
                sub_graph.configuration = 'C_1'
            elif s_types == [1, 0, 0, 0, 0]:
                sub_graph.configuration = 'D_1'
            elif s_types == [0, 0, 2, 0, 0]:
                sub_graph.configuration = 'E_1'

        elif sub_graph.class_ == 4:
            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            a = sub_graph.meshes[0].vertex_indices[2]
            b = sub_graph.meshes[1].vertex_indices[2]
            c = sub_graph.meshes[1].vertex_indices[3]
            d = sub_graph.meshes[1].vertex_indices[4]
            s = [(j, a), (a, i), (i, b), (b, c), (c, d), (d, j)]
            s_types = []

            for f in range(0, 6):
                if atomic_graph.vertices[s[f][0]].partner_query(s[f][1]):
                    if atomic_graph.vertices[s[f][1]].partner_query(s[f][0]):
                        s_types.append(0)
                    else:
                        s_types.append(1)
                else:
                    s_types.append(2)

            if s_types == [0, 0, 0, 0, 0, 0]:
                sub_graph.configuration = 'A_1'
            elif s_types == [0, 0, 0, 0, 2, 0]:
                sub_graph.configuration = 'B_1'

        elif sub_graph.class_ == 5:
            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            a = sub_graph.meshes[0].vertex_indices[2]
            b = sub_graph.meshes[0].vertex_indices[3]
            c = sub_graph.meshes[0].vertex_indices[4]
            d = sub_graph.meshes[1].vertex_indices[2]
            s = [(j, a), (a, b), (b, c), (c, i), (i, d), (d, j)]
            s_types = []

            for f in range(0, 6):
                if atomic_graph.vertices[s[f][0]].partner_query(s[f][1]):
                    if atomic_graph.vertices[s[f][1]].partner_query(s[f][0]):
                        s_types.append(0)
                    else:
                        s_types.append(1)
                else:
                    s_types.append(2)

            if s_types == [0, 0, 0, 0, 0, 0]:
                sub_graph.configuration = 'A_1'
            elif s_types == [0, 1, 0, 0, 0, 0]:
                sub_graph.configuration = 'B_1'

        elif sub_graph.class_ == 6:
            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            a = sub_graph.meshes[0].vertex_indices[2]
            b = sub_graph.meshes[0].vertex_indices[3]
            c = sub_graph.meshes[1].vertex_indices[2]
            d = sub_graph.meshes[1].vertex_indices[3]
            s = [(j, a), (a, b), (b, i), (i, c), (c, d), (d, j)]
            s_types = []

            for f in range(0, 6):
                if atomic_graph.vertices[s[f][0]].partner_query(s[f][1]):
                    if atomic_graph.vertices[s[f][1]].partner_query(s[f][0]):
                        s_types.append(0)
                    else:
                        s_types.append(1)
                else:
                    s_types.append(2)

            if s_types == [0, 0, 0, 0, 0, 0]:
                sub_graph.configuration = 'A_1'
            elif s_types == [0, 0, 0, 0, 0, 1]:
                sub_graph.configuration = 'B_1'
            elif s_types == [2, 0, 0, 0, 0, 0]:
                sub_graph.configuration = 'B_2'
            elif s_types == [2, 0, 0, 0, 0, 1]:
                sub_graph.configuration = 'C_1'
            elif s_types == [0, 0, 0, 2, 0, 0]:
                sub_graph.configuration = 'D_1'
            elif s_types == [0, 0, 1, 0, 0, 0]:
                sub_graph.configuration = 'D_2'
            elif s_types == [0, 0, 1, 2, 0, 0]:
                sub_graph.configuration = 'E_1'


def weak_resolve(graph_obj, sub_graphs, ui_obj=None):
    for sub_graph in sub_graphs:
        separation_cut_off = 350

        if sub_graph.class_ == 1:

            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            a = sub_graph.meshes[0].vertex_indices[2]
            b = sub_graph.meshes[1].vertex_indices[2]

            if sub_graph.configuration == 'A_1':
                k = graph_obj.weak_remove_edge(i, j, aggressive=True)
                if not k == -1:
                    if graph_obj.get_projected_separation(i, k) < separation_cut_off:
                        if graph_obj.permute_j_k(i, j, k):
                            if ui_obj is not None:
                                ui_obj.gs_atomic_graph.perturb_edge(i, j, k, permute_data=False, center_view=True)
            elif sub_graph.configuration == 'B_1':
                pass

            elif sub_graph.configuration == 'B_2':
                pass

            elif sub_graph.configuration == 'C_1':
                if graph_obj.vertices[i].partner_query(j):
                    if graph_obj.get_projected_separation(i, a) < separation_cut_off:
                        if graph_obj.permute_j_k(i, j, a):
                            if ui_obj is not None:
                                ui_obj.gs_atomic_graph.perturb_edge(i, j, a, permute_data=False, center_view=True)
            elif sub_graph.configuration == 'C_2':
                if graph_obj.vertices[i].partner_query(j):
                    if graph_obj.get_projected_separation(i, b) < separation_cut_off:
                        if graph_obj.permute_j_k(i, j, b):
                            if ui_obj is not None:
                                ui_obj.gs_atomic_graph.perturb_edge(i, j, b, permute_data=False, center_view=True)
            elif sub_graph.configuration == 'D_1':
                if graph_obj.vertices[i].partner_query(j):
                    if graph_obj.get_projected_separation(i, a) < separation_cut_off:
                        if graph_obj.permute_j_k(i, j, a):
                            if ui_obj is not None:
                                ui_obj.gs_atomic_graph.perturb_edge(i, j, a, permute_data=False, center_view=True)
            elif sub_graph.configuration == 'D_2':
                if graph_obj.vertices[i].partner_query(j):
                    if graph_obj.get_projected_separation(i, b) < separation_cut_off:
                        if graph_obj.permute_j_k(i, j, b):
                            if ui_obj is not None:
                                ui_obj.gs_atomic_graph.perturb_edge(i, j, b, permute_data=False, center_view=True)
            elif sub_graph.configuration == 'E_1' or sub_graph.configuration == 'E_2':
                pass
            elif sub_graph.configuration == 'F_1' or sub_graph.configuration == 'F_2':
                k = graph_obj.weak_remove_edge(i, j, aggressive=False)
                if not k == -1:
                    if graph_obj.get_projected_separation(i, k) < separation_cut_off:
                        if graph_obj.permute_j_k(i, j, k):
                            if ui_obj is not None:
                                ui_obj.gs_atomic_graph.perturb_edge(i, j, k, permute_data=False, center_view=True)

            elif sub_graph.configuration == 'G_1':
                pass

            elif sub_graph.configuration == 'G_2':
                pass

            elif sub_graph.configuration == 'H_1':
                if graph_obj.vertices[i].partner_query(j):
                    if graph_obj.get_projected_separation(i, a) < separation_cut_off:
                        if graph_obj.permute_j_k(i, j, a):
                            if ui_obj is not None:
                                ui_obj.gs_atomic_graph.perturb_edge(i, j, a, permute_data=False, center_view=True)
            elif sub_graph.configuration == 'H_2':
                if graph_obj.vertices[i].partner_query(j):
                    if graph_obj.get_projected_separation(i, b) < separation_cut_off:
                        if graph_obj.permute_j_k(i, j, b):
                            if ui_obj is not None:
                                ui_obj.gs_atomic_graph.perturb_edge(i, j, b, permute_data=False, center_view=True)
            elif sub_graph.configuration == 'I_1':
                if graph_obj.get_projected_separation(i, b) < separation_cut_off:
                    if graph_obj.permute_j_k(i, j, b):
                        if ui_obj is not None:
                            ui_obj.gs_atomic_graph.perturb_edge(i, j, b, permute_data=False, center_view=True)

        elif sub_graph.class_ == 2:

            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            a = sub_graph.meshes[0].vertex_indices[2]
            b = sub_graph.meshes[0].vertex_indices[3]
            c = sub_graph.meshes[1].vertex_indices[2]

            if sub_graph.configuration == 'A_1':
                pass

            elif sub_graph.configuration == 'B_1' or sub_graph.configuration == 'C_1':
                if graph_obj.vertices[i].partner_query(j) and graph_obj.vertices[a].partner_query(b):
                    if graph_obj.get_projected_separation(i, a) < separation_cut_off:
                        if graph_obj.permute_j_k(i, j, a):
                            if ui_obj is not None:
                                ui_obj.gs_atomic_graph.perturb_edge(i, j, a, permute_data=False, center_view=True)

            elif sub_graph.configuration == 'D_1':
                pass

            elif sub_graph.configuration == 'E_1':
                pass

        elif sub_graph.class_ == 3:

            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            a = sub_graph.meshes[0].vertex_indices[2]
            b = sub_graph.meshes[1].vertex_indices[2]
            c = sub_graph.meshes[1].vertex_indices[3]

            if sub_graph.configuration == 'A_1':
                pass

            elif sub_graph.configuration == 'B_1' or sub_graph.configuration == 'C_1':
                if graph_obj.vertices[i].partner_query(j) and graph_obj.vertices[c].partner_query(b):
                    if graph_obj.get_projected_separation(i, c) < separation_cut_off:
                        if graph_obj.permute_j_k(i, j, c):
                            if ui_obj is not None:
                                ui_obj.gs_atomic_graph.perturb_edge(i, j, c, permute_data=False, center_view=True)

            elif sub_graph.configuration == 'D_1':
                pass

            elif sub_graph.configuration == 'E_1':
                pass

        elif sub_graph.class_ == 4:

            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            d = sub_graph.meshes[1].vertex_indices[4]

            if sub_graph.configuration == 'A_1':
                if graph_obj.get_projected_separation(i, d) < separation_cut_off:
                    if graph_obj.permute_j_k(i, j, d):
                        if ui_obj is not None:
                            ui_obj.gs_atomic_graph.perturb_edge(i, j, d, permute_data=False, center_view=True)

            elif sub_graph.configuration == 'B_1':
                if graph_obj.get_projected_separation(i, d) < separation_cut_off:
                    if graph_obj.permute_j_k(i, j, d):
                        if ui_obj is not None:
                            ui_obj.gs_atomic_graph.perturb_edge(i, j, d, permute_data=False, center_view=True)

        elif sub_graph.class_ == 5:

            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            a = sub_graph.meshes[0].vertex_indices[2]

            if sub_graph.configuration == 'A_1':
                if graph_obj.get_projected_separation(i, a) < separation_cut_off:
                    if graph_obj.permute_j_k(i, j, a):
                        if ui_obj is not None:
                            ui_obj.gs_atomic_graph.perturb_edge(i, j, a, permute_data=False, center_view=True)

            elif sub_graph.configuration == 'B_1':
                if graph_obj.get_projected_separation(i, a) < separation_cut_off:
                    if graph_obj.permute_j_k(i, j, a):
                        if ui_obj is not None:
                            ui_obj.gs_atomic_graph.perturb_edge(i, j, a, permute_data=False, center_view=True)

        elif sub_graph.class_ == 6:

            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            a = sub_graph.meshes[0].vertex_indices[2]
            b = sub_graph.meshes[0].vertex_indices[3]
            c = sub_graph.meshes[1].vertex_indices[2]
            d = sub_graph.meshes[1].vertex_indices[3]

            if sub_graph.configuration == 'A_1':
                k = graph_obj.weak_preserve_edge(i, j)
                if not k == -1:
                    if graph_obj.get_projected_separation(i, j) < separation_cut_off:
                        if graph_obj.permute_j_k(j, k, i):
                            if ui_obj is not None:
                                ui_obj.gs_atomic_graph.perturb_edge(j, k, i, permute_data=False, center_view=True)

            elif sub_graph.configuration == 'B_1' or sub_graph.configuration == 'B_2':
                k = graph_obj.weak_preserve_edge(i, j)
                if not k == -1:
                    if graph_obj.get_projected_separation(i, j) < separation_cut_off:
                        if graph_obj.permute_j_k(j, k, i):
                            if ui_obj is not None:
                                ui_obj.gs_atomic_graph.perturb_edge(j, k, i, permute_data=False, center_view=True)

            elif sub_graph.configuration == 'C_1':
                k = graph_obj.weak_preserve_edge(i, j)
                if not k == -1:
                    if graph_obj.get_projected_separation(i, j) < separation_cut_off:
                        if graph_obj.permute_j_k(j, k, i):
                            if ui_obj is not None:
                                ui_obj.gs_atomic_graph.perturb_edge(j, k, i, permute_data=False, center_view=True)


def strong_resolve(graph_obj, sub_graphs, ui_obj=None):
    for sub_graph in sub_graphs:
        if sub_graph.class_ == 1:
            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            a = sub_graph.meshes[0].vertex_indices[2]
            b = sub_graph.meshes[1].vertex_indices[2]

            if sub_graph.configuration == 'A_1':
                graph_obj.strong_remove_arc(i, j)
            elif sub_graph.configuration == 'B_1':
                pass
            elif sub_graph.configuration == 'B_2':
                pass
            elif sub_graph.configuration == 'C_1':
                pass
            elif sub_graph.configuration == 'C_2':
                pass
            elif sub_graph.configuration == 'D_1':
                pass
            elif sub_graph.configuration == 'D_2':
                pass
            elif sub_graph.configuration == 'E_1' or sub_graph.configuration == 'E_2':
                pass
            elif sub_graph.configuration == 'F_1' or sub_graph.configuration == 'F_2':
                graph_obj.strong_remove_arc(i, j)
            elif sub_graph.configuration == 'G_1':
                pass
            elif sub_graph.configuration == 'G_2':
                pass
            elif sub_graph.configuration == 'H_1':
                pass
            elif sub_graph.configuration == 'H_2':
                pass

        elif sub_graph.class_ == 2:
            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            a = sub_graph.meshes[0].vertex_indices[2]
            b = sub_graph.meshes[0].vertex_indices[3]

            if sub_graph.configuration == 'A_1':
                pass
            elif sub_graph.configuration == 'B_1':
                if graph_obj.vertices[i].partner_query(j) and graph_obj.vertices[a].partner_query(b):
                    graph_obj.permute_j_k(i, j, a)
                    graph_obj.permute_j_k(a, b, i)

        elif sub_graph.class_ == 3:
            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            b = sub_graph.meshes[1].vertex_indices[2]
            c = sub_graph.meshes[1].vertex_indices[3]

            if sub_graph.configuration == 'A_1':
                pass
            elif sub_graph.configuration == 'B_1':
                if graph_obj.vertices[i].partner_query(j) and graph_obj.vertices[c].partner_query(b):
                    graph_obj.permute_j_k(i, j, c)
                    graph_obj.permute_j_k(c, b, i)

        elif sub_graph.class_ == 4:
            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            d = sub_graph.meshes[1].vertex_indices[4]

            if sub_graph.configuration == 'A_1':
                if graph_obj.vertices[i].partner_query(j):
                    graph_obj.permute_j_k(i, j, d)

        elif sub_graph.class_ == 5:
            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]
            a = sub_graph.meshes[0].vertex_indices[2]

            if sub_graph.configuration == 'A_1':
                if graph_obj.vertices[i].partner_query(j):
                    graph_obj.permute_j_k(i, j, a)

        elif sub_graph.class_ == 6:
            i = sub_graph.meshes[0].vertex_indices[0]
            j = sub_graph.meshes[0].vertex_indices[1]

            if sub_graph.configuration == 'A_1' or sub_graph.configuration == 'B_1' or \
                    sub_graph.configuration == 'B_2' or sub_graph.configuration == 'C_1' or \
                    sub_graph.configuration == 'D_1' or \
                    sub_graph.configuration == 'D_2' or sub_graph.configuration == 'E_1':
                graph_obj.strong_preserve_arc(i, j)


