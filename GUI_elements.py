# By Haakon Tvedt @ NTNU
# Contributors:
"""Module container for high-level custom GUI-elements"""

# Program imports:
import GUI_custom_components
import GUI_settings
import GUI_tooltips
import GUI
import utils
import data_module
# External imports:
import numpy as np
import copy
from PyQt5 import QtWidgets, QtGui, QtCore
import pathlib
import os
import logging
# Instantiate logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


# ----------
# Custom QGraphicsScene and QGraphicsView classes:
# ----------


class RawImage(QtWidgets.QGraphicsScene):

    def __init__(self, *args, ui_obj=None, background=None):
        """Initialize a custom QtWidgets.QGraphicsScene object for **raw image**."""

        super().__init__(*args)

        self.ui_obj = ui_obj
        self.background_image = background
        if self.background_image is not None:
            self.addPixmap(self.background_image)
        self.scale_bar = None
        if self.ui_obj.project_instance is not None:
            self.scale_bar = GUI_custom_components.ScaleBar(length=2, scale=self.ui_obj.project_instance.scale, r=self.ui_obj.project_instance.r, height=self.ui_obj.project_instance.im_height)
            self.addItem(self.scale_bar)
            if self.ui_obj.control_window.widgets['overlay']['chb_raw_image']['widget'].isChecked():
                self.scale_bar.show()
            else:
                self.scale_bar.hide()


class StaticImage(QtWidgets.QGraphicsScene):

    def __init__(self, *args, ui_obj=None, background=None):
        """Initialize a custom QtWidgets.QGraphicsScene object for **raw image**."""

        super().__init__(*args)

        self.ui_obj = ui_obj
        self.background_image = background
        if self.background_image is not None:
            self.addPixmap(self.background_image)


class AtomicPositions(QtWidgets.QGraphicsScene):

    def __init__(self, *args, ui_obj=None, background=None):
        """Initialize a custom QtWidgets.QGraphicsScene object for **atomic positions**."""

        super().__init__(*args)

        self.ui_obj = ui_obj
        self.interactive_position_objects = []
        self.background_image = background
        if self.background_image is not None:
            self.addPixmap(self.background_image)

    def re_draw(self):
        """Redraw contents."""
        self.interactive_position_objects = []
        if self.ui_obj.project_instance is not None:
            for vertex in self.ui_obj.project_instance.graph.vertices:
                self.interactive_position_objects.append(GUI_custom_components.InteractivePosColumn(
                    ui_obj=self.ui_obj,
                    vertex=vertex,
                    r=self.ui_obj.overlay_settings['overlay_radii'] * vertex.r
                ))
                self.addItem(self.interactive_position_objects[-1])
                if not self.ui_obj.control_window.widgets['column_detection']['chb_toggle_positions']['widget'].isChecked():
                    self.interactive_position_objects[-1].hide()
                else:
                    self.interactive_position_objects[-1].show()
                if vertex.void:
                    self.interactive_position_objects[-1].hide()


class HeatMap(QtWidgets.QGraphicsScene):

    def __init__(self, *args, ui_obj=None, background=None, kernel_size=1.0, step_size=10, attribute='normalized_peak_gamma', measure_type='variance', kernel_type='square', legend=True, title='None'):
        super().__init__(*args)

        self.ui_obj = ui_obj
        self.background_image = background
        self.kernel_size = kernel_size
        self.step_size = step_size
        self.measure_type = measure_type
        self.attribute = attribute
        self.kernel_type = kernel_type
        self.legend = None
        self.show_legend = legend
        self.margin = 20
        self.width = 40
        self.heat_map = None
        self.title = title
        if self.background_image is not None:
            self.addPixmap(self.background_image)

    def make(self):
        if self.ui_obj.project_instance is not None:
            self.heat_map = self.ui_obj.project_instance.make_heat_map(
                self.attribute,
                kernel_size=self.kernel_size,
                kernel_step_size=self.step_size,
                measure_type=self.measure_type,
                kernel_type=self.kernel_type,
                title=self.title
            )
            min_ = self.heat_map.min()
            max_ = self.heat_map.max()
            heat_map = utils.normalize_array(self.heat_map, 1)
            utils.im_out_static(heat_map.astype(np.float64), 'Images\Outputs\Buffers\\heat_map.png')
            self.addPixmap(QtGui.QPixmap('Images\Outputs\Buffers\\heat_map.png'))

            self.legend = GUI_custom_components.HeatLegend(
                ui_obj=self.ui_obj,
                min=min_,
                max=max_
            )
            self.addItem(self.legend)
            if self.show_legend:
                self.legend.show()
            else:
                self.legend.hide()


class OverlayComposition(QtWidgets.QGraphicsScene):

    def __init__(self, *args, ui_obj=None, background=None):
        """Initialize a custom QtWidgets.QGraphicsScene object for **overlay composition**."""

        super().__init__(*args)

        self.ui_obj = ui_obj
        self.interactive_overlay_objects = []
        self.background_image = background
        self.pixmap = None

        if self.ui_obj.project_instance is not None:
            self.species_dict = copy.deepcopy(self.ui_obj.project_instance.species_dict)
            self.scale = self.ui_obj.project_instance.scale
        else:
            self.species_dict = copy.deepcopy(GUI.core.Project.default_species_dict)
            self.scale = 1

        if self.ui_obj.overlay_settings['display_mode'] == 'atomic_species':
            for value in self.species_dict['atomic_species'].values():
                value['pixel_radii'] = (value['atomic_radii'] - 15) / self.scale
                value['pen'] = QtGui.QPen(QtGui.QColor(
                    value['color'][0],
                    value['color'][1],
                    value['color'][2]
                ))
                value['pen'].setWidth(30 / self.scale)
                value['zeta_0_brush'] = QtGui.QBrush(QtGui.QColor(
                    value['color'][0],
                    value['color'][1],
                    value['color'][2]
                ))
                value['zeta_1_brush'] = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        else:
            for value in self.species_dict['advanced_species'].values():
                value['pixel_radii'] = (self.species_dict['atomic_species'][value['atomic_species']]['atomic_radii'] - 15) / self.scale
                value['pen'] = QtGui.QPen(QtGui.QColor(
                    value['color'][0],
                    value['color'][1],
                    value['color'][2]
                ))
                value['pen'].setWidth(30 / self.scale)
                value['zeta_0_brush'] = QtGui.QBrush(QtGui.QColor(
                    value['color'][0],
                    value['color'][1],
                    value['color'][2]
                ))
                value['zeta_1_brush'] = QtGui.QBrush(QtGui.QColor(0, 0, 0))

        if background is not None:
            self.pixmap = self.addPixmap(self.background_image)
            if self.ui_obj.control_window.widgets['overlay']['chb_raw_image']['widget'].isChecked():
                self.pixmap.show()
            else:
                self.pixmap.hide()
        if self.ui_obj.control_window.widgets['overlay']['chb_black_background']['widget'].isChecked():
            self.setBackgroundBrush(GUI_settings.brush_black)
        else:
            if GUI_settings.theme == 'dark':
                self.setBackgroundBrush(GUI_settings.background_brush)
            else:
                self.setBackgroundBrush(GUI_settings.brush_white)
        self.scale_bar = None
        if self.ui_obj.project_instance is not None:
            self.scale_bar = GUI_custom_components.ScaleBar(
                length=self.ui_obj.overlay_settings['scalebar_length'],
                scale=self.ui_obj.project_instance.scale,
                r=self.ui_obj.project_instance.r,
                height=self.ui_obj.project_instance.im_height,
                unit=self.ui_obj.overlay_settings['scalebar_unit']
            )
            self.addItem(self.scale_bar)
            self.legend = GUI_custom_components.Legend(
                ui_obj=self.ui_obj
            )
            self.addItem(self.legend)
            self.scale_bar.setZValue(2)
            self.legend.setZValue(3)
            if self.ui_obj.control_window.widgets['overlay']['chb_scalebar']['widget'].isChecked():
                self.scale_bar.show()
            else:
                self.scale_bar.hide()
            if self.ui_obj.control_window.widgets['overlay']['chb_legend']['widget'].isChecked():
                self.legend.show()
            else:
                self.legend.hide()

    def re_draw(self):
        """Redraw contents."""
        self.interactive_overlay_objects = []
        if self.ui_obj.project_instance is not None:
            for vertex in self.ui_obj.project_instance.graph.vertices:
                if self.ui_obj.overlay_settings['display_mode'] == 'atomic_species':
                    species = vertex.atomic_species
                else:
                    species = vertex.advanced_species
                self.interactive_overlay_objects.append(GUI_custom_components.InteractiveOverlayColumn(
                    ui_obj=self.ui_obj,
                    vertex=vertex,
                    r=self.ui_obj.overlay_settings['overlay_radii'] * self.species_dict[self.ui_obj.overlay_settings['display_mode']][species]['pixel_radii'],
                    pen=self.species_dict[self.ui_obj.overlay_settings['display_mode']][species]['pen'],
                    brush=self.species_dict[self.ui_obj.overlay_settings['display_mode']][species]['zeta_{}_brush'.format(vertex.zeta)]
                ))
                self.addItem(self.interactive_overlay_objects[-1])
                if vertex.show_in_overlay:
                    self.interactive_overlay_objects[-1].show()
                else:
                    self.interactive_overlay_objects[-1].hide()
                if vertex.void:
                    self.interactive_overlay_objects[-1].hide()

    def re_draw_vertex(self, i):
        self.removeItem(self.interactive_overlay_objects[i])
        vertex = self.ui_obj.project_instance.graph.vertices[i]
        if self.ui_obj.overlay_settings['display_mode'] == 'atomic_species':
            species = vertex.atomic_species
        else:
            species = vertex.advanced_species
        self.interactive_overlay_objects[i] = GUI_custom_components.InteractiveOverlayColumn(
            ui_obj=self.ui_obj,
            vertex=vertex,
            r=self.ui_obj.overlay_settings['overlay_radii'] * self.species_dict[self.ui_obj.overlay_settings['display_mode']][species]['pixel_radii'],
            pen=self.species_dict[self.ui_obj.overlay_settings['display_mode']][species]['pen'],
            brush=self.species_dict[self.ui_obj.overlay_settings['display_mode']][species]['zeta_{}_brush'.format(vertex.zeta)]
        )
        self.addItem(self.interactive_overlay_objects[i])
        if vertex.show_in_overlay:
            self.interactive_overlay_objects[-1].show()
        else:
            self.interactive_overlay_objects[-1].hide()
        if vertex.void:
            self.interactive_overlay_objects[-1].hide()


class AtomicGraph(QtWidgets.QGraphicsScene):

    def __init__(self, *args, ui_obj=None, background=None, scale_factor=1, mode='district'):
        """Initialize a custom QtWidgets.QGraphicsScene object for **atomic graphs**."""

        super().__init__(*args)
        self.ui_obj = ui_obj
        self.scale_factor = scale_factor
        self.mode = mode
        self.lbl_progress = QtWidgets.QLabel('')
        self.interactive_vertex_objects = []
        self.arcs = []
        self.mesh_details = []
        self.background_image = background
        if self.background_image is not None:
            self.addPixmap(self.background_image)
        if GUI_settings.theme == 'dark':
            self.setBackgroundBrush(GUI_settings.background_brush)
        self.re_draw()

    def perturb_edge(self, i, j, k, permute_data=True, center_view=False):
        """Finds the edge from i to j, and makes it point from i to k."""
        if permute_data:
            if self.mode == 'district':
                self.ui_obj.project_instance.graph.permute_j_k(i, j, k)
            elif self.mode == 'zeta':
                self.ui_obj.project_instance.graph.permute_zeta_j_k(i, j, k)
            else:
                logger.warning('Unknown mode. Using district!')
                self.ui_obj.project_instance.graph.permute_j_k(i, j, k)
        self.redraw_neighbourhood(i)
        if center_view:
            pass
            # self.ui_obj.btn_snap_trigger(column_index=i)
            # self.ui_obj.column_selected(i)
            # print('Hello')

    def eval_style(self, i, m):
        vertex_a = self.ui_obj.project_instance.graph.vertices[i]
        vertex_b = self.ui_obj.project_instance.graph.vertices[self.arcs[i][m].j]
        consistent = vertex_b.partner_query(vertex_a.i)
        if vertex_a.zeta == vertex_b.zeta:
            dislocation = True
        else:
            dislocation = False
        return consistent, dislocation

    def re_draw(self):
        """Redraw contents."""
        if self.ui_obj.project_instance is not None:
            self.re_draw_edges()
            self.re_draw_vertices()
            self.re_draw_mesh_details()

    def re_draw_vertices(self):
        """Redraws all column elements."""
        self.interactive_vertex_objects = []
        for vertex in self.ui_obj.project_instance.graph.vertices:
            self.interactive_vertex_objects.append(GUI_custom_components.InteractiveGraphColumn(
                ui_obj=self.ui_obj,
                vertex=vertex,
                scale_factor=self.scale_factor,
                r=self.ui_obj.overlay_settings['overlay_radii'] * vertex.r
            ))
            if not vertex.void:
                self.addItem(self.interactive_vertex_objects[-1])

    def re_draw_edges(self):
        """Redraws all edge elements."""
        for inner_edges in self.arcs:
            for edge_item in inner_edges:
                self.removeItem(edge_item)
        self.arcs = []

        for vertex_a in self.ui_obj.project_instance.graph.vertices:
            self.arcs.append([])
            if self.mode == 'district':
                partners = self.ui_obj.project_instance.graph.get_vertex_objects_from_indices(
                    vertex_a.partners
                )
                out_semi_partners = self.ui_obj.project_instance.graph.get_vertex_objects_from_indices(
                    vertex_a.out_semi_partners
                )
            elif self.mode == 'zeta':
                partners = self.ui_obj.project_instance.graph.get_vertex_objects_from_indices(
                    vertex_a.local_zeta_map['partners']
                )
                out_semi_partners = self.ui_obj.project_instance.graph.get_vertex_objects_from_indices(
                    vertex_a.local_zeta_map['out_semi_partners']
                )
            else:
                logger.warning('Unknown mode. Using district!')
                partners = self.ui_obj.project_instance.graph.get_vertex_objects_from_indices(
                    vertex_a.partners
                )
                out_semi_partners = self.ui_obj.project_instance.graph.get_vertex_objects_from_indices(
                    vertex_a.out_semi_partners
                )

            for vertex_b in partners:
                p1 = vertex_a.im_pos()
                p2 = vertex_b.im_pos()
                if vertex_a.zeta == vertex_b.zeta:
                    co_planar = True
                else:
                    co_planar = False
                self.arcs[-1].append(GUI_custom_components.Arrow(
                    i=vertex_a.i,
                    j=vertex_b.i,
                    p1=p1,
                    p2=p2,
                    r=self.ui_obj.project_instance.r,
                    scale_factor=self.scale_factor,
                    dual_arc=True,
                    co_planar=co_planar
                ))
                self.addItem(self.arcs[-1][-1])

            for vertex_b in out_semi_partners:
                p1 = vertex_a.im_pos()
                p2 = vertex_b.im_pos()
                self.arcs[-1].append(GUI_custom_components.Arrow(
                    i=vertex_a.i,
                    j=vertex_b.i,
                    p1=p1,
                    p2=p2,
                    r=self.ui_obj.project_instance.r,
                    scale_factor=self.scale_factor,
                    dual_arc=False,
                    co_planar=False
                ))
                self.addItem(self.arcs[-1][-1])

    def re_draw_mesh_details(self):
        """Redraws all mesh details"""
        if self.mode == 'district':
            for mesh_detail in self.mesh_details:
                self.removeItem(mesh_detail)
            self.mesh_details = []

            for mesh in self.ui_obj.project_instance.graph.meshes:
                detail = GUI_custom_components.MeshDetail(mesh=mesh)
                self.addItem(detail)
                self.mesh_details.append(detail)
                if self.ui_obj.control_window.widgets['atomic_graph']['chb_show_mesh_vertex_order']['widget'].isChecked():
                    detail.show()
                else:
                    detail.hide()

    def redraw_star(self, i):
        for edge_item in self.arcs[i]:
            edge_item.hide()
            self.removeItem(edge_item)
        self.arcs[i] = []
        vertex_a = self.ui_obj.project_instance.graph.vertices[i]
        if self.mode == 'district':
            partners = self.ui_obj.project_instance.graph.get_vertex_objects_from_indices(
                vertex_a.partners
            )
            out_semi_partners = self.ui_obj.project_instance.graph.get_vertex_objects_from_indices(
                vertex_a.out_semi_partners
            )
        elif self.mode == 'zeta':
            partners = self.ui_obj.project_instance.graph.get_vertex_objects_from_indices(
                vertex_a.local_zeta_map['partners']
            )
            out_semi_partners = self.ui_obj.project_instance.graph.get_vertex_objects_from_indices(
                vertex_a.local_zeta_map['out_semi_partners']
            )
        else:
            logger.warning('Unknown mode. Using district!')
            partners = self.ui_obj.project_instance.graph.get_vertex_objects_from_indices(
                vertex_a.partners
            )
            out_semi_partners = self.ui_obj.project_instance.graph.get_vertex_objects_from_indices(
                vertex_a.out_semi_partners
            )
        for vertex_b in partners:
            p1 = vertex_a.im_pos()
            p2 = vertex_b.im_pos()
            if vertex_a.zeta == vertex_b.zeta:
                co_planar = True
            else:
                co_planar = False
            self.arcs[i].append(GUI_custom_components.Arrow(
                i=vertex_a.i,
                j=vertex_b.i,
                p1=p1,
                p2=p2,
                r=self.ui_obj.project_instance.r,
                scale_factor=self.scale_factor,
                dual_arc=True,
                co_planar=co_planar
            ))
            self.addItem(self.arcs[i][-1])
        for vertex_b in out_semi_partners:
            p1 = vertex_a.im_pos()
            p2 = vertex_b.im_pos()
            self.arcs[i].append(GUI_custom_components.Arrow(
                i=vertex_a.i,
                j=vertex_b.i,
                p1=p1,
                p2=p2,
                r=self.ui_obj.project_instance.r,
                scale_factor=self.scale_factor,
                dual_arc=False,
                co_planar=False
            ))
            self.addItem(self.arcs[i][-1])

    def redraw_neighbourhood(self, i):
        for neighbours in self.ui_obj.project_instance.graph.vertices[i].district + [i]:
            self.redraw_star(neighbours)


class AtomicSubGraph(QtWidgets.QGraphicsScene):

    def __init__(self, *args, ui_obj=None, background=None, sub_graph=None, scale_factor=1):
        """Initialize a custom QtWidgets.QGraphicsScene object for **atomic sub-graphs**."""

        super().__init__(*args)

        self.ui_obj = ui_obj
        self.scale_factor = scale_factor
        self.interactive_vertex_objects = []
        self.edges = []
        self.vectors = []
        self.labels = []
        self.background_image = background
        self.report = ''
        if self.background_image is not None:
            self.addPixmap(self.background_image)
        if GUI_settings.theme == 'dark':
            self.setBackgroundBrush(GUI_settings.background_brush)
        self.sub_graph = sub_graph
        if sub_graph is not None:
            self.re_draw()

    def re_draw(self):
        """Redraw contents."""
        if self.ui_obj.project_instance is not None:
            self.re_draw_vertices()
            self.re_draw_edges()
            self.label_angles(angle_vectors=False)
            self.relay_summary()

    def re_draw_vertices(self):
        """Redraws all column elements."""
        for vertex in self.sub_graph.vertices:
            self.interactive_vertex_objects.append(GUI_custom_components.InteractiveGraphColumn(self.ui_obj, vertex.i, vertex.r, self.scale_factor))
            self.interactive_vertex_objects[-1].setFlag(QtWidgets.QGraphicsItem.ItemIsSelectable, False)
            self.interactive_vertex_objects[-1].setFlag(QtWidgets.QGraphicsItem.ItemIsPanel, True)
            self.addItem(self.interactive_vertex_objects[-1])

    def re_draw_edges(self):
        """Redraws all edge elements."""
        for edge in self.sub_graph.arcs:
            consistent = edge.dual_arc
            dislocation = edge.co_planar
            p1 = edge.vertex_a.im_pos()
            p2 = edge.vertex_b.im_pos()
            self.edges.append(GUI_custom_components.Arrow(i=edge.vertex_a.i, j=edge.vertex_b.i, p1=p1, p2=p2, r=self.ui_obj.project_instance.r, scale_factor=self.scale_factor, dual_arc=consistent, co_planar=dislocation))
            self.addItem(self.edges[-1])

    def label_angles(self, angle_vectors=False):
        """Label all sub-graph angles."""
        self.report += 'Sub-graph centered on vertex {}:----------\n'.format(self.sub_graph.vertex_indices[0])
        for m, mesh in enumerate(self.sub_graph.meshes):
            self.report += 'Mesh {}:\n'.format(m)
            self.report += '    Number of corners: {}\n'.format(str(mesh.order))
            self.report += '    Sum of angles: {}\n'.format(str(sum(mesh.angles)))
            self.report += '    Variance of angles: {}\n'.format(utils.variance(mesh.angles))
            self.report += '    Symmetry prob vector from central angle: {}\n'.format(str([0, 0, 0]))
            self.report += '    corners: {}\n'.format(mesh.vertex_indices)
            self.report += '    Angles:\n'
            for i, corner in enumerate(mesh.vertices):
                self.report += '        a{}{} = {}\n'.format(m, i, mesh.angles[i])
                p1 = corner.im_pos()
                p1 = (p1[0], p1[1])
                p2 = (p1[0] + 0.5 * corner.r * mesh.angle_vectors[i][0], p1[1] + 0.5 * corner.r * mesh.angle_vectors[i][1])

                if angle_vectors:
                    self.vectors.append(GUI_custom_components.Arrow(p1, p2, corner.r, self.scale_factor, False, False))
                    self.addItem(self.vectors[-1].arrow[0])

                angle_text = QtWidgets.QGraphicsSimpleTextItem()
                angle_text.setText('a{}{}'.format(m, i))
                angle_text.setFont(GUI_settings.font_tiny)
                rect = angle_text.boundingRect()
                angle_text.setPos(self.scale_factor * p2[0] - 0.5 * rect.width(), self.scale_factor * p2[1] - 0.5 * rect.height())
                self.addItem(angle_text)

    def relay_summary(self):
        logger.info(self.report)


class AntiGraph(QtWidgets.QGraphicsScene):

    def __init__(self, *args, ui_obj=None, background=None, scale_factor=1, graph=None):
        """Initialize a custom QtWidgets.QGraphicsScene object for **atomic graphs**."""

        super().__init__(*args)

        self.ui_obj = ui_obj
        self.scale_factor = scale_factor
        self.graph = graph
        self.interactive_vertex_objects = []
        self.edges = []
        self.background_image = background
        if self.background_image is not None:
            self.addPixmap(self.background_image)
        if GUI_settings.theme == 'dark':
            self.setBackgroundBrush(GUI_settings.background_brush)
        self.re_draw()

    def re_draw(self):
        """Redraw contents."""
        if self.graph is not None:
            self.re_draw_edges()
            self.re_draw_vertices()

    def re_draw_vertices(self):
        """Redraws all column elements."""
        self.interactive_vertex_objects = []
        for vertex in self.graph.vertices:
            self.interactive_vertex_objects.append(GUI_custom_components.InteractiveGraphColumn(
                ui_obj=self.ui_obj,
                vertex=vertex,
                scale_factor=self.scale_factor,
                r=self.ui_obj.overlay_settings['overlay_radii'] * vertex.r
            ))
            self.addItem(self.interactive_vertex_objects[-1])

    def re_draw_edges(self):
        """Redraws all edge elements."""
        self.graph.map_arcs()
        for edge_item in self.edges:
            self.removeItem(edge_item)
        self.edges = []
        for edge in self.graph.arcs:
            p1 = edge.vertex_a.im_pos()
            p1 = (p1[0], p1[1])
            p2 = edge.vertex_b.im_pos()
            p2 = (p2[0], p2[1])
            real_distance = edge.spatial_separation
            hard_sphere_distance = edge.hard_sphere_separation
            max_shift = 70
            if hard_sphere_distance > real_distance:
                if hard_sphere_distance - real_distance > max_shift:
                    difference = max_shift
                else:
                    difference = hard_sphere_distance - real_distance
                multiplier = - difference / max_shift + 1
                color = (multiplier * 255, multiplier * 255, 255, 255)
            else:
                if real_distance - hard_sphere_distance > max_shift:
                    difference = max_shift
                else:
                    difference = real_distance - hard_sphere_distance
                multiplier = - difference / max_shift + 1
                color = (255, multiplier * 255, multiplier * 255, 255)

            self.edges.append(GUI_custom_components.DistanceArrow(color=color, i=edge.vertex_a.i, j=edge.vertex_b.i,
                                                                  p1=p1, p2=p2, r=self.ui_obj.project_instance.r,
                                                                  scale_factor=self.scale_factor))
            self.addItem(self.edges[-1])
            if edge.vertex_a.zeta == 1 and edge.vertex_b.zeta == 1:
                pass
            else:
                self.edges[-1].setZValue(-2)

    def toggle_level_0(self, on):
        for vertex in self.graph.vertices:
            if vertex.zeta == 0:
                if on:
                    self.interactive_vertex_objects[vertex.i].show()
                else:
                    self.interactive_vertex_objects[vertex.i].hide()
        for m, edge in enumerate(self.graph.arcs):
            if edge.vertex_a.zeta == 0 and edge.vertex_b.zeta == 0:
                if on:
                    self.edges[m].show()
                else:
                    self.edges[m].hide()

    def toggle_level_1(self, on):
        for vertex in self.graph.vertices:
            if vertex.zeta == 1:
                if on:
                    self.interactive_vertex_objects[vertex.i].show()
                else:
                    self.interactive_vertex_objects[vertex.i].hide()
        for m, edge in enumerate(self.graph.arcs):
            if edge.vertex_a.zeta == 1 and edge.vertex_b.zeta == 1:
                if on:
                    self.edges[m].show()
                else:
                    self.edges[m].hide()


class ZoomGraphicsView(QtWidgets.QGraphicsView):
    """An adaptation of QtWidgets.QGraphicsView that supports zooming"""

    def __init__(self, parent=None, ui_obj=None, trigger_func=None, tab_index=0, scale=1):
        super(ZoomGraphicsView, self).__init__(parent)
        self.ui_obj = ui_obj
        self.trigger_func = trigger_func
        self.tab_index = tab_index
        self.scaling_factor = scale

    def translate_w(self, amount, relay=True):

        # Set Anchors
        self.setTransformationAnchor(QtWidgets.QGraphicsView.NoAnchor)
        self.setResizeAnchor(QtWidgets.QGraphicsView.NoAnchor)

        self.translate(0, amount)

        if self.ui_obj.lock_views and relay:
            for i, zgv in enumerate(self.ui_obj.gv_list):
                if not i == self.tab_index:
                    zgv.translate_w(amount, relay=False)

    def translate_s(self, amount, relay=True):

        # Set Anchors
        self.setTransformationAnchor(QtWidgets.QGraphicsView.NoAnchor)
        self.setResizeAnchor(QtWidgets.QGraphicsView.NoAnchor)

        self.translate(0, -amount)

        if self.ui_obj.lock_views and relay:
            for i, zgv in enumerate(self.ui_obj.gv_list):
                if not i == self.tab_index:
                    zgv.translate_s(amount, relay=False)

    def translate_d(self, amount, relay=True):

        # Set Anchors
        self.setTransformationAnchor(QtWidgets.QGraphicsView.NoAnchor)
        self.setResizeAnchor(QtWidgets.QGraphicsView.NoAnchor)

        self.translate(-amount, 0)

        if self.ui_obj.lock_views and relay:
            for i, zgv in enumerate(self.ui_obj.gv_list):
                if not i == self.tab_index:
                    zgv.translate_d(amount, relay=False)

    def translate_a(self, amount, relay=True):

        # Set Anchors
        self.setTransformationAnchor(QtWidgets.QGraphicsView.NoAnchor)
        self.setResizeAnchor(QtWidgets.QGraphicsView.NoAnchor)

        self.translate(amount, 0)

        if self.ui_obj.lock_views and relay:
            for i, zgv in enumerate(self.ui_obj.gv_list):
                if not i == self.tab_index:
                    zgv.translate_a(amount, relay=False)

    def wheelEvent(self, event, relay=True):

        modifier = QtWidgets.QApplication.keyboardModifiers()
        if modifier == QtCore.Qt.ShiftModifier:

            # Zoom Factor
            zoom_in_factor = 1.25
            zoom_out_factor = 1 / zoom_in_factor

            # Set Anchors
            self.setTransformationAnchor(QtWidgets.QGraphicsView.NoAnchor)
            self.setResizeAnchor(QtWidgets.QGraphicsView.NoAnchor)

            # Save the scene pos
            oldPos = self.mapToScene(event.pos())

            # Zoom
            if event.angleDelta().y() > 0:
                zoomFactor = zoom_in_factor
            else:
                zoomFactor = zoom_out_factor
            self.scale(zoomFactor, zoomFactor)

            # Get the new position
            newPos = self.mapToScene(event.pos())

            # Move scene to old position
            delta = newPos - oldPos
            self.translate(delta.x(), delta.y())

        else:

            super(ZoomGraphicsView, self).wheelEvent(event)

        if self.ui_obj.lock_views and relay:
            for i, zgv in enumerate(self.ui_obj.gv_list):
                if not i == self.tab_index:
                    zgv.wheelEvent(event, relay=False)

    def keyPressEvent(self, event):
        if self.trigger_func is not None:
            self.trigger_func(event.key())


# ----------
# Custom Main window tools:
# ----------


class TerminalTextEdit(QtWidgets.QPlainTextEdit):

    def __init__(self, *args):
        super().__init__(*args)

        self.setReadOnly(True)
        self.setWordWrapMode(QtGui.QTextOption.NoWrap)
        self.setLineWrapMode(QtWidgets.QPlainTextEdit.NoWrap)


class TerminalHandler(logging.Handler):

    def __init__(self):
        super().__init__()

        self.widget = TerminalTextEdit()

        self.setLevel(logging.INFO)
        self.setFormatter(logging.Formatter('%(name)s:%(levelname)s:%(funcName)s:%(message)s'))

    def set_mode(self, debug=False):

        if debug:
            self.setLevel(logging.DEBUG)
            self.setFormatter(logging.Formatter('%(name)s: %(levelname)s: %(funcName)s: %(message)s'))
        else:
            self.setLevel(logging.INFO)
            self.setFormatter(logging.Formatter('%(name)s: %(message)s'))

    def emit(self, record):
        msg = self.format(record)
        self.widget.appendPlainText(msg)
        QtWidgets.QApplication.processEvents()


class Terminal(QtWidgets.QWidget):

    def __init__(self, *args, obj=None):
        super().__init__(*args)

        self.ui_obj = obj

        self.handler = TerminalHandler()

        self.btn_save_log = GUI_custom_components.SmallButton('Save log', self, trigger_func=self.ui_obj.btn_save_log_trigger)
        self.btn_clear_log = GUI_custom_components.SmallButton('Clear log', self, trigger_func=self.ui_obj.btn_clear_log_trigger)
        self.btn_eval = GUI_custom_components.SmallButton('Eval', self, trigger_func=self.eval_input)
        self.btn_eval.setDisabled(True)

        self.chb_input = QtWidgets.QCheckBox('Enable input (advanced)')
        self.chb_input.setChecked(False)
        self.chb_input.toggled.connect(self.set_input)

        self.terminal_btns_layout = QtWidgets.QHBoxLayout()
        self.terminal_btns_layout.addWidget(self.btn_save_log)
        self.terminal_btns_layout.addWidget(self.btn_clear_log)
        self.terminal_btns_layout.addWidget(self.chb_input)
        self.terminal_btns_layout.addStretch()

        self.eval_text = QtWidgets.QLineEdit()
        self.eval_text.setDisabled(True)
        self.input_layout = QtWidgets.QHBoxLayout()
        self.input_layout.addWidget(self.eval_text)
        self.input_layout.addWidget(self.btn_eval)

        self.terminal_display_layout = QtWidgets.QVBoxLayout()
        self.terminal_display_layout.addLayout(self.terminal_btns_layout)
        self.terminal_display_layout.addWidget(self.handler.widget)
        self.terminal_display_layout.addLayout(self.input_layout)

        self.setLayout(self.terminal_display_layout)

        # Set tooltips:
        self.mode_tooltip(self.ui_obj.menu.toggle_tooltips_action.isChecked())

    def mode_tooltip(self, on):
        if on:
            self.btn_save_log.setToolTip('Save the log contents to text-file.')
            self.btn_clear_log.setToolTip('Empty the contents of the log.')
        else:
            self.btn_save_log.setToolTip('')
            self.btn_clear_log.setToolTip('')

    def set_input(self, state):
        self.eval_text.setDisabled(not state)
        self.btn_eval.setDisabled(not state)

    def eval_input(self):
        if self.chb_input.isChecked():
            self.ui_obj.btn_eval_trigger(self.eval_text.text())


class ControlWindow(QtWidgets.QWidget):

    def __init__(self, *args, obj=None):
        super().__init__(*args)

        self.ui_obj = obj

        # Prob. Vector graphic:
        self.height = 130

        self.probGraphicView = QtWidgets.QGraphicsView()

        self.probGraphicView.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.probGraphicView.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.probGraphicView.setMinimumHeight(self.height)
        self.probGraphicView.setMaximumHeight(self.height)
        if GUI_settings.theme == 'dark':
            self.probGraphicView.setBackgroundBrush(GUI_settings.background_brush)

        self.probGraphicLayout = QtWidgets.QHBoxLayout()
        self.probGraphicLayout.addWidget(self.probGraphicView)
        self.probGraphicLayout.addStretch()

        self.draw_histogram()

        # -------------------------------------------------------------------------
        # Define control widgets:
        # -------------------------------------------------------------------------

        self.widgets = {
            'project': {
                'lbl_filename': {
                    'widget': GUI_custom_components.Label(text='Filename'),
                    'value_attribute': 'project_instance.filename_full[1]',
                    'widget_type': 'lbl'},
                'lbl_location': {
                    'widget': GUI_custom_components.Label(text='Location'),
                    'value_attribute': 'project_instance.filename_full[0]',
                    'widget_type': 'lbl'},
                'sbtn_active_model': {
                    'widget': GUI_custom_components.SetButtonLayout(obj=self, trigger_func=self.sbtn_active_model_trigger, text='Associated model', value='default_model'),
                    'value_attribute': 'project_instance.graph.active_model',
                    'widget_type': 'sbtn'},
                'btn_import': {
                    'widget': GUI_custom_components.MediumButton('Import', self, trigger_func=self.btn_import_trigger),
                    'widget_type': 'btn'},
                'btn_open_project': {
                    'widget': GUI_custom_components.MediumButton('Open', self, trigger_func=self.btn_open_project_trigger),
                    'widget_type': 'btn'},
                'btn_save_project': {
                    'widget': GUI_custom_components.MediumButton('Save', self, trigger_func=self.btn_save_project_trigger),
                    'widget_type': 'btn'},
                'btn_close_project': {
                    'widget': GUI_custom_components.MediumButton('Close', self, trigger_func=self.btn_close_project_trigger),
                    'widget_type': 'btn'},
                'btn_species_dict': {
                    'widget': GUI_custom_components.MediumButton('Species dict', self, trigger_func=self.btn_species_dict_trigger),
                    'widget_type': 'btn'},
                'btn_export': {
                    'widget': GUI_custom_components.MediumButton('Export', self, trigger_func=self.btn_export_trigger),
                    'widget_type': 'btn'},
                'btn_show_stats': {
                    'widget': GUI_custom_components.MediumButton('Stats', self, trigger_func=self.btn_show_stats_trigger),
                    'widget_type': 'btn'}
            }, 'image': {
                'lbl_image_width': {
                    'widget': GUI_custom_components.Label(text='Image width (pixels)'),
                    'widget_type': 'lbl',
                    'value_attribute': 'project_instance.im_width'},
                'lbl_image_height': {
                    'widget': GUI_custom_components.Label(text='Image height (pixels)'),
                    'widget_type': 'lbl',
                    'value_attribute': 'project_instance.im_height'},
                'lbl_upsampled': {
                    'widget': GUI_custom_components.Label(text='Automatic bilinear upsampling performed'),
                    'widget_type': 'lbl',
                    'value_attribute': 'project_instance.has_been_up_sampled'},
                'chb_lock_views': {
                    'widget': GUI_custom_components.CheckBox(label='Lock views', trigger_func=self.chb_lock_views_trigger),
                    'widget_type': 'chb',
                    'value_attribute': 'lock_views'},
                'btn_show_source': {
                    'widget': GUI_custom_components.MediumButton('Source', self, trigger_func=self.btn_show_source_trigger),
                    'widget_type': 'btn'},
                'btn_align_views': {
                    'widget': GUI_custom_components.MediumButton('Align views', self, trigger_func=self.btn_align_views_trigger),
                    'widget_type': 'btn'},
            }, 'column_detection': {
                'lbl_atomic_radii': {
                    'widget': GUI_custom_components.Label(text='Approx atomic radii, r (pixels)'),
                    'widget_type': 'lbl',
                    'value_attribute': 'project_instance.r'},
                'lbl_overhead': {
                    'widget': GUI_custom_components.Label(text='Overhead, o (pixels)'),
                    'widget_type': 'lbl',
                    'value_attribute': 'project_instance.overhead'},
                'lbl_pixel_average': {
                    'widget': GUI_custom_components.Label(text='Average pixel value'),
                    'widget_type': 'lbl',
                    'value_attribute': 'project_instance.pixel_average'},
                'lbl_search_matrix_peak': {
                    'widget': GUI_custom_components.Label(text='Search matrix peak'),
                    'widget_type': 'lbl',
                    'value_attribute': 'project_instance.search_mat.max()'},
                'lbl_num_detected_columns': {
                    'widget': GUI_custom_components.Label(text='Number of detected columns'),
                    'widget_type': 'lbl',
                    'value_attribute': 'project_instance.num_columns'},
                'chb_plot_results': {
                    'widget': GUI_custom_components.CheckBox(label='Plot column detection summary'),
                    'widget_type': 'chb'},
                'chb_automatic_threshold': {
                    'widget': GUI_custom_components.CheckBox(label='Use automatic thresholding (experimental)'),
                    'widget_type': 'chb'},
                'chb_toggle_positions': {
                    'widget': GUI_custom_components.CheckBox(label='Show position overlays', trigger_func=self.chb_toggle_positions_trigger, default_state=True),
                    'widget_type': 'chb'},
                'sbtn_threshold': {
                    'widget': GUI_custom_components.SetButtonLayout(obj=self, trigger_func=self.sbtn_threshold_trigger, text='Detection threshold value, T'),
                    'widget_type': 'sbtn',
                    'value_attribute': 'project_instance.threshold'},
                'sbtn_search_size': {
                    'widget': GUI_custom_components.SetButtonLayout(obj=self, trigger_func=self.sbtn_search_size_trigger, text='Search size'),
                    'widget_type': 'sbtn',
                    'value_attribute': 'project_instance.search_size'},
                'sbtn_scale': {
                    'widget': GUI_custom_components.SetButtonLayout(obj=self, trigger_func=self.sbtn_scale_trigger, text='Scale, s (pm / pixel)'),
                    'widget_type': 'sbtn',
                    'value_attribute': 'project_instance.scale'},
                'btn_start_column_detection': {
                    'widget': GUI_custom_components.MediumButton('Start', self, trigger_func=self.btn_start_column_detection_trigger),
                    'widget_type': 'btn'},
                'btn_reset_column_detection': {
                    'widget': GUI_custom_components.MediumButton('Reset', self, trigger_func=self.btn_reset_column_detection_trigger),
                    'widget_type': 'btn'},
                'btn_redraw_search_mat': {
                    'widget': GUI_custom_components.MediumButton('Redraw search mat', self, trigger_func=self.btn_redraw_search_mat_trigger),
                    'widget_type': 'btn'}
            }, 'column_characterization': {
                'chb_show_graphics': {
                    'widget': GUI_custom_components.CheckBox(label='Show graphic updates (slow)'),
                    'widget_type': 'chb'},
                'lbl_alloy': {
                    'widget': GUI_custom_components.Label(text='Alloy'),
                    'widget_type': 'lbl',
                    'value_attribute': 'project_instance.get_alloy_string()'},
                'btn_start_column_characterization': {
                    'widget': GUI_custom_components.MediumButton('Start', self, trigger_func=self.btn_start_column_characterization_trigger),
                    'widget_type': 'btn'},
                'btn_reset_column_characterization': {
                    'widget': GUI_custom_components.MediumButton('Reset', self, trigger_func=self.btn_reset_column_characterization_trigger),
                    'widget_type': 'btn'},
                'btn_invert_zeta': {
                    'widget': GUI_custom_components.MediumButton('Invert lvl', self, trigger_func=self.btn_invert_zeta_trigger),
                    'widget_type': 'btn'}
            }, 'selected_column': {
                'sbtn_find_column': {
                    'widget': GUI_custom_components.SetButtonLayout(obj=self, trigger_func=self.sbtn_find_column_trigger, text='Column index'),
                    'widget_type': 'sbtn',
                    'value_attribute': 'vertex.i'},
                'sbtn_atomic_species': {
                    'widget': GUI_custom_components.SetButtonLayout(obj=self, trigger_func=self.sbtn_atomic_species_trigger, text='Atomic species'),
                    'widget_type': 'sbtn',
                    'value_attribute': 'vertex.atomic_species'},
                'sbtn_advanced_species': {
                    'widget': GUI_custom_components.SetButtonLayout(obj=self, trigger_func=self.sbtn_advanced_species_trigger, text='Advanced species'),
                    'widget_type': 'sbtn',
                    'value_attribute': 'vertex.advanced_species'},
                'sbtn_zeta': {
                    'widget': GUI_custom_components.SetButtonLayout(obj=self, trigger_func=self.sbtn_zeta_trigger, text='Zeta'),
                    'widget_type': 'sbtn',
                    'value_attribute': 'vertex.zeta'},
                'lbl_im_coor_x': {
                    'widget': GUI_custom_components.Label(text='x (pixels)'),
                    'widget_type': 'lbl',
                    'value_attribute': 'vertex.im_coor_x'},
                'lbl_im_coor_y': {
                    'widget': GUI_custom_components.Label(text='y (pixels)'),
                    'widget_type': 'lbl',
                    'value_attribute': 'vertex.im_coor_y'},
                'lbl_spatial_coor_x': {
                    'widget': GUI_custom_components.Label(text='x (pm)'),
                    'widget_type': 'lbl',
                    'value_attribute': 'vertex.spatial_coor_x'},
                'lbl_spatial_coor_y': {
                    'widget': GUI_custom_components.Label(text='y (pm)'),
                    'widget_type': 'lbl',
                    'value_attribute': 'vertex.spatial_coor_y'},
                'lbl_peak_gamma': {
                    'widget': GUI_custom_components.Label(text='Peak gamma'),
                    'widget_type': 'lbl',
                    'value_attribute': 'vertex.peak_gamma'},
                'lbl_avg_gamma': {
                    'widget': GUI_custom_components.Label(text='Avg gamma'),
                    'widget_type': 'lbl',
                    'value_attribute': 'vertex.avg_gamma'},
                'lbl_area_gamma': {
                    'widget': GUI_custom_components.Label(text='Area gamma'),
                    'widget_type': 'lbl',
                    'value_attribute': 'vertex.area_gamma'},
                'lbl_normalized_peak_gamma': {
                    'widget': GUI_custom_components.Label(text='Normalized peak gamma'),
                    'widget_type': 'lbl',
                    'value_attribute': 'vertex.normalized_peak_gamma'},
                'lbl_normalized_avg_gamma': {
                    'widget': GUI_custom_components.Label(text='Normalized avg gamma'),
                    'widget_type': 'lbl',
                    'value_attribute': 'vertex.normalized_avg_gamma'},
                'lbl_normalized_area_gamma': {
                    'widget': GUI_custom_components.Label(text='Normalized area gamma'),
                    'widget_type': 'lbl',
                    'value_attribute': 'vertex.normalized_area_gamma'},
                'lbl_alpha_min': {
                    'widget': GUI_custom_components.Label(text='Alpha min'),
                    'widget_type': 'lbl',
                    'value_attribute': 'vertex.alpha_min'},
                'lbl_alpha_max': {
                    'widget': GUI_custom_components.Label(text='Alpha max'),
                    'widget_type': 'lbl',
                    'value_attribute': 'vertex.alpha_min'},
                'lbl_theta_min': {
                    'widget': GUI_custom_components.Label(text='Theta min'),
                    'widget_type': 'lbl',
                    'value_attribute': 'vertex.theta_min'},
                'lbl_theta_max': {
                    'widget': GUI_custom_components.Label(text='Theta max'),
                    'widget_type': 'lbl',
                    'value_attribute': 'vertex.theta_max'},
                'lbl_theta_mean': {
                    'widget': GUI_custom_components.Label(text='Theta mean'),
                    'widget_type': 'lbl',
                    'value_attribute': 'vertex.theta_angle_mean'},
                'lbl_theta_angle_variance': {
                    'widget': GUI_custom_components.Label(text='Theta variance'),
                    'widget_type': 'lbl',
                    'value_attribute': 'vertex.theta_angle_variance'},
                'chb_precipitate_column': {
                    'widget': GUI_custom_components.CheckBox(label='Precipitate column', trigger_func=self.chb_precipitate_column_trigger),
                    'widget_type': 'chb',
                    'value_attribute': 'vertex.is_in_precipitate'},
                'chb_show_in_overlay': {
                    'widget': GUI_custom_components.CheckBox(label='Show in overlay', trigger_func=self.chb_show_in_overlay_trigger),
                    'widget_type': 'chb',
                    'value_attribute': 'vertex.show_in_overlay'},
                'mchb_move': {
                    'widget': GUI_custom_components.MoveControls(text='Enable move', trigger_chb=self.mchb_move_trigger, trigger_1=self.btn_set_position_trigger, trigger_2=self.btn_cancel_move_trigger),
                    'widget_type': 'mchb'},
                'btn_delete_column': {
                    'widget': GUI_custom_components.MediumButton('Delete', self, trigger_func=self.btn_delete_column_trigger),
                    'widget_type': 'btn'},
                'btn_print_details': {
                    'widget': GUI_custom_components.MediumButton('Print', self, trigger_func=self.btn_print_details_trigger),
                    'widget_type': 'btn'},
                'btn_snap': {
                    'widget': GUI_custom_components.MediumButton('Snap', self, trigger_func=self.btn_snap_trigger),
                    'widget_type': 'btn'},
                'btn_deselect': {
                    'widget': GUI_custom_components.MediumButton('Deselect', self, trigger_func=self.btn_deselect_trigger),
                    'widget_type': 'btn'},
                'btn_new_column': {
                    'widget': GUI_custom_components.MediumButton('New', self, trigger_func=self.btn_new_column_trigger),
                    'widget_type': 'btn'}
            }, 'atomic_graph': {
                'lbl_order': {
                    'widget': GUI_custom_components.Label(text='Order'),
                    'widget_type': 'lbl',
                    'value_attribute': 'project_instance.graph.order'},
                'lbl_size': {
                    'widget': GUI_custom_components.Label(text='Size'),
                    'widget_type': 'lbl',
                    'value_attribute': 'project_instance.graph.size'},
                'lbl_volume': {
                    'widget': GUI_custom_components.Label(text='Volume'),
                    'widget_type': 'lbl',
                    'value_attribute': 'project_instance.graph.volume'},
                'lbl_chi': {
                    'widget': GUI_custom_components.Label(text='Chi'),
                    'widget_type': 'lbl',
                    'value_attribute': 'project_instance.graph.chi'},
                'lbl_zeta_chi': {
                    'widget': GUI_custom_components.Label(text='Zeta-chi'),
                    'widget_type': 'lbl',
                    'value_attribute': 'project_instance.graph.zeta_chi'},
                'lbl_average_degree': {
                    'widget': GUI_custom_components.Label(text='Average degree'),
                    'widget_type': 'lbl',
                    'value_attribute': 'project_instance.graph.avg_degree'},
                'lbl_average_zeta_degree': {
                    'widget': GUI_custom_components.Label(text='Average zeta degree'),
                    'widget_type': 'lbl',
                    'value_attribute': 'project_instance.graph.avg_zeta_degree'},
                'chb_permute_mode': {
                    'widget': GUI_custom_components.CheckBox(label='Enable permute mode'),
                    'widget_type': 'chb'},
                'chb_enable_ruler': {
                    'widget': GUI_custom_components.CheckBox(label='Enable ruler'),
                    'widget_type': 'chb'},
                'chb_show_mesh_vertex_order': {
                    'widget': GUI_custom_components.CheckBox(label='Show mesh vertex order', trigger_func=self.chb_show_mesh_vertex_order),
                    'widget_type': 'chb'},
                'chb_show_zeta_0': {
                    'widget': GUI_custom_components.CheckBox(label='Show zeta=0 anti-graph', trigger_func=self.chb_show_zeta_0_trigger, default_state=True),
                    'widget_type': 'chb'},
                'chb_show_zeta_1': {
                    'widget': GUI_custom_components.CheckBox(label='Show zeta=1 anti-graph', trigger_func=self.chb_show_zeta_1_trigger, default_state=True),
                    'widget_type': 'chb'},
                'btn_refresh_graph': {
                    'widget': GUI_custom_components.MediumButton('Refresh', self,  trigger_func=self.btn_refresh_graph_trigger),
                    'widget_type': 'btn'},
                'btn_refresh_mesh': {
                    'widget': GUI_custom_components.MediumButton('Refresh mesh', self, trigger_func=self.btn_refresh_mesh_trigger),
                    'widget_type': 'btn'},
                'btn_print_hsd': {
                    'widget': GUI_custom_components.MediumButton('Print H.S.D.', self, trigger_func=self.btn_print_hsd_trigger),
                    'widget_type': 'btn'},
                'btn_gen_sub': {
                    'widget': GUI_custom_components.MediumButton('Build sub-graph', self, trigger_func=self.btn_gen_sub_trigger),
                    'widget_type': 'btn'},
                'btn_build_anti_graph': {
                    'widget': GUI_custom_components.MediumButton('Build anti-graph', self, trigger_func=self.btn_build_anti_graph_trigger),
                    'widget_type': 'btn'},
                'btn_match_zeta': {
                    'widget': GUI_custom_components.MediumButton('Match zeta', self, trigger_func=self.btn_match_zeta_trigger),
                    'widget_type': 'btn'}
            }, 'heat_maps': {
                'sbtn_current_map': {
                    'widget': GUI_custom_components.SetButtonLayout(text='Currently displayed heat-map', value='None', obj=self, trigger_func=self.sbtn_current_map_trigger),
                    'widget_type': 'sbtn'},
                'btn_calc_heat': {
                    'widget': GUI_custom_components.MediumButton('Make heat-map', self, trigger_func=self.btn_calc_heat_trigger),
                    'widget_type': 'btn'},
                'btn_delete_heat': {
                    'widget': GUI_custom_components.MediumButton('Delete heat-map', self, trigger_func=self.btn_delete_heat_trigger),
                    'widget_type': 'btn'},
                'chb_legend': {
                    'widget': GUI_custom_components.CheckBox(label='legend', trigger_func=self.chb_heat_legend_trigger, default_state=True),
                    'widget_type': 'chb'}
            }, 'data': {
                'btn_project_plots': {
                    'widget': GUI_custom_components.MediumButton('Project plots', self, trigger_func=self.btn_project_plots_trigger),
                    'widget_type': 'btn'},
                'btn_model_plots': {
                    'widget': GUI_custom_components.MediumButton('Model plots', self, trigger_func=self.btn_model_plots_trigger),
                    'widget_type': 'btn'},
                'btn_calc_models': {
                    'widget': GUI_custom_components.MediumButton('Calculate model', self, trigger_func=self.btn_calc_models_trigger),
                    'widget_type': 'btn'}
            }, 'overlay': {
                'btn_configure': {
                    'widget': GUI_custom_components.MediumButton('Configure', self, trigger_func=self.btn_configure_trigger),
                    'widget_type': 'btn'},
                'btn_show_all': {
                    'widget': GUI_custom_components.MediumButton('Show all', self, trigger_func=self.btn_show_all_trigger),
                    'widget_type': 'btn'},
                'btn_hide_all': {
                    'widget': GUI_custom_components.MediumButton('Hide all', self, trigger_func=self.btn_hide_all_trigger),
                    'widget_type': 'btn'},
                'chb_raw_image': {
                    'widget': GUI_custom_components.CheckBox(label='Image', trigger_func=self.chb_raw_image_trigger, default_state=True),
                    'widget_type': 'chb'},
                'chb_black_background': {
                    'widget': GUI_custom_components.CheckBox(label='Black background', trigger_func=self.chb_black_background_trigger, default_state=True),
                    'widget_type': 'chb'},
                'chb_edge_columns': {
                    'widget': GUI_custom_components.CheckBox(label='Edge columns', trigger_func=self.chb_edge_columns_trigger, default_state=True),
                    'widget_type': 'chb'},
                'chb_precipitate': {
                    'widget': GUI_custom_components.CheckBox(label='Precipitate', trigger_func=self.chb_precipitate_trigger, default_state=True),
                    'widget_type': 'chb'},
                'chb_matrix': {
                    'widget': GUI_custom_components.CheckBox(label='Matrix', trigger_func=self.chb_matrix_trigger, default_state=True),
                    'widget_type': 'chb'},
                'chb_legend': {
                    'widget': GUI_custom_components.CheckBox(label='Legend', trigger_func=self.chb_legend_trigger, default_state=True),
                    'widget_type': 'chb'},
                'chb_scalebar': {
                    'widget': GUI_custom_components.CheckBox(label='Scalebar', trigger_func=self.chb_scalebar_trigger, default_state=True),
                    'widget_type': 'chb'},
                'chb_zeta_0': {
                    'widget': GUI_custom_components.CheckBox(label='Zeta=0 plane', trigger_func=self.chb_zeta_0_trigger, default_state=True),
                    'widget_type': 'chb'},
                'chb_zeta_1': {
                    'widget': GUI_custom_components.CheckBox(label='Zeta=1 plane', trigger_func=self.chb_zeta_1_trigger, default_state=True),
                    'widget_type': 'chb'}
            }, 'debug': {
                'sbtn_starting_index': {
                    'widget': GUI_custom_components.SetButtonLayout(obj=self, trigger_func=self.sbtn_starting_index_trigger, text='Default starting index'),
                    'widget_type': 'sbtn',
                    'value_attribute': 'project_instance.starting_index'},
                'btn_test': {
                    'widget': GUI_custom_components.MediumButton('Test', self, trigger_func=self.btn_test_trigger),
                    'widget_type': 'btn'},
                'btn_crash': {
                    'widget': GUI_custom_components.MediumButton('Crash program', self, trigger_func=self.btn_crash_trigger),
                    'widget_type': 'btn'},
                'btn_test_data_manager': {
                    'widget': GUI_custom_components.MediumButton('Test data-manager', self, trigger_func=self.btn_test_data_manager_trigger),
                    'widget_type': 'btn'},
                'btn_show_prediction': {
                    'widget': GUI_custom_components.MediumButton('Show model prediction', self, trigger_func=self.btn_show_prediction_trigger),
                    'widget_type': 'btn'},
                'btn_check_integrity': {
                    'widget': GUI_custom_components.MediumButton('Check integrity', self, trigger_func=self.btn_check_integrity_trigger),
                    'widget_type': 'btn'}
            }
        }

        # Add dynamic widgets:
        if self.ui_obj.project_instance is None:
            dict_ = GUI.core.Project.default_species_dict
        else:
            dict_ = self.ui_obj.project_instance.species_dict
        for species in dict_[self.ui_obj.overlay_settings['display_mode']]:
            self.widgets['overlay']['chb_{}'.format(species)] = {
                'widget': GUI_custom_components.CheckBox(label='{}'.format(species), trigger_func=self.chb_dynamic_trigger, default_state=True),
                'widget_type': 'chb'
            }
            GUI_tooltips.tooltips['overlay']['chb_{}'.format(species)] = 'Toggle {} columns'.format(species)
        self.species_states = {}
        for species in dict_[self.ui_obj.overlay_settings['display_mode']]:
            self.species_states['chb_{}'.format(species)] = self.widgets['overlay']['chb_{}'.format(species)]['widget'].isChecked()

        # -------------------------------------------------------------------------
        # Layouts
        # -------------------------------------------------------------------------

        self.groups = {}

        for group in self.widgets:
            self.groups[group] = GUI_custom_components.GroupBox(group.replace('_', ' '))
            num_buttons = 0
            btn_layout = QtWidgets.QVBoxLayout()
            btn_inner_layouts = []
            widget_layout = QtWidgets.QVBoxLayout()
            combined_layout = QtWidgets.QVBoxLayout()
            for widget in self.widgets[group]:
                if self.widgets[group][widget]['widget_type'] == 'btn':
                    if np.mod(num_buttons, 3) == 0:
                        btn_inner_layouts.append(QtWidgets.QHBoxLayout())
                    btn_inner_layouts[-1].addWidget(self.widgets[group][widget]['widget'])
                    num_buttons += 1
                else:
                    if self.widgets[group][widget]['widget_type'] == 'sbtn' or self.widgets[group][widget]['widget_type'] == 'mchb':
                        widget_layout.addLayout(self.widgets[group][widget]['widget'])
                    else:
                        widget_layout.addWidget(self.widgets[group][widget]['widget'])
            for inner_layout in btn_inner_layouts:
                inner_layout.addStretch()
                btn_layout.addLayout(inner_layout)
            combined_layout.addLayout(btn_layout)
            combined_layout.addLayout(widget_layout)
            self.groups[group].setLayout(combined_layout)

        top_layout = QtWidgets.QVBoxLayout()
        for group in self.groups:
            top_layout.addWidget(self.groups[group])

        self.setLayout(top_layout)

        # Set tooltips:
        self.mode_tooltip(self.ui_obj.menu.toggle_tooltips_action.isChecked())

    def mode_tooltip(self, on):
        if on:
            for group in self.widgets:
                for widget in self.widgets[group]:
                    if widget in GUI_tooltips.tooltips[group]:
                        self.widgets[group][widget]['widget'].setToolTip(GUI_tooltips.tooltips[group][widget])
        else:
            for group in self.widgets:
                for widget in self.widgets[group]:
                    if widget in GUI_tooltips.tooltips[group]:
                        self.widgets[group][widget]['widget'].setToolTip('')

    def mode_move(self, on):
        if self.ui_obj.project_instance is not None and not self.ui_obj.selected_column == -1:
            self.widgets['selected_column']['mchb_move']['widget'].blockSignals(True)
            self.widgets['selected_column']['mchb_move']['widget'].setChecked(on)
            self.widgets['selected_column']['mchb_move']['widget'].blockSignals(False)
            for i, pos_obj in enumerate(self.ui_obj.gs_atomic_positions.interactive_position_objects):
                if not i == self.ui_obj.selected_column:
                    pos_obj.setFlag(QtWidgets.QGraphicsItem.ItemIsSelectable, not on)
                else:
                    pos_obj.setFlag(QtWidgets.QGraphicsItem.ItemIsMovable, on)
                self.ui_obj.gs_overlay_composition.interactive_overlay_objects[i].setFlag(QtWidgets.QGraphicsItem.ItemIsSelectable, not on)
                self.ui_obj.gs_atomic_graph.interactive_vertex_objects[i].setFlag(QtWidgets.QGraphicsItem.ItemIsSelectable, not on)
            for group in self.widgets:
                for tag in self.widgets[group]:
                    if not tag == 'mchb_move':
                        self.widgets[group][tag]['widget'].setDisabled(on)
            self.widgets['selected_column']['mchb_move']['widget'].setDisabled(not on)

    def draw_histogram(self):

        box_width = 15  # bin width
        box_seperation = 10  # bin separation
        box_displacement = 25  # width + separation

        species = set()
        colors = {}
        if self.ui_obj.project_instance is None:
            species.add('Si')
            colors['Si'] = (255, 0, 0)
            species.add('Cu')
            colors['Cu'] = (255, 255, 0)
            species.add('Al')
            colors['Al'] = (0, 255, 0)
            species.add('Mg')
            colors['Mg'] = (138, 43, 226)
            species.add('Un')
            colors['Un'] = (0, 0, 255)
        else:
            for key, value in self.ui_obj.project_instance.species_dict[self.ui_obj.overlay_settings['display_mode']].items():
                species.add(key)
                colors[key] = value['color']

        species = list(species)
        species.sort()
        box_heights = []
        if not self.ui_obj.selected_column == -1:
            for atomic_species in species:
                if self.ui_obj.overlay_settings['display_mode'] == 'atomic_species':
                    box_heights.append(int(100 * self.ui_obj.project_instance.graph.vertices[self.ui_obj.selected_column].probability_vector[atomic_species]))
                else:
                    box_heights.append(int(100 * self.ui_obj.project_instance.graph.vertices[self.ui_obj.selected_column].advanced_probability_vector[atomic_species]))

        else:
            box_heights = [0.0] * len(species)

        width = len(species) * box_width + (len(species) + 1) * box_seperation
        self.probGraphicView.setMinimumWidth(width)
        self.probGraphicView.setMinimumWidth(width)
        box = QtWidgets.QGraphicsRectItem(0, 0, width, self.height)
        box.setPen(GUI_settings.pen_boarder)
        box.hide()

        probGraphicScene = QtWidgets.QGraphicsScene()
        probGraphicScene.addItem(box)

        for i, box_height in enumerate(box_heights):
            x_box = box_seperation + i * box_displacement
            species_box = QtWidgets.QGraphicsRectItem(0, 0, box_width, box_height)
            species_box.setBrush(QtGui.QBrush(QtGui.QColor(
                colors[species[i]][0],
                colors[species[i]][1],
                colors[species[i]][2]
            )))
            species_box.setX(x_box)
            species_box.setY(100 - box_height + 10)
            species_text = QtWidgets.QGraphicsSimpleTextItem()
            species_text.setText(species[i])
            species_text.setFont(GUI_settings.font_tiny)
            species_text.setX(x_box + box_width / 2 - species_text.boundingRect().width() / 2)
            species_text.setY(100 + 4 + 10)
            species_number = QtWidgets.QGraphicsSimpleTextItem()
            species_number.setText(str(box_height / 100))
            species_number.setFont(GUI_settings.font_tiny)
            species_number.setX(x_box - 1)
            species_number.setY(100 - box_height)

            probGraphicScene.addItem(species_box)
            probGraphicScene.addItem(species_text)
            probGraphicScene.addItem(species_number)

        if GUI_settings.theme == 'dark':
            probGraphicScene.palette().setColor(QtGui.QPalette.Text, QtCore.Qt.white)

        self.probGraphicView.setScene(probGraphicScene)

    def select_column(self):

        i = self.ui_obj.selected_column

        if i == -1:
            self.deselect_column()
        else:

            vertex = self.ui_obj.project_instance.graph.vertices[i]

            for key, widget_dict in self.widgets['selected_column'].items():
                if widget_dict['widget_type'] == 'sbtn' or widget_dict['widget_type'] == 'lbl':
                    self.widgets['selected_column'][key]['widget'].set_value(
                        eval('{}'.format(widget_dict['value_attribute']))
                    )
                elif widget_dict['widget_type'] == 'chb':
                    self.widgets['selected_column'][key]['widget'].set_state_no_trigger(
                        eval('{}'.format(widget_dict['value_attribute']))
                    )

            self.draw_histogram()

    def deselect_column(self):

        self.draw_histogram()

        for key, widget_dict in self.widgets['selected_column'].items():
            if widget_dict['widget_type'] == 'sbtn' or widget_dict['widget_type'] == 'lbl':
                self.widgets['selected_column'][key]['widget'].set_value(
                    ''
                )
            elif widget_dict['widget_type'] == 'chb':
                self.widgets['selected_column'][key]['widget'].set_state_no_trigger(
                    eval('{}'.format(widget_dict['widget'].default_state))
                )

    def update_display(self):

        if self.ui_obj.project_instance is not None:

            self.widgets['project']['sbtn_active_model']['widget'].set_value(self.ui_obj.project_instance.graph.active_model)
            if self.ui_obj.savefile:
                self.widgets['project']['lbl_filename']['widget'].set_value(os.path.split(self.ui_obj.savefile)[1])
                self.widgets['project']['lbl_location']['widget'].set_value(os.path.split(self.ui_obj.savefile)[0])
            else:
                self.widgets['project']['lbl_filename']['widget'].set_value('')
                self.widgets['project']['lbl_location']['widget'].set_value('')
            self.widgets['heat_maps']['sbtn_current_map']['widget'].set_value(self.ui_obj.gs_heat.title)

            for group in self.widgets:
                if not group == 'project' and not group == 'selected_column':
                    for key, widget_dict in self.widgets[group].items():
                        if 'value_attribute' in widget_dict and not widget_dict['widget_type'] == 'chb':
                            self.widgets[group][key]['widget'].set_value(
                                eval('self.ui_obj.{}'.format(widget_dict['value_attribute']))
                            )

            if self.ui_obj.selected_column == -1:
                self.deselect_column()
            else:
                self.select_column()

        else:

            self.empty_display()

    def empty_display(self):

        self.deselect_column()

        self.widgets['project']['sbtn_active_model']['widget'].set_value('default_model')
        self.widgets['project']['lbl_filename']['widget'].set_value('')
        self.widgets['project']['lbl_location']['widget'].set_value('')

        for group in self.widgets:
            if not group == 'project' and not group == 'selected_column':
                for key, widget_dict in self.widgets[group].items():
                    if 'value_attribute' in widget_dict and not widget_dict['widget_type'] == 'chb':
                        self.widgets[group][key]['widget'].set_value('')

    def keyPressEvent(self, event):

        self.ui_obj.keyPressEvent(event)

    def sbtn_active_model_trigger(self):
        if self.ui_obj.project_instance is not None:
            model_filename = QtWidgets.QFileDialog.getOpenFileName(self, "Load model", '', "")
            if model_filename[0]:
                self.ui_obj.project_instance.graph.active_model = model_filename[0]
            self.widgets['project']['sbtn_active_model']['widget'].set_value(self.ui_obj.project_instance.graph.active_model)

    def btn_import_trigger(self):
        ImportWizard(ui_obj=self.ui_obj)

    def btn_open_project_trigger(self):
        self.ui_obj.menu_open_trigger()

    def btn_save_project_trigger(self):
        if self.ui_obj.project_instance is not None:
            self.ui_obj.menu_save_trigger()

    def btn_close_project_trigger(self):
        if self.ui_obj.project_instance is not None:
            self.ui_obj.menu_close_trigger()

    def btn_species_dict_trigger(self):
        EditSpeciesDict(ui_obj=self.ui_obj)

    def btn_export_trigger(self):
        ExportWizard(ui_obj=self.ui_obj)

    def chb_lock_views_trigger(self, state):
        if state:
            self.btn_align_views_trigger()
            self.ui_obj.lock_views = True
        else:
            self.ui_obj.lock_views = False

    def btn_show_stats_trigger(self):
        if self.ui_obj.project_instance is not None:
            self.ui_obj.project_instance.report(supress_log=False)

    def btn_show_source_trigger(self):
        if self.ui_obj.project_instance is not None:
            logger.info(self.ui_obj.project_instance.filename_full)

    def btn_align_views_trigger(self):
        tab = self.ui_obj.tabs.currentIndex()
        coor = self.ui_obj.gv_list[tab].mapToScene(self.ui_obj.gv_list[tab].viewport().rect().center()) / self.ui_obj.gv_list[tab].scaling_factor
        transform = self.ui_obj.gv_list[tab].transform()
        transform.setMatrix(
            transform.m11() * self.ui_obj.gv_list[tab].scaling_factor,
            transform.m12(),
            transform.m13(),
            transform.m21(),
            transform.m22() * self.ui_obj.gv_list[tab].scaling_factor,
            transform.m23(),
            transform.m31(),
            transform.m32(),
            transform.m33()
        )
        for i, gv in enumerate(self.ui_obj.gv_list):
            gv.resetTransform()
            transform.setMatrix(
                transform.m11() / gv.scaling_factor,
                transform.m12(),
                transform.m13(),
                transform.m21(),
                transform.m22() / gv.scaling_factor,
                transform.m23(),
                transform.m31(),
                transform.m32(),
                transform.m33()
            )
            gv.setTransform(transform)
            transform.setMatrix(
                transform.m11() * gv.scaling_factor,
                transform.m12(),
                transform.m13(),
                transform.m21(),
                transform.m22() * gv.scaling_factor,
                transform.m23(),
                transform.m31(),
                transform.m32(),
                transform.m33()
            )
            gv.centerOn(coor * gv.scaling_factor)

    def chb_toggle_positions_trigger(self, state):
        if self.ui_obj.project_instance is not None:
            if self.ui_obj.project_instance.num_columns > 0:
                for position_graphic in self.ui_obj.gs_atomic_positions.interactive_position_objects:
                    if state:
                        position_graphic.show()
                    else:
                        position_graphic.hide()

    def sbtn_threshold_trigger(self):
        if self.ui_obj.project_instance is not None:
            threshold, ok_pressed = QtWidgets.QInputDialog.getDouble(self, "Set", "Threshold value (decimal between 0 and 1):", self.ui_obj.project_instance.threshold, 0, 1, 5)
            if ok_pressed:
                self.ui_obj.project_instance.threshold = threshold
                self.widgets['column_detection']['sbtn_threshold']['widget'].set_value(self.ui_obj.project_instance.threshold)

    def sbtn_search_size_trigger(self):
        if self.ui_obj.project_instance is not None:
            search_size, ok_pressed = QtWidgets.QInputDialog.getInt(self, "Set", "Search size:",
                                                                    self.ui_obj.project_instance.search_size, 0, 100000, 100)
            if ok_pressed:
                self.ui_obj.project_instance.search_size = search_size
                self.widgets['column_detection']['sbtn_search_size']['widget'].set_value(self.ui_obj.project_instance.search_size)

    def sbtn_scale_trigger(self):
        if self.ui_obj.project_instance is not None:
            scale, ok_pressed = QtWidgets.QInputDialog.getDouble(self, "Set", "Image scale (pm/pixel):", self.ui_obj.project_instance.scale, 0, 10000, 4)
            if ok_pressed:
                self.ui_obj.sys_message('Working...')
                self.ui_obj.project_instance.scale = scale
                self.ui_obj.project_instance.r = int(100 / scale)
                self.ui_obj.project_instance.overhead = int(6 * (self.ui_obj.project_instance.r / 10))
                self.widgets['column_detection']['sbtn_scale']['widget'].set_value(self.ui_obj.project_instance.scale)
                self.widgets['column_detection']['lbl_atomic_radii']['widget'].set_value(self.ui_obj.project_instance.r)
                self.widgets['column_detection']['lbl_overhead']['widget'].set_value(self.ui_obj.project_instance.overhead)
                self.ui_obj.project_instance.redraw_search_mat()
                self.ui_obj.project_instance.graph.scale = scale
                for vertex in self.ui_obj.project_instance.graph.vertices:
                    vertex.r = self.ui_obj.project_instance.r
                    vertex.scale = scale
                    vertex.spatial_coor_x = vertex.im_coor_x * scale
                    vertex.spatial_coor_y = vertex.im_coor_y * scale
                    vertex.spatial_coor_z = vertex.zeta * 0.5 * vertex.al_lattice_const
                self.ui_obj.update_display()
                self.ui_obj.sys_message('Ready.')

    def btn_start_column_detection_trigger(self):
        if self.ui_obj.project_instance is not None:
            items = ('Threshold', 'Search size', 'other')
            item, ok_pressed = QtWidgets.QInputDialog.getItem(self, "Set", "Search type", items, 0, False)
            if ok_pressed and item:
                self.ui_obj.sys_message('Working, column detection might take several minutes...')
                if item == 'Search size':
                    self.ui_obj.project_instance.column_detection('s', plot=self.widgets['column_detection']['chb_plot_results']['widget'].isChecked())
                elif item == 'Threshold':
                    self.ui_obj.project_instance.column_detection('t', plot=self.widgets['column_detection']['chb_plot_results']['widget'].isChecked())
                else:
                    self.ui_obj.project_instance.column_detection('o', plot=self.widgets['column_detection']['chb_plot_results']['widget'].isChecked())
                self.ui_obj.update_display()
                self.ui_obj.sys_message('Ready.')

    def btn_reset_column_detection_trigger(self):
        if self.ui_obj.project_instance is not None:
            self.ui_obj.project_instance.reset_graph()

    def btn_redraw_search_mat_trigger(self):
        if self.ui_obj.project_instance is not None:
            self.ui_obj.sys_message('Working...')
            self.ui_obj.project_instance.redraw_search_mat()
            self.ui_obj.update_search_matrix()
            self.ui_obj.sys_message('Ready.')

    def btn_start_column_characterization_trigger(self):
        if self.ui_obj.project_instance is not None:

            strings = [
                '0 - Full column characterization algorithm',
                '1 - Basic mappings...',
                '2 - ...The rest',
                '3 - Spatial mapping',
                '4 - Identify edge columns',
                '5 - Not in use',
                '6 - Zeta analysis',
                '7 - Apply alpha model',
                '8 - Particle detection',
                '9 - Calculate normalized gamma',
                '10 - Evaluate sub-species',
                '11 - Apply composite model',
                '12 - Reset probability vectors',
                '13 - Reset user-set columns',
                '14 - Search for intersections',
                '15 - Not in use',
                '16 - Map vertex connectivity',
                '17 - Map vertex in-connectivity',
                '18 - Untangling algorithm',
                '19 - Weak untangling',
                '20 - Not in use',
                '21 - Not in use',
                '22 - Not in use',
                '23 - Plot gamma'
            ]

            string, ok_pressed = QtWidgets.QInputDialog.getItem(self, "Set", "Search step", strings, 0, False)
            if ok_pressed and strings:
                self.ui_obj.sys_message('Analyzing... This may take a long time...')
                choice = -1
                for k in range(0, len(strings)):
                    if string == strings[k]:
                        choice = k
                if not choice == -1:
                    if self.ui_obj.selected_column == -1:
                        if self.ui_obj.project_instance.starting_index is not None:
                            starting_column = self.ui_obj.project_instance.starting_index
                        else:
                            starting_column = 0
                    else:
                        starting_column = self.ui_obj.selected_column

                    if self.widgets['column_characterization']['chb_show_graphics']['widget'].isChecked():
                        self.ui_obj.project_instance.column_characterization(starting_column, choice, ui_obj=self)
                    else:
                        self.ui_obj.project_instance.column_characterization(starting_column, choice)
                    self.ui_obj.update_display()
                else:
                    logger.error('Invalid selection. Was not able to start column detection.')

    def btn_reset_column_characterization_trigger(self):
        pass

    def btn_invert_zeta_trigger(self):
        if self.ui_obj.project_instance is not None:
            self.ui_obj.project_instance.graph.invert_levels()
            self.ui_obj.update_central_widget()
            self.select_column()

    def chb_precipitate_column_trigger(self, state):
        if self.ui_obj.project_instance is not None and not self.ui_obj.selected_column == -1:
            self.ui_obj.project_instance.graph.vertices[self.ui_obj.selected_column].is_in_precipitate = state

    def chb_show_in_overlay_trigger(self, state):
        if self.ui_obj.project_instance is not None and not self.ui_obj.selected_column == -1:
            self.ui_obj.project_instance.graph.vertices[self.ui_obj.selected_column].show_in_overlay = state
            self.ui_obj.gs_atomic_positions.interactive_position_objects[self.ui_obj.selected_column].set_style()
            self.ui_obj.gs_overlay_composition.interactive_overlay_objects[self.ui_obj.selected_column].set_style()

    def mchb_move_trigger(self, state):
        self.ui_obj.sys_message('Working...')
        self.mode_move(state)
        self.ui_obj.sys_message('Ready.')

    def sbtn_find_column_trigger(self):
        if self.ui_obj.project_instance is not None:
            index, ok_pressed = QtWidgets.QInputDialog.getInt(self, "Set", "Find column by index:", 0, 0, 100000, 1)
            if ok_pressed:
                if index < self.ui_obj.project_instance.num_columns:
                    self.ui_obj.gs_atomic_positions.interactive_position_objects[index].mouseReleaseEvent(
                        QtWidgets.QGraphicsEllipseItem.mouseReleaseEvent)

    def sbtn_atomic_species_trigger(self):
        if self.ui_obj.project_instance is not None and not self.ui_obj.selected_column == -1:
            items = []
            for item in self.ui_obj.project_instance.species_dict['atomic_species']:
                items.append(item)
            item, ok_pressed = QtWidgets.QInputDialog.getItem(self, "Set", "Species", items, 0, False)
            if ok_pressed and item:
                self.ui_obj.set_atomic_species(item)
                self.ui_obj.project_instance.graph.vertices[self.ui_obj.selected_column].is_set_by_user = True

    def sbtn_advanced_species_trigger(self):
        if self.ui_obj.project_instance is not None and not self.ui_obj.selected_column == -1:
            items = []
            for item in self.ui_obj.project_instance.species_dict['advanced_species']:
                items.append(item)
            item, ok_pressed = QtWidgets.QInputDialog.getItem(self, "Set", "Species", items, 0, False)
            if ok_pressed and item:
                self.ui_obj.set_advanced_species(item)
                self.ui_obj.project_instance.graph.vertices[self.ui_obj.selected_column].is_set_by_user = True

    def sbtn_zeta_trigger(self):
        if self.ui_obj.project_instance is not None and not self.ui_obj.selected_column == -1:
            items = ('0', '1')
            item, ok_pressed = QtWidgets.QInputDialog.getItem(self, "Set", "zeta", items, 0, False)
            if ok_pressed and item:
                if item == '0':
                    zeta = 0
                elif item == '1':
                    zeta = 1
                else:
                    zeta = 0
                self.ui_obj.set_zeta(zeta)

    def btn_cancel_move_trigger(self):
        self.mode_move(False)
        self.ui_obj.update_central_widget()

    def btn_set_position_trigger(self):
        if self.ui_obj.project_instance is not None:
            self.ui_obj.sys_message('Working...')
            (x, y) = self.ui_obj.gs_atomic_positions.interactive_position_objects[self.ui_obj.selected_column].get_vertex_pos()
            vertex = self.ui_obj.project_instance.graph.vertices[self.ui_obj.selected_column]
            vertex.im_coor_x = x
            vertex.im_coor_y = y
            vertex.spatial_coor_x = x * self.ui_obj.project_instance.scale
            vertex.spatial_coor_y = y * self.ui_obj.project_instance.scale
            vertex.avg_gamma, vertex.peak_gamma = utils.circular_average(
                self.ui_obj.project_instance.im_mat,
                int(x),
                int(y),
                self.ui_obj.project_instance.r
            )
            self.ui_obj.project_instance.graph.vertex_moved(vertex.i)
            self.mode_move(False)
            self.ui_obj.update_central_widget()
            self.ui_obj.sys_message('Ready.')

    def btn_delete_column_trigger(self):
        if self.ui_obj.project_instance is not None and not self.ui_obj.selected_column == -1:
            ok_pressed = QtWidgets.QMessageBox.question(self, 'Confirm', 'Are you sure you wish to delete this column?',
                                                        QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                                                        QtWidgets.QMessageBox.No)
            if ok_pressed == QtWidgets.QMessageBox.Yes:
                self.ui_obj.sys_message('Working...')
                i = self.ui_obj.selected_column
                self.btn_deselect_trigger()
                self.ui_obj.project_instance.graph.remove_vertex(i)
                self.ui_obj.project_instance.num_columns -= 1
                self.ui_obj.update_central_widget()
                logger.info('Vertex removed. remember to redraw the search matrix before continuing column detection!')
                self.ui_obj.sys_message('Ready.')

    def btn_print_details_trigger(self):
        if self.ui_obj.project_instance is not None and not self.ui_obj.selected_column == -1:
            self.ui_obj.project_instance.vertex_report(self.ui_obj.selected_column)

    def btn_snap_trigger(self):
        if self.ui_obj.project_instance is not None and not self.ui_obj.selected_column == -1:
            i = self.ui_obj.selected_column
            tab = self.ui_obj.tabs.currentIndex()
            coor = QtCore.QPointF(self.ui_obj.project_instance.graph.vertices[i].im_coor_x, self.ui_obj.project_instance.graph.vertices[i].im_coor_y)
            transform = QtGui.QTransform()
            transform.scale(7.5 * self.ui_obj.gv_list[tab].scaling_factor, 7.5 * self.ui_obj.gv_list[tab].scaling_factor)
            for i, gv in enumerate(self.ui_obj.gv_list):
                gv.resetTransform()
                transform.setMatrix(
                    transform.m11() / gv.scaling_factor,
                    transform.m12(),
                    transform.m13(),
                    transform.m21(),
                    transform.m22() / gv.scaling_factor,
                    transform.m23(),
                    transform.m31(),
                    transform.m32(),
                    transform.m33()
                )
                gv.setTransform(transform)
                transform.setMatrix(
                    transform.m11() * gv.scaling_factor,
                    transform.m12(),
                    transform.m13(),
                    transform.m21(),
                    transform.m22() * gv.scaling_factor,
                    transform.m23(),
                    transform.m31(),
                    transform.m32(),
                    transform.m33()
                )
                gv.centerOn(coor * gv.scaling_factor)

    def btn_deselect_trigger(self):
        self.ui_obj.column_selected(-1)

    def btn_new_column_trigger(self):
        if self.ui_obj.project_instance is not None:
            self.ui_obj.sys_message('Working...')
            new_vertex = GUI.core.graph_2.Vertex(
                self.ui_obj.project_instance.graph.order,
                self.ui_obj.project_instance.im_width / 2 + 0.1,
                self.ui_obj.project_instance.im_height / 2 + 0.1,
                self.ui_obj.project_instance.r,
                self.ui_obj.project_instance.scale,
                parent_graph=self.ui_obj.project_instance.graph
            )
            self.ui_obj.project_instance.graph.add_vertex(new_vertex)
            self.ui_obj.project_instance.num_columns += 1
            index = new_vertex.i
            self.ui_obj.update_central_widget()
            self.ui_obj.column_selected(index)
            self.ui_obj.sys_message('Setting move mode to on...')
            self.mode_move(True)

    def chb_show_mesh_vertex_order(self):
        if self.ui_obj.project_instance is not None:
            self.ui_obj.sys_message('Working...')
            self.ui_obj.gs_atomic_graph.re_draw_mesh_details()
            self.ui_obj.sys_message('Ready.')

    def chb_show_zeta_0_trigger(self, state):
        if self.ui_obj.project_instance is not None:
            if self.ui_obj.gs_anti_graph.graph is not None:
                self.ui_obj.gs_anti_graph.toggle_level_0(state)

    def chb_show_zeta_1_trigger(self, state):
        if self.ui_obj.project_instance is not None:
            if self.ui_obj.gs_anti_graph.graph is not None:
                self.ui_obj.gs_anti_graph.toggle_level_1(state)

    def btn_refresh_graph_trigger(self):
        if self.ui_obj.project_instance is not None:
            self.ui_obj.sys_message('Working...')
            self.ui_obj.project_instance.graph.refresh_graph()
            self.ui_obj.update_display()
            self.ui_obj.sys_message('Ready.')

    def btn_refresh_mesh_trigger(self):
        if self.ui_obj.project_instance is not None:
            self.ui_obj.sys_message('Working...')
            self.ui_obj.project_instance.graph.map_meshes()
            self.ui_obj.update_display()
            self.ui_obj.sys_message('Ready.')

    def btn_print_hsd_trigger(self):
        if self.ui_obj.project_instance is not None:
            interatomicHSD = []
            for atomic_species_1, details_1 in self.ui_obj.project_instance.species_dict['atomic_species'].items():
                for atomic_species_2, details_2 in self.ui_obj.project_instance.species_dict['atomic_species'].items():
                    interatomicHSD.append({
                        'pair': '{}<=>{}'.format(atomic_species_1, atomic_species_2),
                        'HSD': details_1['atomic_radii'] + details_2['atomic_radii']
                    })
            string = 'Interatomic hard-sphere distances (pm):\n'
            for item in interatomicHSD:
                string += '    {}: {}\n'.format(item['pair'], item['HSD'])
            logger.info(string)
        else:
            interatomicHSD = []
            for atomic_species_1, details_1 in GUI.core.Project.default_species_dict['atomic_species'].items():
                for atomic_species_2, details_2 in GUI.core.Project.default_species_dict['atomic_species'].items():
                    interatomicHSD.append({
                        'pair': '{}<=>{}'.format(atomic_species_1, atomic_species_2),
                        'HSD': details_1['atomic_radii'] + details_2['atomic_radii']
                    })
            string = 'Interatomic hard-sphere distances (pm):\n'
            for item in interatomicHSD:
                string += '    {}: {}\n'.format(item['pair'], item['HSD'])
            logger.info(string)

    def btn_gen_sub_trigger(self):
        if self.ui_obj.project_instance is not None:
            SubGraphWizard(ui_obj=self.ui_obj)

    def btn_build_anti_graph_trigger(self):
        if self.ui_obj.project_instance is not None:
            if self.ui_obj.project_instance.num_columns > 0:
                if len(self.ui_obj.project_instance.graph.vertices[0].district) > 0:
                    self.ui_obj.sys_message('Working...')
                    anti_graph = self.ui_obj.project_instance.graph.get_anti_graph()
                    self.ui_obj.gs_anti_graph = AntiGraph(ui_obj=self.ui_obj, scale_factor=2, graph=anti_graph)
                    self.ui_obj.gv_anti_graph.setScene(self.ui_obj.gs_anti_graph)
                    logger.info('Got anti-graph!')
                    self.ui_obj.sys_message('Ready')

    def btn_match_zeta_trigger(self):
        if self.ui_obj.project_instance is not None:
            self.ui_obj.sys_message('Working...')
            self.ui_obj.project_instance.graph.match_zeta_graph()
            self.ui_obj.update_display()
            self.ui_obj.sys_message('Ready.')

    def btn_project_plots_trigger(self):
        if self.ui_obj.project_instance is not None:
            PlotProject(ui_obj=self.ui_obj)

    def btn_model_plots_trigger(self):
        if self.ui_obj.project_instance is not None:
            PlotModels(ui_obj=self.ui_obj, model=self.ui_obj.project_instance.graph.active_model)
        else:
            PlotModels(ui_obj=self.ui_obj, model='default_model')

    def btn_calc_models_trigger(self):
        CalcModels(ui_obj=self.ui_obj)

    def btn_calc_heat_trigger(self):
        if self.ui_obj.project_instance is not None:
            HeatMapWizard(ui_obj=self.ui_obj)

    def btn_delete_heat_trigger(self):
        if self.ui_obj.project_instance is not None:
            items = ['None']
            for map_ in self.ui_obj.project_instance.maps:
                items.append(map_['title'])
            item, ok_pressed = QtWidgets.QInputDialog.getItem(self, "Set", "Heat map", items, 0, False)
            if ok_pressed and item:
                if item == 'None':
                    pass
                else:
                    for i, map_ in enumerate(self.ui_obj.project_instance.maps):
                        if map_['title'] == item:
                            del self.ui_obj.project_instance.maps[i]
                            if self.widgets['heat_maps']['sbtn_current_map']['widget'].value == item:
                                self.ui_obj.gs_heat = HeatMap(ui_obj=self, background=self.ui_obj.no_graphic)
                                self.ui_obj.gv_heat.setScene(self.ui_obj.gs_heat)
                                item = 'None'
                                self.widgets['heat_maps']['sbtn_current_map']['widget'].set_value(item)
                            break
                    else:
                        pass

    def chb_heat_legend_trigger(self, state):
        if not self.widgets['heat_maps']['sbtn_current_map']['widget'].value == 'None':
            self.ui_obj.gs_heat.show_legend = state
            if state:
                self.ui_obj.gs_heat.legend.show()
            else:
                self.ui_obj.gs_heat.legend.hide()

    def sbtn_current_map_trigger(self):
        if self.ui_obj.project_instance is not None:
            items = ['None']
            for map_ in self.ui_obj.project_instance.maps:
                items.append(map_['title'])
            item, ok_pressed = QtWidgets.QInputDialog.getItem(self, "Set", "Heat map", items, 0, False)
            if ok_pressed and item:
                if item == 'None':
                    self.ui_obj.gs_heat = HeatMap(ui_obj=self, background=self.ui_obj.no_graphic)
                    self.ui_obj.gv_heat.setScene(self.ui_obj.gs_heat)
                else:
                    for map_ in self.ui_obj.project_instance.maps:
                        if map_['title'] == item:
                            if map_ is not None:
                                self.ui_obj.gs_heat = HeatMap(
                                    ui_obj=self.ui_obj,
                                    kernel_size=map_['kernel']['size'],
                                    step_size=map_['kernel']['step_size'],
                                    attribute=map_['attribute'],
                                    measure_type=map_['measure']
                                )
                                self.ui_obj.gs_heat.heat_map = utils.normalize_array(map_['heat_mat'])
                                min_ = self.ui_obj.gs_heat.heat_map.min()
                                max_ = self.ui_obj.gs_heat.heat_map.max()
                                heat_map = utils.normalize_array(self.ui_obj.gs_heat.heat_map, 1)
                                utils.im_out_static(heat_map.astype(np.float64), 'Images\Outputs\Buffers\\heat_map.png')
                                self.ui_obj.gs_heat.addPixmap(QtGui.QPixmap('Images\Outputs\Buffers\\heat_map.png'))

                                self.ui_obj.gs_heat.legend = GUI_custom_components.HeatLegend(
                                    ui_obj=self.ui_obj,
                                    min=min_,
                                    max=max_
                                )
                                self.ui_obj.gs_heat.addItem(self.ui_obj.gs_heat.legend)
                                if self.ui_obj.gs_heat.show_legend:
                                    self.ui_obj.gs_heat.legend.show()
                                else:
                                    self.ui_obj.gs_heat.legend.hide()
                                self.ui_obj.gv_heat.setScene(self.ui_obj.gs_heat)
                            break
                    else:
                        self.ui_obj.gs_heat = HeatMap(ui_obj=self, background=self.ui_obj.no_graphic)
                        self.ui_obj.gv_heat.setScene(self.ui_obj.gs_heat)
                        item = 'None'
                self.widgets['heat_maps']['sbtn_current_map']['widget'].set_value(item)

    def btn_configure_trigger(self):
        CustomizeOverlay(ui_obj=self.ui_obj)

    def btn_show_all_trigger(self):
        if self.ui_obj.project_instance is not None:
            self.ui_obj.sys_message('Working...')
            for vertex in self.ui_obj.project_instance.graph.vertices:
                vertex.show_in_overlay = True
            for graphic_item in self.ui_obj.gs_overlay_composition.interactive_overlay_objects:
                graphic_item.set_style()
            for channel_dict in self.widgets['overlay'].values():
                if channel_dict['widget_type'] == 'chb':
                    channel_dict['widget'].set_state_no_trigger(True)
            self.ui_obj.sys_message('Ready.')

    def btn_hide_all_trigger(self):
        if self.ui_obj.project_instance is not None:
            self.ui_obj.sys_message('Working...')
            for vertex in self.ui_obj.project_instance.graph.vertices:
                vertex.show_in_overlay = False
            for graphic_item in self.ui_obj.gs_overlay_composition.interactive_overlay_objects:
                graphic_item.set_style()
            for channel_dict in self.widgets['overlay'].values():
                if channel_dict['widget_type'] == 'chb':
                    channel_dict['widget'].set_state_no_trigger(False)
            self.ui_obj.sys_message('Ready.')

    def sbtn_starting_index_trigger(self):
        if self.ui_obj.project_instance is not None and not self.ui_obj.selected_column == -1:
            self.ui_obj.project_instance.starting_index = self.ui_obj.selected_column
            self.widgets['debug']['sbtn_starting_index']['widget'].set_value(
                '{}'.format(self.ui_obj.project_instance.starting_index)
            )

    def btn_test_trigger(self):
        pass

    def btn_crash_trigger(self):
        raise IndexError

    def btn_test_data_manager_trigger(self):
        GUI.test_module.test_statistical_basics()

    def btn_show_prediction_trigger(self):
        if self.ui_obj.project_instance is not None and not self.ui_obj.selected_column == -1:
            vertex = self.ui_obj.project_instance.graph.vertices[self.ui_obj.selected_column]
            model = data_module.VertexDataManager.load(self.ui_obj.project_instance.graph.active_model)
            attributes = model.attribute_keys
            self.ui_obj.project_instance.graph.calc_vertex_parameters(self.ui_obj.selected_column)
            data_line = {}
            for attribute in attributes:
                data_line[attribute] = eval('vertex.{}'.format(attribute))
            if model.category_key == 'advanced_species':
                species_list = self.ui_obj.project_instance.graph.get_advanced_species_list()
            else:
                species_list = self.ui_obj.project_instance.graph.get_atomic_species_list()
            prediction = model.calc_prediction(data_line, species_list)
            prediction = utils.normalize_dict(prediction, 1)
            alpha_data = {}
            for key in data_line:
                if key == 'alpha_max' or key == 'alpha_min':
                    alpha_data[key] = data_line[key]
            alpha_prediction = model.calc_prediction(alpha_data, species_list)

            string = 'Model prediction for vertex {}:\n'.format(vertex.i)
            string += '    Attributes:\n'
            for key, value in data_line.items():
                string += '        {}: {}\n'.format(key, value)
            string += '    Alpha prediction:\n'
            for key, value in alpha_prediction.items():
                string += '        {}: {}\n'.format(key, value)
            string += '    Full prediction:\n'
            for key, value in prediction.items():
                string += '        {}: {}\n'.format(key, value)
            logger.info(string)

    def chb_raw_image_trigger(self, state):
        if self.ui_obj.project_instance is not None:
            self.ui_obj.sys_message('Working...')
            if state:
                self.ui_obj.gs_overlay_composition.pixmap.show()
            else:
                self.ui_obj.gs_overlay_composition.pixmap.hide()
            self.ui_obj.sys_message('Ready.')

    def chb_black_background_trigger(self, state):
        if state:
            self.ui_obj.gs_overlay_composition.setBackgroundBrush(GUI_settings.brush_black)
        else:
            if GUI_settings.theme == 'dark':
                self.ui_obj.gs_overlay_composition.setBackgroundBrush(GUI_settings.background_brush)
            else:
                self.ui_obj.gs_overlay_composition.setBackgroundBrush(GUI_settings.brush_white)

    def chb_edge_columns_trigger(self, state):
        if self.ui_obj.project_instance is not None and self.ui_obj.project_instance.num_columns > 0:
            self.ui_obj.sys_message('Working...')
            for vertex in self.ui_obj.project_instance.graph.vertices:
                if vertex.is_edge_column:
                    vertex.show_in_overlay = state
            for graphic_item in self.ui_obj.gs_overlay_composition.interactive_overlay_objects:
                graphic_item.set_style()
            self.ui_obj.sys_message('Ready.')

    def chb_precipitate_trigger(self, state):
        if self.ui_obj.project_instance is not None and self.ui_obj.project_instance.num_columns > 0:
            self.ui_obj.sys_message('Working...')
            for vertex in self.ui_obj.project_instance.graph.vertices:
                if vertex.is_in_precipitate:
                    vertex.show_in_overlay = state
            for graphic_item in self.ui_obj.gs_overlay_composition.interactive_overlay_objects:
                graphic_item.set_style()
            self.ui_obj.sys_message('Ready.')

    def chb_matrix_trigger(self, state):
        if self.ui_obj.project_instance is not None and self.ui_obj.project_instance.num_columns > 0:
            self.ui_obj.sys_message('Working...')
            for vertex in self.ui_obj.project_instance.graph.vertices:
                if not vertex.is_in_precipitate:
                    vertex.show_in_overlay = state
            for graphic_item in self.ui_obj.gs_overlay_composition.interactive_overlay_objects:
                graphic_item.set_style()
            self.ui_obj.sys_message('Ready.')

    def chb_legend_trigger(self, state):
        if self.ui_obj.project_instance is not None:
            if state:
                self.ui_obj.gs_overlay_composition.legend.show()
            else:
                self.ui_obj.gs_overlay_composition.legend.hide()

    def chb_scalebar_trigger(self, state):
        if self.ui_obj.project_instance is not None:
            if state:
                self.ui_obj.gs_raw_image.scale_bar.show()
                self.ui_obj.gs_overlay_composition.scale_bar.show()
            else:
                self.ui_obj.gs_raw_image.scale_bar.hide()
                self.ui_obj.gs_overlay_composition.scale_bar.hide()

    def chb_zeta_0_trigger(self, state):
        if self.ui_obj.project_instance is not None and self.ui_obj.project_instance.num_columns > 0:
            self.ui_obj.sys_message('Working...')
            for vertex in self.ui_obj.project_instance.graph.vertices:
                if vertex.zeta == 0:
                    vertex.show_in_overlay = state
            for graphic_item in self.ui_obj.gs_overlay_composition.interactive_overlay_objects:
                graphic_item.set_style()
            self.ui_obj.sys_message('Ready.')

    def chb_zeta_1_trigger(self, state):
        if self.ui_obj.project_instance is not None and self.ui_obj.project_instance.num_columns > 0:
            self.ui_obj.sys_message('Working...')
            for vertex in self.ui_obj.project_instance.graph.vertices:
                if vertex.zeta == 1:
                    vertex.show_in_overlay = state
            for graphic_item in self.ui_obj.gs_overlay_composition.interactive_overlay_objects:
                graphic_item.set_style()
            self.ui_obj.sys_message('Ready.')

    def chb_dynamic_trigger(self):
        if self.ui_obj.project_instance is not None and self.ui_obj.project_instance.num_columns > 0:
            self.ui_obj.sys_message('Working...')
            dict_ = self.ui_obj.project_instance.species_dict
            for species in dict_[self.ui_obj.overlay_settings['display_mode']]:
                if not self.widgets['overlay']['chb_{}'.format(species)]['widget'].isChecked() == self.species_states['chb_{}'.format(species)]:
                    for vertex in self.ui_obj.project_instance.graph.vertices:
                        if eval('vertex.{}'.format(self.ui_obj.overlay_settings['display_mode'])) == species:
                            vertex.show_in_overlay = self.widgets['overlay']['chb_{}'.format(species)]['widget'].isChecked()
                    for graphic_item in self.ui_obj.gs_overlay_composition.interactive_overlay_objects:
                        graphic_item.set_style()
                    self.species_states['chb_{}'.format(species)] = self.widgets['overlay']['chb_{}'.format(species)]['widget'].isChecked()
                    break
            self.ui_obj.sys_message('Ready.')

    def btn_check_integrity_trigger(self):
        if self.ui_obj.project_instance is not None:
            self.ui_obj.project_instance.graph.check_integrity()


class MenuBar:

    def __init__(self, bar_obj, ui_obj):

        self.bar_obj = bar_obj
        self.ui_obj = ui_obj

        # Increase the visibility of separators when 'dark' theme!
        self.bar_obj.setStyleSheet("""
                            QMenu::separator {
                                height: 1px;
                                background: grey;
                                margin-left: 10px;
                                margin-right: 5px;
                            }
                        """)

        # Main headers
        self.file = self.bar_obj.addMenu('File')
        self.view = self.bar_obj.addMenu('View')
        self.debug = self.bar_obj.addMenu('Debug')
        self.plugins_ = self.bar_obj.addMenu('Plugins')
        self.help = self.bar_obj.addMenu('Help')

        # Create actions for menus
        # - file
        new_action = QtWidgets.QAction('New', self.ui_obj)
        open_action = QtWidgets.QAction('Open', self.ui_obj)
        save_action = QtWidgets.QAction('Save', self.ui_obj)
        close_action = QtWidgets.QAction('Close', self.ui_obj)
        exit_action = QtWidgets.QAction('Exit', self.ui_obj)
        # - view
        update_display_action = QtWidgets.QAction('Update display', self.ui_obj)
        # - Debug
        self.advanced_debug_mode_action = QtWidgets.QAction('Advanced debug mode', self.ui_obj)
        self.advanced_debug_mode_action.setCheckable(True)
        self.advanced_debug_mode_action.blockSignals(True)
        self.advanced_debug_mode_action.setChecked(False)
        self.advanced_debug_mode_action.blockSignals(False)
        add_mark_action = QtWidgets.QAction('Add mark to terminal', self.ui_obj)
        reset_flags_action = QtWidgets.QAction('Reset all flags', self.ui_obj)
        set_control_file_action = QtWidgets.QAction('Set control instance', self.ui_obj)
        display_deviations_action = QtWidgets.QAction('Display deviation stats', self.ui_obj)
        run_validation_test_action = QtWidgets.QAction('Run algorithm benchmark', self.ui_obj)
        test_consistency_action = QtWidgets.QAction('Reset levels', self.ui_obj)
        ad_hoc_action = QtWidgets.QAction('Ad Hoc functionality', self.ui_obj)
        # - Plugins
        self.plugin_actions = []
        plugin_paths = []
        for plugin in pathlib.Path('plugins/').glob('*.py'):
            plugin_paths.append(plugin)
            self.plugin_actions.append(QtWidgets.QAction(os.path.splitext(plugin.name)[0], self.ui_obj))
        # - help
        self.toggle_tooltips_action = QtWidgets.QAction('Show tooltips', self.ui_obj)
        self.toggle_tooltips_action.setCheckable(True)
        self.toggle_tooltips_action.setChecked(GUI_settings.tooltips)
        set_theme_action = QtWidgets.QAction('Set theme', self.ui_obj)
        there_is_no_help_action = QtWidgets.QAction('HJALP!', self.ui_obj)

        # Add actions to menus
        # - file
        self.file.addAction(new_action)
        self.file.addAction(open_action)
        self.file.addAction(save_action)
        self.file.addAction(close_action)
        self.file.addAction(exit_action)
        # - View
        self.view.addAction(update_display_action)
        self.view.addSeparator()
        # - Debug
        self.debug.addAction(self.advanced_debug_mode_action)
        self.debug.addAction(add_mark_action)
        self.debug.addAction(reset_flags_action)
        self.debug.addAction(set_control_file_action)
        self.debug.addAction(run_validation_test_action)
        self.debug.addAction(display_deviations_action)
        self.debug.addAction(test_consistency_action)
        self.debug.addAction(ad_hoc_action)
        # - Plugins
        for plugin in self.plugin_actions:
            self.plugins_.addAction(plugin)
        # - Help
        self.help.addAction(self.toggle_tooltips_action)
        self.help.addSeparator()
        self.help.addAction(set_theme_action)
        self.help.addAction(there_is_no_help_action)

        # Events
        # - file
        new_action.triggered.connect(self.ui_obj.menu_new_trigger)
        open_action.triggered.connect(self.ui_obj.menu_open_trigger)
        save_action.triggered.connect(self.ui_obj.menu_save_trigger)
        close_action.triggered.connect(self.ui_obj.menu_close_trigger)
        exit_action.triggered.connect(self.ui_obj.menu_exit_trigger)
        # - view
        update_display_action.triggered.connect(self.ui_obj.menu_update_display)
        # - debug
        self.advanced_debug_mode_action.triggered.connect(self.ui_obj.menu_toggle_debug_mode_trigger)
        add_mark_action.triggered.connect(self.ui_obj.menu_add_mark_trigger)
        reset_flags_action.triggered.connect(self.ui_obj.menu_clear_flags_trigger)
        set_control_file_action.triggered.connect(self.ui_obj.menu_set_control_file_trigger)
        run_validation_test_action.triggered.connect(self.ui_obj.menu_run_benchmark_trigger)
        display_deviations_action.triggered.connect(self.ui_obj.menu_display_deviations_trigger)
        test_consistency_action.triggered.connect(self.ui_obj.menu_test_consistency_trigger)
        ad_hoc_action.triggered.connect(self.ui_obj.menu_ad_hoc_trigger)
        # - Plugins (lol, this is so hacky... There must be a better way! - Raymond Hettinger)
        with open('plugin_modules.py', 'w') as imp:
            for path in plugin_paths:
                imp.writelines('import plugins.{}\n'.format(os.path.splitext(path.name)[0]))
        import plugin_modules
        self.plugin_instances = []
        for path, action in zip(plugin_paths, self.plugin_actions):
            module_name = os.path.splitext(path.name)[0]
            plugin_instance = eval('plugin_modules.plugins.{}.Bridge(self.ui_obj)'.format(module_name))
            self.plugin_instances.append(plugin_instance)
            action.triggered.connect(plugin_instance.trigger)
        # - help
        self.toggle_tooltips_action.triggered.connect(self.ui_obj.menu_toggle_tooltips_trigger)
        set_theme_action.triggered.connect(self.ui_obj.menu_set_theme_trigger)
        there_is_no_help_action.triggered.connect(self.ui_obj.menu_there_is_no_help_trigger)

# -----------
# wizards
# -----------


class ImportWizard(QtWidgets.QDialog):

    def __init__(self, *args, ui_obj=None):
        super().__init__(*args)

        self.ui_obj = ui_obj

        self.setWindowTitle('Import')

        self.btn_next = QtWidgets.QPushButton('Next')
        self.btn_next.clicked.connect(self.btn_next_trigger)
        self.btn_cancel = QtWidgets.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(self.btn_cancel_trigger)

        self.cmb_export = QtWidgets.QComboBox()
        self.cmb_export.addItem('Digital Micrograph 3 (.dm3)')
        self.cmb_export.addItem('TIFF')
        self.cmb_export.addItem('AtoMap')
        self.lbl_export = QtWidgets.QLabel('Select export type: ')

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.cmb_export)
        btn_layout = QtWidgets.QHBoxLayout()
        btn_layout.addStretch()
        btn_layout.addWidget(self.btn_cancel)
        btn_layout.addWidget(self.btn_next)
        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        self.setLayout(layout)
        self.exec_()

    def btn_next_trigger(self):
        self.close()
        if self.cmb_export.currentText() == 'Digital Micrograph 3 (.dm3)':
            self.ui_obj.import_trigger('dm3')
        elif self.cmb_export.currentText() == 'TIFF':
            message = QtWidgets.QMessageBox()
            message.setText('Import format not yet supported!')
            message.exec_()
        elif self.cmb_export.currentText() == 'AtoMap':
            scale, ok_pressed = QtWidgets.QInputDialog.getDouble(self, "Set scale", "Please set the image scale (in pm/pixel) with\nat least 5 significant figures:", 6.0000, 0, 1000, 4)
            if ok_pressed:
                self.ui_obj.import_trigger('AtoMap', scale=scale)
        else:
            message = QtWidgets.QMessageBox()
            message.setText('Import format not supported!')
            message.exec_()

    def btn_cancel_trigger(self):
        self.close()


class ExportWizard(QtWidgets.QDialog):

    def __init__(self, *args, ui_obj=None):
        super().__init__(*args)

        self.ui_obj = ui_obj

        self.setWindowTitle('Export wizard')

        self.btn_next = QtWidgets.QPushButton('Next')
        self.btn_next.clicked.connect(self.btn_next_trigger)
        self.btn_cancel = QtWidgets.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(self.btn_cancel_trigger)

        self.cmb_export = QtWidgets.QComboBox()
        self.cmb_export.addItem('Comma Separated Values (.csv)')
        self.cmb_export.addItem('Scalable Vector Graphics (.svg)')
        self.cmb_export.addItem('Portable Network Graphics (.png)')
        self.lbl_export = QtWidgets.QLabel('Select export type: ')

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.cmb_export)
        btn_layout = QtWidgets.QHBoxLayout()
        btn_layout.addStretch()
        btn_layout.addWidget(self.btn_cancel)
        btn_layout.addWidget(self.btn_next)
        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        self.setLayout(layout)
        self.exec_()

    def btn_next_trigger(self):
        self.close()
        if self.cmb_export.currentText() == 'Comma Separated Values (.csv)':
            ExportCSV(ui_obj=self.ui_obj)
        elif self.cmb_export.currentText() == 'Scalable Vector Graphics (.svg)':
            ExportSVG(ui_obj=self.ui_obj)
        elif self.cmb_export.currentText() == 'Portable Network Graphics (.png)':
            ExportPNG(ui_obj=self.ui_obj)
        else:
            message = QtWidgets.QMessageBox()
            message.setText('Export format not supported!')
            message.exec_()

    def btn_cancel_trigger(self):
        self.close()


class ExportCSV(QtWidgets.QDialog):

    def __init__(self, *args, ui_obj=None):
        super().__init__(*args)

        self.ui_obj = ui_obj

        self.setWindowTitle('Export wizard')

        self.page_index = 0
        self.num_frames = 4

        self.btn_next = QtWidgets.QPushButton('Next')
        self.btn_next.clicked.connect(self.btn_next_trigger)
        self.btn_back = QtWidgets.QPushButton('Back')
        self.btn_back.clicked.connect(self.btn_back_trigger)
        self.btn_cancel = QtWidgets.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(self.btn_cancel_trigger)

        self.btn_layout = QtWidgets.QHBoxLayout()
        self.stack_layout = QtWidgets.QStackedLayout()
        self.top_layout = QtWidgets.QVBoxLayout()

        # Frame 0
        self.lbl_file = QtWidgets.QLabel('Export data from current project or enter a list of projects')
        self.rbtn_current_project = QtWidgets.QRadioButton('Current project')
        self.rbtn_list_of_projects = QtWidgets.QRadioButton('Enter list of projects files')
        self.lst_files = QtWidgets.QListWidget()
        self.btn_add_files = QtWidgets.QPushButton('Add files')

        # Frame 1
        self.list_1 = QtWidgets.QListWidget()
        self.list_2 = QtWidgets.QListWidget()
        self.btn_list_1_up = QtWidgets.QPushButton('Move up')
        self.btn_list_1_down = QtWidgets.QPushButton('Move down')
        self.btn_list_2_up = QtWidgets.QPushButton('Move up')
        self.btn_list_2_down = QtWidgets.QPushButton('Move down')
        self.btn_add = QtWidgets.QPushButton('Add')
        self.btn_remove = QtWidgets.QPushButton('Remove')
        self.lbl_included_data = QtWidgets.QLabel('Included data-attributes:')
        self.lbl_available_data = QtWidgets.QLabel('Available data-attributes:')

        # Frame 2
        self.lbl_filter = QtWidgets.QLabel('Set exclusionn filter: (Columns with checked properties will not be included)')
        self.chb_edge_columns = QtWidgets.QCheckBox('Exclude edge columns')
        self.chb_matrix_columns = QtWidgets.QCheckBox('Exclude matrix columns')
        self.chb_particle_columns = QtWidgets.QCheckBox('Exclude particle columns')
        self.chb_hidden_columns = QtWidgets.QCheckBox('Exclude columns that are set as hidden in the overlay')
        self.chb_flag_1 = QtWidgets.QCheckBox('Exclude columns where flag 1 is set to True')
        self.chb_flag_2 = QtWidgets.QCheckBox('Exclude columns where flag 2 is set to True')
        self.chb_flag_3 = QtWidgets.QCheckBox('Exclude columns where flag 3 is set to True')
        self.chb_flag_4 = QtWidgets.QCheckBox('Exclude columns where flag 4 is set to True')

        # Frame 4
        self.chb_recalculate_graphs = QtWidgets.QCheckBox('Recalculate graph data before compiling data (might be very slow)')

        self.set_layout()
        self.exec_()

    def set_layout(self):
        self.btn_layout.addStretch()
        self.btn_layout.addWidget(self.btn_cancel)
        self.btn_layout.addWidget(self.btn_back)
        self.btn_layout.addWidget(self.btn_next)
        self.btn_layout.addStretch()

        self.set_frame_0_layout()
        self.set_frame_1_layout()
        self.set_frame_2_layout()
        self.set_frame_3_layout()

        self.top_layout.addLayout(self.stack_layout)
        self.top_layout.addLayout(self.btn_layout)

        self.setLayout(self.top_layout)

    def set_frame_0_layout(self):
        h_layout = QtWidgets.QHBoxLayout()
        v_layout = QtWidgets.QVBoxLayout()

        v_layout.addStretch()
        v_layout.addWidget(self.lbl_file)
        v_layout.addWidget(self.rbtn_current_project)
        v_layout.addWidget(self.rbtn_list_of_projects)
        v_layout.addWidget(self.lst_files)
        v_layout.addWidget(self.btn_add_files)
        v_layout.addStretch()

        h_layout.addStretch()
        h_layout.addLayout(v_layout)
        h_layout.addStretch()

        self.rbtn_current_project.setChecked(True)
        self.rbtn_list_of_projects.setChecked(False)
        self.lst_files.setDisabled(True)
        self.btn_add_files.setDisabled(True)

        self.lst_files.setMinimumWidth(500)

        self.rbtn_current_project.toggled.connect(self.rbtn_current_project_trigger)
        self.btn_add_files.clicked.connect(self.btn_add_files_trigger)

        widget = QtWidgets.QWidget()
        widget.setLayout(h_layout)
        self.stack_layout.addWidget(widget)

    def set_frame_1_layout(self):
        self.btn_list_1_up.clicked.connect(self.btn_list_1_up_trigger)
        self.btn_list_1_down.clicked.connect(self.btn_list_1_down_trigger)
        self.btn_list_2_up.clicked.connect(self.btn_list_2_up_trigger)
        self.btn_list_2_down.clicked.connect(self.btn_list_2_down_trigger)
        self.btn_add.clicked.connect(self.btn_add_item_trigger)
        self.btn_remove.clicked.connect(self.btn_remove_item_trigger)

        self.list_2.addItems(
            list(GUI.core.graph_2.Vertex.nominal_attributes) + list(GUI.core.graph_2.Vertex.numerical_attributes)
        )

        h_layout = QtWidgets.QHBoxLayout()
        v_layout_1 = QtWidgets.QVBoxLayout()
        v_layout_2 = QtWidgets.QVBoxLayout()
        v_layout_3 = QtWidgets.QVBoxLayout()
        v_layout_4 = QtWidgets.QVBoxLayout()
        v_layout_5 = QtWidgets.QVBoxLayout()

        v_layout_1.addStretch()
        v_layout_1.addWidget(self.btn_list_1_up)
        v_layout_1.addWidget(self.btn_list_1_down)
        v_layout_1.addStretch()

        v_layout_2.addWidget(self.lbl_included_data)
        v_layout_2.addWidget(self.list_1)

        v_layout_3.addStretch()
        v_layout_3.addWidget(self.btn_add)
        v_layout_3.addWidget(self.btn_remove)
        v_layout_3.addStretch()

        v_layout_4.addWidget(self.lbl_available_data)
        v_layout_4.addWidget(self.list_2)

        v_layout_5.addStretch()
        v_layout_5.addWidget(self.btn_list_2_up)
        v_layout_5.addWidget(self.btn_list_2_down)
        v_layout_5.addStretch()

        h_layout.addLayout(v_layout_1)
        h_layout.addLayout(v_layout_2)
        h_layout.addLayout(v_layout_3)
        h_layout.addLayout(v_layout_4)
        h_layout.addLayout(v_layout_5)

        widget = QtWidgets.QWidget()
        widget.setLayout(h_layout)
        self.stack_layout.addWidget(widget)

    def set_frame_2_layout(self):
        self.chb_edge_columns.setChecked(True)
        self.chb_matrix_columns.setChecked(False)
        self.chb_particle_columns.setChecked(False)
        self.chb_hidden_columns.setChecked(False)
        self.chb_flag_1.setChecked(False)
        self.chb_flag_2.setChecked(False)
        self.chb_flag_3.setChecked(False)
        self.chb_flag_4.setChecked(False)

        h_layout = QtWidgets.QHBoxLayout()
        v_layout = QtWidgets.QVBoxLayout()

        v_layout.addStretch()
        v_layout.addWidget(self.lbl_filter)
        v_layout.addWidget(self.chb_edge_columns)
        v_layout.addWidget(self.chb_matrix_columns)
        v_layout.addWidget(self.chb_particle_columns)
        v_layout.addWidget(self.chb_hidden_columns)
        v_layout.addWidget(self.chb_flag_1)
        v_layout.addWidget(self.chb_flag_2)
        v_layout.addWidget(self.chb_flag_3)
        v_layout.addWidget(self.chb_flag_4)
        v_layout.addStretch()

        h_layout.addStretch()
        h_layout.addLayout(v_layout)
        h_layout.addStretch()

        widget = QtWidgets.QWidget()
        widget.setLayout(h_layout)
        self.stack_layout.addWidget(widget)

    def set_frame_3_layout(self):
        self.chb_recalculate_graphs.setChecked(False)
        v_layout = QtWidgets.QVBoxLayout()
        v_layout.addWidget(self.chb_recalculate_graphs)
        v_layout.addStretch()

        widget = QtWidgets.QWidget()
        widget.setLayout(v_layout)
        self.stack_layout.addWidget(widget)

    def btn_list_2_up_trigger(self):
        if self.list_2.currentItem() is not None:
            text = self.list_2.currentItem().text()
            index = self.list_2.currentRow()
            if not index == 0:
                self.list_2.takeItem(self.list_2.row(self.list_2.currentItem()))
                self.list_2.insertItem(index - 1, text)
                self.list_2.setCurrentRow(index - 1)

    def btn_list_2_down_trigger(self):
        if self.list_2.currentItem() is not None:
            text = self.list_2.currentItem().text()
            index = self.list_2.currentRow()
            if not index == self.list_2.count() - 1:
                self.list_2.takeItem(self.list_2.row(self.list_2.currentItem()))
                self.list_2.insertItem(index + 1, text)
                self.list_2.setCurrentRow(index + 1)

    def btn_list_1_up_trigger(self):
        if self.list_1.currentItem() is not None:
            text = self.list_1.currentItem().text()
            index = self.list_1.currentRow()
            if not index == 0:
                self.list_1.takeItem(self.list_1.row(self.list_1.currentItem()))
                self.list_1.insertItem(index - 1, text)
                self.list_1.setCurrentRow(index - 1)

    def btn_list_1_down_trigger(self):
        if self.list_1.currentItem() is not None:
            text = self.list_1.currentItem().text()
            index = self.list_1.currentRow()
            if not index == self.list_1.count() - 1:
                self.list_1.takeItem(self.list_1.row(self.list_1.currentItem()))
                self.list_1.insertItem(index + 1, text)
                self.list_1.setCurrentRow(index + 1)

    def btn_add_item_trigger(self):
        if self.list_2.currentItem() is not None:
            self.list_1.addItem(self.list_2.currentItem().text())
            self.list_2.takeItem(self.list_2.row(self.list_2.currentItem()))

    def btn_remove_item_trigger(self):
        if self.list_1.currentItem() is not None:
            self.list_2.addItem(self.list_1.currentItem().text())
            self.list_1.takeItem(self.list_1.row(self.list_1.currentItem()))

    def btn_add_files_trigger(self):
        prompt = QtWidgets.QFileDialog()
        prompt.setFileMode(QtWidgets.QFileDialog.ExistingFiles)
        if prompt.exec_():
            filenames = prompt.selectedFiles()
        else:
            filenames = None
        if filenames is not None:
            for file_ in filenames:
                exlusive = True
                for i in range(self.lst_files.count()):
                    if file_ == self.lst_files.item(i).text():
                        exlusive = False
                if exlusive:
                    self.lst_files.addItem(file_)

    def rbtn_current_project_trigger(self, state):
        self.lst_files.setDisabled(state)
        self.btn_add_files.setDisabled(state)

    def btn_next_trigger(self):
        if self.stack_layout.currentIndex() == self.num_frames - 1:
            self.export()
        elif self.stack_layout.currentIndex() == self.num_frames -2:
            self.btn_next.setText('export')
            self.stack_layout.setCurrentIndex(self.stack_layout.currentIndex() + 1)
        else:
            self.stack_layout.setCurrentIndex(self.stack_layout.currentIndex() + 1)

    def btn_back_trigger(self):
        if self.stack_layout.currentIndex() == 0:
            self.close()
            ExportWizard(ui_obj=self.ui_obj)
        elif self.stack_layout.currentIndex() == self.num_frames - 1:
            self.btn_next.setText('Next')
            self.stack_layout.setCurrentIndex(self.stack_layout.currentIndex() - 1)
        else:
            self.stack_layout.setCurrentIndex(self.stack_layout.currentIndex() - 1)

    def btn_cancel_trigger(self):
        self.close()

    def export(self):
        self.close()
        self.ui_obj.sys_message('Working...')

        filename = QtWidgets.QFileDialog.getSaveFileName(self, "Save data as", '', ".csv")
        if filename[0]:
            # Prepare string with filenames
            if self.rbtn_list_of_projects.isChecked():
                files = ''
                for i in range(self.lst_files.count()):
                    if not i == self.lst_files.count() - 1:
                        files += self.lst_files.item(i).text() + '\n'
                    else:
                        files += self.lst_files.item(i).text()
            else:
                files = self.ui_obj.savefile
            # Prepare filter
            filter_ = {
                'exclude_edge_columns': self.chb_edge_columns.isChecked(),
                'exclude_matrix_columns': self.chb_matrix_columns.isChecked(),
                'exclude_particle_columns': self.chb_particle_columns.isChecked(),
                'exclude_hidden_columns': self.chb_hidden_columns.isChecked(),
                'exclude_flag_1_columns': self.chb_flag_1.isChecked(),
                'exclude_flag_2_columns': self.chb_flag_2.isChecked(),
                'exclude_flag_3_columns': self.chb_flag_3.isChecked(),
                'exclude_flag_4_columns': self.chb_flag_4.isChecked()
            }
            # Prepare keys
            export_keys = []
            for i in range(self.list_1.count()):
                export_keys.append(self.list_1.item(i).text())
            # Initate manager
            manager = data_module.VertexDataManager(
                files,
                filter_=filter_,
                save_filename=filename[0],
                recalc=self.chb_recalculate_graphs.isChecked(),
                attr_keys=export_keys)
            manager.collect_data()
            manager.export_csv(filename[0] + filename[1], export_keys)
            GUI.logger.info('Successfully exported data to {}'.format(filename[0]))
        self.ui_obj.sys_message('Ready.')


class ExportSVG(QtWidgets.QDialog):

    def __init__(self, *args, ui_obj=None):
        super().__init__(*args)

        self.ui_obj = ui_obj

        self.setWindowTitle('Export wizard')

        self.page_index = 0
        self.num_frames = 3

        self.btn_next = QtWidgets.QPushButton('Next')
        self.btn_next.clicked.connect(self.btn_next_trigger)
        self.btn_back = QtWidgets.QPushButton('Back')
        self.btn_back.clicked.connect(self.btn_back_trigger)
        self.btn_cancel = QtWidgets.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(self.btn_cancel_trigger)

        self.btn_layout = QtWidgets.QHBoxLayout()
        self.stack_layout = QtWidgets.QStackedLayout()
        self.top_layout = QtWidgets.QVBoxLayout()

        # Frame 0
        self.lbl_column_grouping = QtWidgets.QLabel('Group vertices by: ')
        self.cmb_grouping = QtWidgets.QComboBox()
        self.chb_graph = QtWidgets.QCheckBox('Include graph: ')
        self.cmb_graph = QtWidgets.QComboBox()
        self.chb_raw = QtWidgets.QCheckBox('Include image')
        self.lbl_raw = QtWidgets.QLabel('*(This option will include a rectangle in place of the image.\nImport the PNG image with inkscape and fit it to the rectangle..)')
        self.box_radii = QtWidgets.QDoubleSpinBox()

        # Frame 1
        self.lbl_filter = QtWidgets.QLabel('Set exclusionn filter: (Columns with checked properties will not be included)')
        self.chb_edge_columns = QtWidgets.QCheckBox('Exclude edge columns')
        self.chb_matrix_columns = QtWidgets.QCheckBox('Exclude matrix columns')
        self.chb_particle_columns = QtWidgets.QCheckBox('Exclude particle columns')
        self.chb_hidden_columns = QtWidgets.QCheckBox('Exclude columns that are set as hidden in the overlay')
        self.chb_flag_1 = QtWidgets.QCheckBox('Exclude columns where flag 1 is set to True')
        self.chb_flag_2 = QtWidgets.QCheckBox('Exclude columns where flag 2 is set to True')
        self.chb_flag_3 = QtWidgets.QCheckBox('Exclude columns where flag 3 is set to True')
        self.chb_flag_4 = QtWidgets.QCheckBox('Exclude columns where flag 4 is set to True')

        # Frame 2
        self.chb_recalculate_graphs = QtWidgets.QCheckBox('Recalculate graph data before compiling data (might be very slow)')

        self.set_layout()
        self.exec_()

    def set_layout(self):
        self.btn_layout.addStretch()
        self.btn_layout.addWidget(self.btn_cancel)
        self.btn_layout.addWidget(self.btn_back)
        self.btn_layout.addWidget(self.btn_next)
        self.btn_layout.addStretch()

        self.set_frame_0_layout()
        self.set_frame_1_layout()
        self.set_frame_2_layout()

        self.top_layout.addLayout(self.stack_layout)
        self.top_layout.addLayout(self.btn_layout)

        self.setLayout(self.top_layout)

    def set_frame_0_layout(self):
        self.cmb_grouping.addItems([
            'atomic_species',
            'advanced_species'
        ])
        self.cmb_graph.addItems([
            'atomic_graph',
            'zeta_graph'
        ])
        self.box_radii.setMinimum(0.00)
        self.box_radii.setMaximum(1.00)
        self.box_radii.setValue(self.ui_obj.overlay_settings['overlay_radii'])
        self.box_radii.setSingleStep(0.1)

        v_layout = QtWidgets.QVBoxLayout()

        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.lbl_column_grouping)
        layout.addWidget(self.cmb_grouping)

        v_layout.addLayout(layout)

        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(QtWidgets.QLabel('Radius factor (* atomic radii): '))
        layout.addWidget(self.box_radii)

        v_layout.addLayout(layout)

        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.chb_graph)
        layout.addWidget(self.cmb_graph)

        v_layout.addLayout(layout)

        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.chb_raw)
        layout.addWidget(self.lbl_raw)

        v_layout.addLayout(layout)
        v_layout.addStretch()

        widget = QtWidgets.QWidget()
        widget.setLayout(v_layout)
        self.stack_layout.addWidget(widget)

    def set_frame_1_layout(self):
        self.chb_edge_columns.setChecked(True)
        self.chb_matrix_columns.setChecked(False)
        self.chb_particle_columns.setChecked(False)
        self.chb_hidden_columns.setChecked(False)
        self.chb_flag_1.setChecked(False)
        self.chb_flag_2.setChecked(False)
        self.chb_flag_3.setChecked(False)
        self.chb_flag_4.setChecked(False)

        h_layout = QtWidgets.QHBoxLayout()
        v_layout = QtWidgets.QVBoxLayout()

        v_layout.addStretch()
        v_layout.addWidget(self.lbl_filter)
        v_layout.addWidget(self.chb_edge_columns)
        v_layout.addWidget(self.chb_matrix_columns)
        v_layout.addWidget(self.chb_particle_columns)
        v_layout.addWidget(self.chb_hidden_columns)
        v_layout.addWidget(self.chb_flag_1)
        v_layout.addWidget(self.chb_flag_2)
        v_layout.addWidget(self.chb_flag_3)
        v_layout.addWidget(self.chb_flag_4)
        v_layout.addStretch()

        h_layout.addStretch()
        h_layout.addLayout(v_layout)
        h_layout.addStretch()

        widget = QtWidgets.QWidget()
        widget.setLayout(h_layout)
        self.stack_layout.addWidget(widget)

    def set_frame_2_layout(self):
        self.chb_recalculate_graphs.setChecked(False)
        v_layout = QtWidgets.QVBoxLayout()
        v_layout.addWidget(self.chb_recalculate_graphs)
        v_layout.addStretch()

        widget = QtWidgets.QWidget()
        widget.setLayout(v_layout)
        self.stack_layout.addWidget(widget)

    def btn_next_trigger(self):
        if self.stack_layout.currentIndex() == self.num_frames - 1:
            self.export()
        elif self.stack_layout.currentIndex() == self.num_frames -2:
            self.btn_next.setText('export')
            self.stack_layout.setCurrentIndex(self.stack_layout.currentIndex() + 1)
        else:
            self.stack_layout.setCurrentIndex(self.stack_layout.currentIndex() + 1)

    def btn_back_trigger(self):
        if self.stack_layout.currentIndex() == 0:
            self.close()
            ExportWizard(ui_obj=self.ui_obj)
        elif self.stack_layout.currentIndex() == self.num_frames - 1:
            self.btn_next.setText('Next')
            self.stack_layout.setCurrentIndex(self.stack_layout.currentIndex() - 1)
        else:
            self.stack_layout.setCurrentIndex(self.stack_layout.currentIndex() - 1)

    def btn_cancel_trigger(self):
        self.close()

    def export(self):
        self.close()
        if self.ui_obj.project_instance is not None:
            self.ui_obj.sys_message('Working...')
            # Prepare filter
            filter_ = {
                'exclude_edge_columns': self.chb_edge_columns.isChecked(),
                'exclude_matrix_columns': self.chb_matrix_columns.isChecked(),
                'exclude_particle_columns': self.chb_particle_columns.isChecked(),
                'exclude_hidden_columns': self.chb_hidden_columns.isChecked(),
                'exclude_flag_1_columns': self.chb_flag_1.isChecked(),
                'exclude_flag_2_columns': self.chb_flag_2.isChecked(),
                'exclude_flag_3_columns': self.chb_flag_3.isChecked(),
                'exclude_flag_4_columns': self.chb_flag_4.isChecked()
            }
            filename = QtWidgets.QFileDialog.getSaveFileName(self, "Save scalable vector graphics as", '.svg', "")
            if filename[0]:
                if self.chb_graph.isChecked():
                    graph_type = self.cmb_graph.currentText()
                else:
                    graph_type = None
                self.ui_obj.project_instance.export_svg(
                    filename,
                    filter_,
                    species_type=self.cmb_grouping.currentText(),
                    graph=graph_type,
                    image=self.chb_raw.isChecked(),
                    radii_factor=self.box_radii.value()
                )
            self.ui_obj.sys_message('Ready.')
        else:
            message = QtWidgets.QMessageBox()
            message.setText('Open a project first!')
            message.exec_()


class ExportPNG(QtWidgets.QDialog):

    def __init__(self, *args, ui_obj=None):
        super().__init__(*args)

        self.ui_obj = ui_obj

        self.setWindowTitle('Export wizard')

        self.page_index = 0

        self.btn_next = QtWidgets.QPushButton('Export')
        self.btn_next.clicked.connect(self.btn_next_trigger)
        self.btn_back = QtWidgets.QPushButton('Back')
        self.btn_back.clicked.connect(self.btn_back_trigger)
        self.btn_cancel = QtWidgets.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(self.btn_cancel_trigger)

        self.lbl_choose_export = QtWidgets.QLabel('Select what to export:')
        self.rbtn_raw_with = QtWidgets.QRadioButton('Raw image (with scalebar)')
        self.rbtn_raw_without = QtWidgets.QRadioButton('Raw image (without scalebar)')
        self.rbtn_columns = QtWidgets.QRadioButton('Image with column positions overlayed')
        self.rbtn_overlay = QtWidgets.QRadioButton('Atomic overlay (as it is currently displayed)')
        self.rbtn_atomic_graph = QtWidgets.QRadioButton('Atomic graph')
        self.rbtn_zeta_graph = QtWidgets.QRadioButton('Zeta graph')
        self.rbtn_heat_map = QtWidgets.QRadioButton('Currently displayed heat map (without legend)')
        self.rbtn_heat_map_with = QtWidgets.QRadioButton('Currently displayed heat map (with legend)')

        self.set_layout()
        self.exec_()

    def set_layout(self):
        btn_layout = QtWidgets.QHBoxLayout()
        btn_layout.addStretch()
        btn_layout.addWidget(self.btn_back)
        btn_layout.addWidget(self.btn_cancel)
        btn_layout.addWidget(self.btn_next)
        btn_layout.addStretch()

        v_layout = QtWidgets.QVBoxLayout()
        v_layout.addWidget(self.lbl_choose_export)
        v_layout.addWidget(self.rbtn_raw_with)
        v_layout.addWidget(self.rbtn_raw_without)
        v_layout.addWidget(self.rbtn_columns)
        v_layout.addWidget(self.rbtn_overlay)
        v_layout.addWidget(self.rbtn_atomic_graph)
        v_layout.addWidget(self.rbtn_zeta_graph)
        v_layout.addWidget(self.rbtn_heat_map)
        v_layout.addWidget(self.rbtn_heat_map_with)
        v_layout.addLayout(btn_layout)

        self.setLayout(v_layout)

    def btn_cancel_trigger(self):
        self.close()

    def btn_back_trigger(self):
        self.close()
        ExportWizard(ui_obj=self.ui_obj)

    def btn_next_trigger(self):
        self.export()

    def export(self):
        self.close()
        if self.rbtn_raw_without.isChecked():
            self.ui_obj.gs_raw_image.scale_bar.hide()
            self.ui_obj.export_raw_image_trigger()
        elif self.rbtn_raw_with.isChecked():
            self.ui_obj.gs_raw_image.scale_bar.show()
            self.ui_obj.export_raw_image_trigger()
        elif self.rbtn_columns.isChecked():
            self.ui_obj.export_column_position_image_trigger()
        elif self.rbtn_overlay.isChecked():
            self.ui_obj.export_overlay_image_trigger()
        elif self.rbtn_atomic_graph.isChecked():
            self.ui_obj.export_atomic_graph_trigger()
        elif self.rbtn_zeta_graph.isChecked():
            self.ui_obj.export_zeta_graph_trigger()
        elif self.rbtn_heat_map.isChecked():
            self.ui_obj.gs_heat.legend.hide()
            self.ui_obj.export_heat_trigger()
        elif self.rbtn_heat_map_with.isChecked():
            self.ui_obj.gs_heat.legend.show()
            self.ui_obj.export_heat_trigger()
        else:
            message = QtWidgets.QMessageBox()
            message.setText('Unknown export type!')
            message.exec_()


class CalcModels(QtWidgets.QDialog):

    def __init__(self, *args, ui_obj=None):
        super().__init__(*args)

        self.ui_obj = ui_obj

        self.setWindowTitle('Calculate model parameters wizard')

        self.btn_next = QtWidgets.QPushButton('Next')
        self.btn_next.clicked.connect(self.btn_next_trigger)
        self.btn_back = QtWidgets.QPushButton('Back')
        self.btn_back.clicked.connect(self.btn_back_trigger)
        self.btn_cancel = QtWidgets.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(self.btn_cancel_trigger)

        self.btn_layout = QtWidgets.QHBoxLayout()
        self.stack_layout = QtWidgets.QStackedLayout()
        self.top_layout = QtWidgets.QVBoxLayout()

        self.widget_frame_0 = QtWidgets.QWidget()
        self.widget_frame_1 = QtWidgets.QWidget()
        self.widget_frame_2 = QtWidgets.QWidget()
        self.widget_frame_3 = QtWidgets.QWidget()
        self.widget_frame_4 = QtWidgets.QWidget()

        # Frame 0
        self.lbl_file = QtWidgets.QLabel('Calc model from current project or enter a list of projects')
        self.rbtn_current_project = QtWidgets.QRadioButton('Current project')
        self.rbtn_list_of_projects = QtWidgets.QRadioButton('Enter list of projects files')
        self.lst_files = QtWidgets.QListWidget()
        self.btn_add_files = QtWidgets.QPushButton('Add files')

        # Frame 1
        self.lbl_nominal_data = QtWidgets.QLabel('Select a nominal data type to categorize by: ')
        self.cmb_nominal_data = QtWidgets.QComboBox()

        # Frame 2
        self.list_1 = QtWidgets.QListWidget()
        self.list_2 = QtWidgets.QListWidget()
        self.btn_list_1_up = QtWidgets.QPushButton('Move up')
        self.btn_list_1_up.clicked.connect(self.btn_list_1_up_trigger)
        self.btn_list_1_down = QtWidgets.QPushButton('Move down')
        self.btn_list_1_down.clicked.connect(self.btn_list_1_down_trigger)
        self.btn_list_2_up = QtWidgets.QPushButton('Move up')
        self.btn_list_2_up.clicked.connect(self.btn_list_2_up_trigger)
        self.btn_list_2_down = QtWidgets.QPushButton('Move down')
        self.btn_list_2_down.clicked.connect(self.btn_list_2_down_trigger)
        self.btn_add = QtWidgets.QPushButton('Add')
        self.btn_add.clicked.connect(self.btn_add_item_trigger)
        self.btn_remove = QtWidgets.QPushButton('Remove')
        self.btn_remove.clicked.connect(self.btn_remove_item_trigger)
        self.lbl_included_data = QtWidgets.QLabel('Included attributes:')
        self.lbl_available_data = QtWidgets.QLabel('Available attributes:')

        # Frame 3
        self.lbl_filter = QtWidgets.QLabel('Set exclusionn filter: (Columns with checked properties will not be included)')
        self.chb_edge_columns = QtWidgets.QCheckBox('Exclude edge columns')
        self.chb_matrix_columns = QtWidgets.QCheckBox('Exclude matrix columns')
        self.chb_particle_columns = QtWidgets.QCheckBox('Exclude particle columns')
        self.chb_hidden_columns = QtWidgets.QCheckBox('Exclude columns that are set as hidden in the overlay')
        self.chb_flag_1 = QtWidgets.QCheckBox('Exclude columns where flag 1 is set to True')
        self.chb_flag_2 = QtWidgets.QCheckBox('Exclude columns where flag 2 is set to True')
        self.chb_flag_3 = QtWidgets.QCheckBox('Exclude columns where flag 3 is set to True')
        self.chb_flag_4 = QtWidgets.QCheckBox('Exclude columns where flag 4 is set to True')

        # Frame 4
        self.chb_recalculate_graphs = QtWidgets.QCheckBox('Recalculate graph data before compiling data (might be very slow)')

        self.step = 0
        self.state_list = []

        self.set_layout()
        self.exec_()

    def frame_0_layout(self):
        h_layout = QtWidgets.QHBoxLayout()
        v_layout = QtWidgets.QVBoxLayout()

        v_layout.addStretch()
        v_layout.addWidget(self.lbl_file)
        v_layout.addWidget(self.rbtn_current_project)
        v_layout.addWidget(self.rbtn_list_of_projects)
        v_layout.addWidget(self.lst_files)
        v_layout.addWidget(self.btn_add_files)
        v_layout.addStretch()

        h_layout.addStretch()
        h_layout.addLayout(v_layout)
        h_layout.addStretch()

        self.rbtn_current_project.setChecked(True)
        self.rbtn_list_of_projects.setChecked(False)
        self.lst_files.setDisabled(True)
        self.btn_add_files.setDisabled(True)

        self.lst_files.setMinimumWidth(500)

        self.rbtn_current_project.toggled.connect(self.rbtn_current_project_trigger)
        self.btn_add_files.clicked.connect(self.btn_add_files_trigger)

        self.widget_frame_0.setLayout(h_layout)

    def frame_1_layout(self):

        self.cmb_nominal_data.addItems(['advanced_species', 'atomic_species'])

        h_layout = QtWidgets.QHBoxLayout()
        v_layout = QtWidgets.QVBoxLayout()

        v_layout.addStretch()
        v_layout.addWidget(self.lbl_nominal_data)
        v_layout.addWidget(self.cmb_nominal_data)
        v_layout.addStretch()

        h_layout.addStretch()
        h_layout.addLayout(v_layout)
        h_layout.addStretch()

        self.widget_frame_1.setLayout(h_layout)

    def frame_2_layout(self):
        self.list_1.addItems([
            'alpha_min',
            'alpha_max',
            'theta_angle_mean',
            'normalized_avg_gamma',
            'normalized_peak_gamma'
        ])
        self.list_2.addItems([
            'theta_min',
            'theta_max',
            'normalized_area_gamma',
            'redshift'
        ])

        h_layout = QtWidgets.QHBoxLayout()
        v_layout_1 = QtWidgets.QVBoxLayout()
        v_layout_2 = QtWidgets.QVBoxLayout()
        v_layout_3 = QtWidgets.QVBoxLayout()
        v_layout_4 = QtWidgets.QVBoxLayout()
        v_layout_5 = QtWidgets.QVBoxLayout()

        v_layout_1.addStretch()
        v_layout_1.addWidget(self.btn_list_1_up)
        v_layout_1.addWidget(self.btn_list_1_down)
        v_layout_1.addStretch()

        v_layout_2.addWidget(self.lbl_included_data)
        v_layout_2.addWidget(self.list_1)

        v_layout_3.addStretch()
        v_layout_3.addWidget(self.btn_add)
        v_layout_3.addWidget(self.btn_remove)
        v_layout_3.addStretch()

        v_layout_4.addWidget(self.lbl_available_data)
        v_layout_4.addWidget(self.list_2)

        v_layout_5.addStretch()
        v_layout_5.addWidget(self.btn_list_2_up)
        v_layout_5.addWidget(self.btn_list_2_down)
        v_layout_5.addStretch()

        h_layout.addLayout(v_layout_1)
        h_layout.addLayout(v_layout_2)
        h_layout.addLayout(v_layout_3)
        h_layout.addLayout(v_layout_4)
        h_layout.addLayout(v_layout_5)

        self.widget_frame_2.setLayout(h_layout)

    def frame_3_layout(self):
        self.chb_edge_columns.setChecked(True)
        self.chb_matrix_columns.setChecked(False)
        self.chb_particle_columns.setChecked(False)
        self.chb_hidden_columns.setChecked(False)
        self.chb_flag_1.setChecked(False)
        self.chb_flag_2.setChecked(False)
        self.chb_flag_3.setChecked(False)
        self.chb_flag_4.setChecked(False)

        h_layout = QtWidgets.QHBoxLayout()
        v_layout = QtWidgets.QVBoxLayout()

        v_layout.addStretch()
        v_layout.addWidget(self.lbl_filter)
        v_layout.addWidget(self.chb_edge_columns)
        v_layout.addWidget(self.chb_matrix_columns)
        v_layout.addWidget(self.chb_particle_columns)
        v_layout.addWidget(self.chb_hidden_columns)
        v_layout.addWidget(self.chb_flag_1)
        v_layout.addWidget(self.chb_flag_2)
        v_layout.addWidget(self.chb_flag_3)
        v_layout.addWidget(self.chb_flag_4)
        v_layout.addStretch()

        h_layout.addStretch()
        h_layout.addLayout(v_layout)
        h_layout.addStretch()

        self.widget_frame_3.setLayout(h_layout)

    def frame_4_layout(self):
        self.chb_recalculate_graphs.setChecked(False)
        v_layout = QtWidgets.QVBoxLayout()
        v_layout.addWidget(self.chb_recalculate_graphs)
        v_layout.addStretch()
        self.widget_frame_4.setLayout(v_layout)

    def set_layout(self):
        self.btn_layout.addStretch()
        self.btn_layout.addWidget(self.btn_cancel)
        self.btn_layout.addWidget(self.btn_back)
        self.btn_layout.addWidget(self.btn_next)
        self.btn_layout.addStretch()

        self.frame_0_layout()
        self.stack_layout.addWidget(self.widget_frame_0)

        self.frame_1_layout()
        self.stack_layout.addWidget(self.widget_frame_1)

        self.frame_2_layout()
        self.stack_layout.addWidget(self.widget_frame_2)

        self.frame_3_layout()
        self.stack_layout.addWidget(self.widget_frame_3)

        self.frame_4_layout()
        self.stack_layout.addWidget(self.widget_frame_4)

        self.top_layout.addLayout(self.stack_layout)
        self.top_layout.addLayout(self.btn_layout)

        self.setLayout(self.top_layout)

    def btn_add_files_trigger(self):
        prompt = QtWidgets.QFileDialog()
        prompt.setFileMode(QtWidgets.QFileDialog.ExistingFiles)
        if prompt.exec_():
            filenames = prompt.selectedFiles()
        else:
            filenames = None
        if filenames is not None:
            for file_ in filenames:
                exlusive = True
                for i in range(self.lst_files.count()):
                    if file_ == self.lst_files.item(i).text():
                        exlusive = False
                if exlusive:
                    self.lst_files.addItem(file_)

    def rbtn_current_project_trigger(self, state):
        self.lst_files.setDisabled(state)
        self.btn_add_files.setDisabled(state)

    def btn_list_2_up_trigger(self):
        if self.list_2.currentItem() is not None:
            text = self.list_2.currentItem().text()
            index = self.list_2.currentRow()
            if not index == 0:
                self.list_2.takeItem(self.list_2.row(self.list_2.currentItem()))
                self.list_2.insertItem(index - 1, text)
                self.list_2.setCurrentRow(index - 1)

    def btn_list_2_down_trigger(self):
        if self.list_2.currentItem() is not None:
            text = self.list_2.currentItem().text()
            index = self.list_2.currentRow()
            if not index == self.list_2.count() - 1:
                self.list_2.takeItem(self.list_2.row(self.list_2.currentItem()))
                self.list_2.insertItem(index + 1, text)
                self.list_2.setCurrentRow(index + 1)

    def btn_list_1_up_trigger(self):
        if self.list_1.currentItem() is not None:
            text = self.list_1.currentItem().text()
            index = self.list_1.currentRow()
            if not index == 0:
                self.list_1.takeItem(self.list_1.row(self.list_1.currentItem()))
                self.list_1.insertItem(index - 1, text)
                self.list_1.setCurrentRow(index - 1)

    def btn_list_1_down_trigger(self):
        if self.list_1.currentItem() is not None:
            text = self.list_1.currentItem().text()
            index = self.list_1.currentRow()
            if not index == self.list_1.count() - 1:
                self.list_1.takeItem(self.list_1.row(self.list_1.currentItem()))
                self.list_1.insertItem(index + 1, text)
                self.list_1.setCurrentRow(index + 1)

    def btn_add_item_trigger(self):
        if self.list_2.currentItem() is not None:
            self.list_1.addItem(self.list_2.currentItem().text())
            self.list_2.takeItem(self.list_2.row(self.list_2.currentItem()))

    def btn_remove_item_trigger(self):
        if self.list_1.currentItem() is not None:
            self.list_2.addItem(self.list_1.currentItem().text())
            self.list_1.takeItem(self.list_1.row(self.list_1.currentItem()))

    def chb_custom_categories_trigger(self, state):
        if state:
            self.cmb_categorization.setDisabled(True)
            self.cmb_nominal_data.setDisabled(False)
        else:
            self.cmb_categorization.setDisabled(False)
            self.cmb_nominal_data.setDisabled(True)

    def calc(self):
        self.close()
        self.ui_obj.sys_message('Working...')

        filename = QtWidgets.QFileDialog.getSaveFileName(self, "Save model as", '', "")
        if filename[0]:
            # Prepare string with filenames
            if self.rbtn_list_of_projects.isChecked():
                files = ''
                for i in range(self.lst_files.count()):
                    if not i == self.lst_files.count() - 1:
                        files += self.lst_files.item(i).text() + '\n'
                    else:
                        files += self.lst_files.item(i).text()
            else:
                files = self.ui_obj.savefile
            # Prepare filter
            filter_ = {
                'exclude_edge_columns': self.chb_edge_columns.isChecked(),
                'exclude_matrix_columns': self.chb_matrix_columns.isChecked(),
                'exclude_particle_columns': self.chb_particle_columns.isChecked(),
                'exclude_hidden_columns': self.chb_hidden_columns.isChecked(),
                'exclude_flag_1_columns': self.chb_flag_1.isChecked(),
                'exclude_flag_2_columns': self.chb_flag_2.isChecked(),
                'exclude_flag_3_columns': self.chb_flag_3.isChecked(),
                'exclude_flag_4_columns': self.chb_flag_4.isChecked()
            }
            # Prepare keys
            attr_keys = []
            for i in range(self.list_1.count()):
                attr_keys.append(self.list_1.item(i).text())
            cat_key = self.cmb_nominal_data.currentText()

            # Initiate manager
            manager = data_module.VertexDataManager(
                files,
                attr_keys,
                filter_=filter_,
                save_filename=filename[0],
                recalc=self.chb_recalculate_graphs.isChecked(),
                category_key=cat_key
            )
            manager.process_data()
            manager.save()
            GUI.logger.info('Successfully saved model to {}'.format(filename[0]))
        self.ui_obj.sys_message('Ready.')

    def get_next_frame(self):
        next_frame = 0
        if self.stack_layout.currentIndex() == 0:
            next_frame = 1
        elif self.stack_layout.currentIndex() == 1:
            next_frame = 2
        elif self.stack_layout.currentIndex() == 2:
            next_frame = 3
        elif self.stack_layout.currentIndex() == 3:
            next_frame = 4
            self.btn_next.setText('Calcuate!')
        elif self.stack_layout.currentIndex() == 4:
            next_frame = -2
        else:
            logger.error('Error!')
            self.close()

        return next_frame

    def get_previous_frame(self):
        previous_frame = 0
        if self.stack_layout.currentIndex() == 0:
            previous_frame = -1
        elif self.stack_layout.currentIndex() == 1:
            previous_frame = 0
        elif self.stack_layout.currentIndex() == 2:
            previous_frame = 1
        elif self.stack_layout.currentIndex() == 3:
            previous_frame = 2
        elif self.stack_layout.currentIndex() == 4:
            previous_frame = 3
        else:
            logger.error('Error!')
            self.close()

        return previous_frame

    def btn_next_trigger(self):
        next_frame = self.get_next_frame()
        if next_frame == -2:
            self.calc()
        elif next_frame == -1:
            self.close()
        else:
            self.stack_layout.setCurrentIndex(next_frame)

    def btn_back_trigger(self):
        previous_frame = self.get_previous_frame()
        if previous_frame == -1:
            self.close()
        else:
            self.stack_layout.setCurrentIndex(previous_frame)

    def btn_cancel_trigger(self):
        self.close()


class PlotModels(QtWidgets.QDialog):

    def __init__(self, *args, ui_obj=None, model=None):
        super().__init__(*args)

        self.ui_obj = ui_obj

        self.setWindowTitle('Display model details')
        self.model_name = model

        if model:
            self.model = data_module.VertexDataManager.load(model)
        else:
            self.model = data_module.VertexDataManager.load('default_model')

        self.lbl_viewing_model = QtWidgets.QLabel('')
        self.lbl_categories = QtWidgets.QLabel('')
        self.lbl_attributes = QtWidgets.QLabel('')
        self.lbl_size = QtWidgets.QLabel('')
        self.lbl_attribute_1 = QtWidgets.QLabel('Attribute 1: ')
        self.lbl_attribute_2 = QtWidgets.QLabel('Attribute 2: ')
        self.lbl_pca_categories = QtWidgets.QLabel('Show categorization: ')
        self.lbl_single_attribute = QtWidgets.QLabel('Attribute: ')

        self.btn_quit = QtWidgets.QPushButton('Quit')
        self.btn_quit.clicked.connect(self.btn_quit_trigger)
        self.btn_activate_model = QtWidgets.QPushButton('Set as active')
        self.btn_activate_model.clicked.connect(self.btn_activate_model_trigger)
        self.btn_activate_model.setDisabled(True)
        self.btn_load_model = QtWidgets.QPushButton('Load model')
        self.btn_load_model.clicked.connect(self.btn_load_model_trigger)
        self.btn_dual_plot = QtWidgets.QPushButton('Plot')
        self.btn_dual_plot.clicked.connect(self.btn_dual_plot_trigger)
        self.btn_plot_all = QtWidgets.QPushButton('Plot fitted normal distributions')
        self.btn_plot_all.clicked.connect(self.btn_plot_all_trigger)
        self.btn_plot_z_scores = QtWidgets.QPushButton('Plot z-scores')
        self.btn_plot_z_scores.clicked.connect(self.btn_plot_z_scores_trigger)
        self.btn_plot_pca = QtWidgets.QPushButton('Plot 2 first PC')
        self.btn_plot_pca.clicked.connect(self.btn_plot_pca_trigger)
        self.btn_plot_all_pca = QtWidgets.QPushButton('Plot all PC distributions')
        self.btn_plot_all_pca.clicked.connect(self.btn_plot_all_pca_trigger)
        self.btn_plot_single = QtWidgets.QPushButton('Plot')
        self.btn_plot_single.clicked.connect(self.btn_plot_single_trigger)
        self.btn_print_details = QtWidgets.QPushButton('Print details')
        self.btn_print_details.clicked.connect(self.btn_print_details_trigger)
        self.btn_custom = QtWidgets.QPushButton('Custom')
        self.btn_custom.clicked.connect(self.btn_custom_trigger)

        self.cmb_attribute_1 = QtWidgets.QComboBox()
        self.cmb_attribute_2 = QtWidgets.QComboBox()
        self.cmb_single_attribute = QtWidgets.QComboBox()
        self.cmb_pca_setting = QtWidgets.QComboBox()

        self.cmb_pca_setting.addItem('Show categories')
        self.cmb_pca_setting.addItem('No categories')

        self.set_layout()
        self.exec_()

    def set_layout(self):
        # Fill info into widgets:
        self.lbl_viewing_model.setText('Model: {}'.format(self.model.save_filename))
        if self.ui_obj.project_instance is not None:
            if self.model.save_filename == self.ui_obj.project_instance.graph.active_model:
                self.btn_activate_model.setDisabled(True)
            else:
                self.btn_activate_model.setDisabled(False)
        else:
            self.btn_activate_model.setDisabled(True)
        number_of_files = 0
        for file in self.model.files.splitlines(keepends=False):
            number_of_files += 1
        gen_string = ''
        gen_string += 'Total size, n = {}\n'.format(self.model.uncategorized_normal_dist.n)
        gen_string += 'Data gathered from {} images\n'.format(number_of_files)
        gen_string += 'Filter:\n'
        for key, value in self.model.filter_.items():
            gen_string += '    {}: {}\n'.format(key, value)
        self.lbl_size.setText(gen_string)
        attr_string = ''
        for attribute in self.model.attribute_keys:
            self.cmb_attribute_1.addItem(attribute)
            self.cmb_attribute_2.addItem(attribute)
            self.cmb_single_attribute.addItem(attribute)
            attr_string += '{}\n'.format(attribute)
        self.lbl_attributes.setText('{}'.format(attr_string))
        cat_string = ''
        for c, category_key in enumerate(self.model.category_list):
            cat_string += '{}\t\tn = {}\n'.format(category_key, self.model.composite_model[c].n)
        self.lbl_categories.setText('{}'.format(cat_string))

        # Create group boxes:
        grp_set_model = QtWidgets.QGroupBox('Select model')
        grp_dual_plot = QtWidgets.QGroupBox('Dual plot')
        grp_plot_all = QtWidgets.QGroupBox('Plot all model attributes densities')
        grp_plot_pca = QtWidgets.QGroupBox('Plot 2 first principle components')
        grp_categorization = QtWidgets.QGroupBox('Data categories ({}):'.format(self.model.num_data_categories))
        grp_attributes = QtWidgets.QGroupBox('Data attributes ({}):'.format(self.model.k))
        grp_general = QtWidgets.QGroupBox('General:')
        grp_single_plot = QtWidgets.QGroupBox('Plot attribute model distributions')

        # Create group layouts:
        grp_set_model_layout = QtWidgets.QVBoxLayout()
        grp_set_model_layout.addWidget(self.lbl_viewing_model)
        layout = QtWidgets.QHBoxLayout()
        layout.addStretch()
        layout.addWidget(self.btn_activate_model)
        layout.addStretch()
        grp_set_model_layout.addLayout(layout)
        layout = QtWidgets.QHBoxLayout()
        layout.addStretch()
        layout.addWidget(self.btn_load_model)
        layout.addStretch()
        grp_set_model_layout.addLayout(layout)
        grp_set_model.setLayout(grp_set_model_layout)

        grp_dual_plot_layout = QtWidgets.QVBoxLayout()
        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.lbl_attribute_1)
        layout.addWidget(self.cmb_attribute_1)
        grp_dual_plot_layout.addLayout(layout)
        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.lbl_attribute_2)
        layout.addWidget(self.cmb_attribute_2)
        grp_dual_plot_layout.addLayout(layout)
        layout = QtWidgets.QHBoxLayout()
        layout.addStretch()
        layout.addWidget(self.btn_dual_plot)
        layout.addStretch()
        grp_dual_plot_layout.addLayout(layout)
        grp_dual_plot.setLayout(grp_dual_plot_layout)

        grp_plot_all_layout = QtWidgets.QVBoxLayout()
        layout = QtWidgets.QHBoxLayout()
        layout.addStretch()
        layout.addWidget(self.btn_plot_all)
        layout.addStretch()
        grp_plot_all_layout.addLayout(layout)
        layout = QtWidgets.QHBoxLayout()
        layout.addStretch()
        layout.addWidget(self.btn_plot_z_scores)
        layout.addStretch()
        grp_plot_all_layout.addLayout(layout)
        grp_plot_all.setLayout(grp_plot_all_layout)

        grp_plot_pca_layout = QtWidgets.QVBoxLayout()
        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.lbl_pca_categories)
        layout.addWidget(self.cmb_pca_setting)
        grp_plot_pca_layout.addLayout(layout)
        layout = QtWidgets.QHBoxLayout()
        layout.addStretch()
        layout.addWidget(self.btn_plot_pca)
        layout.addWidget(self.btn_plot_all_pca)
        layout.addStretch()
        grp_plot_pca_layout.addLayout(layout)
        grp_plot_pca.setLayout(grp_plot_pca_layout)

        grp_categorization_layout = QtWidgets.QVBoxLayout()
        grp_categorization_layout.addWidget(self.lbl_categories)
        grp_categorization.setLayout(grp_categorization_layout)

        grp_attributes_layout = QtWidgets.QVBoxLayout()
        grp_attributes_layout.addWidget(self.lbl_attributes)
        grp_attributes.setLayout(grp_attributes_layout)

        grp_single_plot_layout = QtWidgets.QVBoxLayout()
        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.lbl_single_attribute)
        layout.addWidget(self.cmb_single_attribute)
        grp_single_plot_layout.addLayout(layout)
        layout = QtWidgets.QHBoxLayout()
        layout.addStretch()
        layout.addWidget(self.btn_plot_single)
        layout.addWidget(self.btn_custom)
        layout.addStretch()
        grp_single_plot_layout.addLayout(layout)
        grp_single_plot.setLayout(grp_single_plot_layout)

        grp_general_layout = QtWidgets.QVBoxLayout()
        grp_general_layout.addWidget(self.lbl_size)
        grp_general_layout.addWidget(self.btn_print_details)
        grp_general.setLayout(grp_general_layout)

        # Set top layout:
        grid = QtWidgets.QGridLayout()
        grid.addWidget(grp_set_model, 0, 0)
        grid.addWidget(grp_general, 1, 0)
        grid.addWidget(grp_categorization, 2, 0)
        grid.addWidget(grp_attributes, 3, 0)
        grid.addWidget(grp_single_plot, 0, 1)
        grid.addWidget(grp_dual_plot, 1, 1)
        grid.addWidget(grp_plot_all, 2, 1)
        grid.addWidget(grp_plot_pca, 3, 1)

        top_layout = QtWidgets.QVBoxLayout()
        top_layout.addLayout(grid)
        layout = QtWidgets.QHBoxLayout()
        layout.addStretch()
        layout.addWidget(self.btn_quit)
        layout.addStretch()
        top_layout.addLayout(layout)

        self.setLayout(top_layout)

    def btn_activate_model_trigger(self):
        self.ui_obj.project_instance.graph.active_model = self.model.save_filename

    def btn_load_model_trigger(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, "Load model", '', "")
        if filename[0]:
            self.close()
            PlotModels(ui_obj=self.ui_obj, model=filename[0])

    def btn_dual_plot_trigger(self):
        self.model.dual_plot(self.cmb_attribute_1.currentText(), self.cmb_attribute_2.currentText())

    def btn_plot_all_trigger(self):
        self.close()
        self.model.plot_all()
        PlotModels(ui_obj=self.ui_obj, model=self.model_name)

    def btn_plot_z_scores_trigger(self):
        self.model.z_plot()

    def btn_plot_pca_trigger(self):
        if self.cmb_pca_setting.currentText() == 'Show categories':
            self.model.plot_pca(show_category=True)
        else:
            self.model.plot_pca(show_category=False)

    def btn_plot_all_pca_trigger(self):
        self.model.plot_all_pc()

    def btn_plot_single_trigger(self):
        self.model.single_plot(self.cmb_single_attribute.currentText())

    def btn_custom_trigger(self):
        self.model.thesis_plot()

    def btn_print_details_trigger(self):
        logger.info(self.model.report())

    def btn_quit_trigger(self):
        self.close()


class CustomizeOverlay(QtWidgets.QDialog):

    def __init__(self, *args, ui_obj=None):
        super().__init__(*args)

        self.ui_obj = ui_obj

        self.setWindowTitle('Customize overlay')

        # Dialog:
        self.btn_apply = QtWidgets.QPushButton('Apply')
        self.btn_apply.clicked.connect(self.btn_apply_trigger)
        self.btn_cancel = QtWidgets.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(self.btn_cancel_changes_trigger)

        # Overlay type:
        self.rdb_atomic = QtWidgets.QRadioButton('Show atomic species overlay')
        self.rdb_atomic.toggled.connect(self.btn_overlay_scheme_trigger)
        self.rdb_advanced = QtWidgets.QRadioButton('Show advanced species overlay')
        self.rdb_advanced.toggled.connect(self.btn_overlay_scheme_trigger)
        if self.ui_obj.overlay_settings['display_mode'] == 'atomic_species':
            self.rdb_atomic.setChecked(True)
        else:
            self.rdb_advanced.setChecked(True)
        self.lbl_atom_radii = QtWidgets.QLabel('Overlay radius (* atomic_radii): ')
        self.box_radius = QtWidgets.QDoubleSpinBox()
        self.box_radius.setMinimum(0.00)
        self.box_radius.setMaximum(1.00)
        self.box_radius.setSingleStep(0.1)
        self.box_radius.setValue(self.ui_obj.overlay_settings['overlay_radii'])

        # Legend
        self.cmb_backgroundtype = QtWidgets.QComboBox()
        self.cmb_backgroundtype.addItem('Transparent')
        self.cmb_backgroundtype.addItem('Black')
        self.cmb_backgroundtype.addItem('White')
        self.cmb_backgroundtype.addItem('Grey')
        if self.ui_obj.overlay_settings['legend_background'] == 'Transparent':
            index = 0
        elif self.ui_obj.overlay_settings['legend_background'] == 'Black':
            index = 1
        elif self.ui_obj.overlay_settings['legend_background'] == 'White':
            index = 2
        elif self.ui_obj.overlay_settings['legend_background'] == 'Grey':
            index = 3
        else:
            index = 0
        self.cmb_backgroundtype.setCurrentIndex(index)

        # Scalebar
        self.box_length = QtWidgets.QSpinBox()
        self.box_length.setMinimum(0)
        self.box_length.setMaximum(1000)
        self.box_length.setValue(self.ui_obj.overlay_settings['scalebar_length'])
        self.cmb_unit = QtWidgets.QComboBox()
        self.cmb_unit.addItem('pm')
        self.cmb_unit.addItem('nm')
        self.cmb_unit.addItem('$\mu$m')
        self.cmb_unit.addItem('Å')
        self.cmb_unit.setCurrentText(self.ui_obj.overlay_settings['scalebar_unit'])
        self.color_scale = GUI_custom_components.RgbaSelector(
            r=self.ui_obj.overlay_settings['scalebar_color'][0],
            g=self.ui_obj.overlay_settings['scalebar_color'][1],
            b=self.ui_obj.overlay_settings['scalebar_color'][2],
            a=255
        )

        self.set_layout()
        self.exec_()

    def set_layout(self):
        layout = QtWidgets.QVBoxLayout()

        layout_0 = QtWidgets.QHBoxLayout()
        layout_0.addWidget(self.lbl_atom_radii)
        layout_0.addWidget(self.box_radius)

        layout_1 = QtWidgets.QHBoxLayout()
        layout_1.addWidget(QtWidgets.QLabel('Legend background: '))
        layout_1.addWidget(self.cmb_backgroundtype)

        layout_2 = QtWidgets.QHBoxLayout()
        layout_2.addWidget(QtWidgets.QLabel('Scale length: '))
        layout_2.addWidget(self.box_length)

        layout_3 = QtWidgets.QHBoxLayout()
        layout_3.addWidget(QtWidgets.QLabel('Scale unit: '))
        layout_3.addWidget(self.cmb_unit)

        layout_4 = QtWidgets.QHBoxLayout()
        layout_4.addWidget(QtWidgets.QLabel('Scalebar color: '))
        layout_4.addWidget(self.color_scale)

        layout_5 = QtWidgets.QHBoxLayout()
        layout_5.addWidget(self.btn_cancel)
        layout_5.addWidget(self.btn_apply)

        layout.addWidget(self.rdb_advanced)
        layout.addWidget(self.rdb_atomic)
        layout.addLayout(layout_0)
        layout.addLayout(layout_1)
        layout.addLayout(layout_2)
        layout.addLayout(layout_3)
        layout.addLayout(layout_4)
        layout.addLayout(layout_5)

        self.setLayout(layout)

    def btn_overlay_scheme_trigger(self):
        if self.rdb_advanced.isChecked():
            self.rdb_atomic.setChecked(False)
        elif self.rdb_atomic.isChecked():
            self.rdb_advanced.setChecked(False)

    def btn_apply_trigger(self):
        self.close()
        if self.rdb_advanced.isChecked():
            if not self.ui_obj.overlay_settings['display_mode'] == 'advanced_species':
                self.ui_obj.overlay_settings['display_mode'] = 'advanced_species'
                self.ui_obj.control_window.btn_show_all_trigger()
        else:
            if not self.ui_obj.overlay_settings['display_mode'] == 'atomic_species':
                self.ui_obj.overlay_settings['display_mode'] = 'atomic_species'
                self.ui_obj.control_window.btn_show_all_trigger()
        self.ui_obj.overlay_settings['scalebar_length'] = self.box_length.value()
        self.ui_obj.overlay_settings['scalebar_unit'] = self.cmb_unit.currentText()
        self.ui_obj.overlay_settings['scalebar_color'] = self.color_scale.get_triple()
        self.ui_obj.overlay_settings['legend_background'] = self.cmb_backgroundtype.currentText()
        self.ui_obj.overlay_settings['overlay_radii'] = self.box_radius.value()
        self.ui_obj.update_overlay()
        self.ui_obj.redraw_control_window()
        self.ui_obj.control_window.draw_histogram()

    def btn_cancel_changes_trigger(self):
        self.close()


class EditSpeciesDict(QtWidgets.QDialog):

    def __init__(self, *args, ui_obj=None, dict_=None):
        super().__init__(*args)

        self.ui_obj = ui_obj
        if self.ui_obj.project_instance is None:
            if dict_ is None:
                self.dict = copy.deepcopy(GUI.core.Project.default_species_dict)
            else:
                self.dict = copy.deepcopy(dict_)
        else:
            if dict_ is None:
                self.dict = copy.deepcopy(self.ui_obj.project_instance.species_dict)
            else:
                self.dict = copy.deepcopy(dict_)

        self.setWindowTitle('Customize the species dictionary for this image')

        # General widgets
        self.btn_what_is_this = QtWidgets.QPushButton('What is this?')
        self.btn_what_is_this.clicked.connect(self.btn_what_is_this_trigger)
        self.btn_apply = QtWidgets.QPushButton('Apply')
        self.btn_apply.clicked.connect(self.btn_apply_trigger)
        self.btn_set_default = QtWidgets.QPushButton('Restore default')
        self.btn_set_default.clicked.connect(self.btn_set_default_trigger)
        self.btn_cancel = QtWidgets.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(self.btn_cancel_trigger)

        # Left widgets
        self.left_btn_move_up = QtWidgets.QPushButton('Move up')
        self.left_btn_move_up.clicked.connect(self.left_btn_move_up_trigger)
        self.left_btn_move_down = QtWidgets.QPushButton('Move down')
        self.left_btn_move_down.clicked.connect(self.left_btn_move_down_trigger)

        self.left_lst = QtWidgets.QListWidget()
        for advanced_species in self.dict['advanced_species']:
            self.left_lst.addItem(advanced_species)
        self.left_lst.itemClicked.connect(self.left_update_details)
        self.left_lbl_details = QtWidgets.QLabel()
        self.left_update_details()

        self.left_btn_new = QtWidgets.QPushButton('New')
        self.left_btn_new.clicked.connect(self.left_btn_new_trigger)
        self.left_btn_edit = QtWidgets.QPushButton('Edit')
        self.left_btn_edit.clicked.connect(self.left_btn_edit_trigger)
        self.left_btn_remove = QtWidgets.QPushButton('Remove')
        self.left_btn_remove.clicked.connect(self.left_btn_remove_trigger)

        self.left_group = QtWidgets.QGroupBox('Advanced species')

        # Right widgets
        self.right_btn_move_up = QtWidgets.QPushButton('Move up')
        self.right_btn_move_up.clicked.connect(self.right_btn_move_up_trigger)
        self.right_btn_move_down = QtWidgets.QPushButton('Move down')
        self.right_btn_move_down.clicked.connect(self.right_btn_move_down_trigger)

        self.right_lst = QtWidgets.QListWidget()
        for atomic_species in self.dict['atomic_species']:
            self.right_lst.addItem(atomic_species)
        self.right_lst.itemClicked.connect(self.right_update_details)
        self.right_lbl_details = QtWidgets.QLabel()
        self.right_update_details()

        self.right_btn_new = QtWidgets.QPushButton('New')
        self.right_btn_new.clicked.connect(self.right_btn_new_trigger)
        self.right_btn_edit = QtWidgets.QPushButton('Edit')
        self.right_btn_edit.clicked.connect(self.right_btn_edit_trigger)
        self.right_btn_remove = QtWidgets.QPushButton('Remove')
        self.right_btn_remove.clicked.connect(self.right_btn_remove_trigger)

        self.right_group = QtWidgets.QGroupBox('Atomic species')

        self.set_layout()
        self.exec_()

    def set_layout(self):
        layout = QtWidgets.QVBoxLayout()

        temp_layout = QtWidgets.QHBoxLayout()
        temp_layout.addStretch()
        temp_layout.addWidget(self.btn_what_is_this)
        temp_layout.addStretch()

        layout.addLayout(temp_layout)

        left_layout = QtWidgets.QHBoxLayout()

        temp_layout = QtWidgets.QVBoxLayout()
        temp_layout.addStretch()
        temp_layout.addWidget(self.left_btn_move_up)
        temp_layout.addWidget(self.left_btn_move_down)
        temp_layout.addStretch()

        left_layout.addLayout(temp_layout)
        left_layout.addWidget(self.left_lst)

        btn_layout = QtWidgets.QHBoxLayout()
        btn_layout.addWidget(self.left_btn_new)
        btn_layout.addWidget(self.left_btn_edit)
        btn_layout.addWidget(self.left_btn_remove)

        temp_layout = QtWidgets.QVBoxLayout()
        temp_layout.addLayout(btn_layout)
        temp_layout.addStretch()
        temp_layout.addWidget(self.left_lbl_details)

        left_layout.addLayout(temp_layout)
        self.left_group.setLayout(left_layout)

        right_layout = QtWidgets.QHBoxLayout()

        temp_layout = QtWidgets.QVBoxLayout()
        temp_layout.addStretch()
        temp_layout.addWidget(self.right_btn_move_up)
        temp_layout.addWidget(self.right_btn_move_down)
        temp_layout.addStretch()

        right_layout.addLayout(temp_layout)
        right_layout.addWidget(self.right_lst)

        btn_layout = QtWidgets.QHBoxLayout()
        btn_layout.addWidget(self.right_btn_new)
        btn_layout.addWidget(self.right_btn_edit)
        btn_layout.addWidget(self.right_btn_remove)

        temp_layout = QtWidgets.QVBoxLayout()
        temp_layout.addLayout(btn_layout)
        temp_layout.addStretch()
        temp_layout.addWidget(self.right_lbl_details)

        right_layout.addLayout(temp_layout)
        self.right_group.setLayout(right_layout)

        temp_layout = QtWidgets.QHBoxLayout()
        temp_layout.addWidget(self.left_group)
        temp_layout.addWidget(self.right_group)

        layout.addLayout(temp_layout)

        dialog_layout = QtWidgets.QHBoxLayout()
        dialog_layout.addStretch()
        dialog_layout.addWidget(self.btn_apply)
        dialog_layout.addWidget(self.btn_set_default)
        dialog_layout.addWidget(self.btn_cancel)
        dialog_layout.addStretch()

        layout.addLayout(dialog_layout)

        self.setLayout(layout)

    def update_(self):
        self.left_lst.clear()
        self.right_lst.clear()
        for advanced_species in self.dict['advanced_species']:
            self.left_lst.addItem(advanced_species)
        for atomic_species in self.dict['atomic_species']:
            self.right_lst.addItem(atomic_species)
        self.right_update_details()
        self.left_update_details()

    def right_update_details(self):
        headers = ['Atomic radii (pm)', 'Color (R, G, B)']
        fields = ['atomic_radii', 'color']
        if self.right_lst.currentItem() is None:
            display_string = 'Atomic species: \n'
            for header, field in zip(headers, fields):
                display_string += '{}: \n'.format(header)
            self.right_lbl_details.setText(display_string)
        else:
            atomic_species = self.right_lst.currentItem().text()
            display_string = 'Atomic species: {}\n'.format(atomic_species)
            for header, field in zip(headers, fields):
                display_string += '{}: {}\n'.format(header, self.dict['atomic_species'][atomic_species][field])
            self.right_lbl_details.setText(display_string)

    def right_btn_move_up_trigger(self):
        if self.right_lst.currentItem() is not None:
            text = self.right_lst.currentItem().text()
            index = self.right_lst.currentRow()
            if not index == 0:
                self.right_lst.takeItem(self.right_lst.row(self.right_lst.currentItem()))
                self.right_lst.insertItem(index - 1, text)
                self.right_lst.setCurrentRow(index - 1)

    def right_btn_move_down_trigger(self):
        if self.right_lst.currentItem() is not None:
            text = self.right_lst.currentItem().text()
            index = self.right_lst.currentRow()
            if not index == self.right_lst.count() - 1:
                self.right_lst.takeItem(self.right_lst.row(self.right_lst.currentItem()))
                self.right_lst.insertItem(index + 1, text)
                self.right_lst.setCurrentRow(index + 1)

    def right_btn_new_trigger(self):
        AtomicSpeciesNew(parent=self)
        self.update_()

    def right_btn_edit_trigger(self):
        if self.right_lst.currentItem() is not None:
            if self.right_lst.currentItem().text() == 'Un':
                message = QtWidgets.QMessageBox()
                message.setText(
                    'Can not edit the Un atomic species, because it is needed for backup classification.')
                message.exec_()
            else:
                AtomicSpeciesEdit(parent=self)
                self.left_update_details()
                self.right_update_details()

    def right_btn_remove_trigger(self):
        if self.right_lst.currentItem() is not None:
            if self.right_lst.currentItem().text() == 'Un':
                message = QtWidgets.QMessageBox()
                message.setText('Can not remove the Un_1 advanced species, because it is needed for backup classification.')
                message.exec_()
            else:
                for value in self.dict['advanced_species'].values():
                    if value['atomic_species'] == self.right_lst.currentItem().text():
                        message = QtWidgets.QMessageBox()
                        message.setText('Cannot remove an atomic species while there are advanced species referencing it')
                        message.exec_()
                        break
                else:
                    del self.dict['atomic_species'][self.right_lst.currentItem().text()]
                    self.right_lst.takeItem(self.right_lst.row(self.right_lst.currentItem()))

    def left_update_details(self):
        headers = ['Symmetry number (n)', 'maps to atomic species', 'Color (R, G, B)', 'Description']
        fields = ['n', 'atomic_species', 'color', 'description']
        if self.left_lst.currentItem() is None:
            display_string = 'Advanced species: \n'
            for header, field in zip(headers, fields):
                display_string += '{}: \n'.format(header)
            self.left_lbl_details.setText(display_string)
        else:
            advanced_species = self.left_lst.currentItem().text()
            display_string = 'Advanced species: {}\n'.format(advanced_species)
            for header, field in zip(headers, fields):
                display_string += '{}: {}\n'.format(header, self.dict['advanced_species'][advanced_species][field])
            self.left_lbl_details.setText(display_string)

    def left_btn_move_down_trigger(self):
        if self.left_lst.currentItem() is not None:
            text = self.left_lst.currentItem().text()
            index = self.left_lst.currentRow()
            if not index == self.left_lst.count() - 1:
                self.left_lst.takeItem(self.left_lst.row(self.left_lst.currentItem()))
                self.left_lst.insertItem(index + 1, text)
                self.left_lst.setCurrentRow(index + 1)

    def left_btn_move_up_trigger(self):
        if self.left_lst.currentItem() is not None:
            text = self.left_lst.currentItem().text()
            index = self.left_lst.currentRow()
            if not index == 0:
                self.left_lst.takeItem(self.left_lst.row(self.left_lst.currentItem()))
                self.left_lst.insertItem(index - 1, text)
                self.left_lst.setCurrentRow(index - 1)

    def left_btn_new_trigger(self):
        AdvancedSpeciesNew(parent=self)
        self.update_()

    def left_btn_edit_trigger(self):
        if self.left_lst.currentItem() is not None:
            if self.left_lst.currentItem().text() == 'Un_1':
                message = QtWidgets.QMessageBox()
                message.setText(
                    'Can not edit the Un_1 advanced species, because it is needed for backup classification.')
                message.exec_()
            else:
                AdvancedSpeciesEdit(parent=self)
                self.left_update_details()
                self.right_update_details()

    def left_btn_remove_trigger(self):
        if self.left_lst.currentItem() is not None:
            if self.left_lst.currentItem().text() == 'Un_1':
                message = QtWidgets.QMessageBox()
                message.setText(
                    'Can not remove the Un_1 advanced species, because it is needed for backup classification.'
                )
                message.exec_()
            else:
                del self.dict['advanced_species'][self.left_lst.currentItem().text()]
                self.left_lst.takeItem(self.left_lst.row(self.left_lst.currentItem()))
                self.left_update_details()
                self.right_update_details()

    def btn_apply_trigger(self):
        reference_list = set()
        for value in self.dict['advanced_species'].values():
            reference_list.add(value['atomic_species'])
        for key in self.dict['atomic_species']:
            if key not in reference_list:
                message = QtWidgets.QMessageBox()
                message.setText(
                    'The Atomic species {}, is not referenced by any advanced species!'.format(key))
                message.exec_()
                break
        else:
            self.close()
            if self.ui_obj.project_instance is not None:
                self.ui_obj.sys_message('Working...')
                self.ui_obj.project_instance.species_dict = self.dict
                self.ui_obj.project_instance.graph.species_dict = self.dict
                self.ui_obj.project_instance.graph.check_species()
                self.ui_obj.update_overlay()
                self.ui_obj.redraw_control_window()
                self.ui_obj.sys_message('Ready.')
            else:
                logger.info('Cannot change the default categories! Open a project first')

    def btn_cancel_trigger(self):
        self.close()

    def btn_set_default_trigger(self):
        self.dict = copy.deepcopy(GUI.core.Project.default_species_dict)
        self.update_()

    @staticmethod
    def btn_what_is_this_trigger():
        string = 'This dialog allows you to customize the \"species dictionary\" of the current project. ' \
                 'One of the goals of AutomAl 6000 is to ' \
                 'determine the \"atomic species\" of individual columns. To achieve this, AutomAl 6000 will first ' \
                 'classify columns by an \"advanced species\", which then can be mapped back to atomic species. This ' \
                 'entails that each advanced species entry in the species dictionary, must map to an atomic species.'
        message = QtWidgets.QMessageBox()
        message.setText(string)
        message.exec_()


class AtomicSpeciesNew(QtWidgets.QDialog):

    def __init__(self, *args, parent=None):
        super().__init__(*args)

        self.parent = parent

        self.btn_ok = QtWidgets.QPushButton('Ok')
        self.btn_ok.clicked.connect(self.btn_ok_trigger)
        self.btn_cancel = QtWidgets.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(self.btn_cancel_trigger)

        self.species_text = QtWidgets.QLineEdit('New')
        self.radii_box = QtWidgets.QDoubleSpinBox()
        self.radii_box.setMinimum(0.00)
        self.radii_box.setMaximum(500.0)
        self.color_edit = GUI_custom_components.RgbaSelector()
        self.color_edit.r_box.setValue(100)
        self.color_edit.g_box.setValue(100)
        self.color_edit.b_box.setValue(100)

        self.set_layout()
        self.exec_()

    def set_layout(self):
        v_layout = QtWidgets.QVBoxLayout()

        grid = QtWidgets.QGridLayout()
        grid.addWidget(QtWidgets.QLabel('Atomic species:'), 0, 0)
        grid.addWidget(QtWidgets.QLabel('Atomic radii (pm):'), 0, 1)
        grid.addWidget(QtWidgets.QLabel('Color (r, g, b):'), 0, 2)

        grid.addWidget(self.species_text, 1, 0)
        grid.addWidget(self.radii_box, 1, 1)
        grid.addWidget(self.color_edit, 1, 2)

        v_layout.addLayout(grid)

        btn_layout = QtWidgets.QHBoxLayout()
        btn_layout.addStretch()
        btn_layout.addWidget(self.btn_ok)
        btn_layout.addWidget(self.btn_cancel)
        btn_layout.addStretch()

        v_layout.addLayout(btn_layout)

        self.setLayout(v_layout)

    def btn_ok_trigger(self):
        if self.species_text.text() in self.parent.dict['atomic_species']:
            string = 'The chosen species name already exist!'
            message = QtWidgets.QMessageBox()
            message.setText(string)
            message.exec_()
        else:
            self.close()
            self.parent.dict['atomic_species'][self.species_text.text()] = {
                'atomic_radii': self.radii_box.value(),
                'color': self.color_edit.get_triple()
            }

    def btn_cancel_trigger(self):
        self.close()


class AtomicSpeciesEdit(QtWidgets.QDialog):

    def __init__(self, *args, parent=None):
        super().__init__(*args)

        self.parent = parent
        self.atomic_species = self.parent.right_lst.currentItem().text()

        self.btn_ok = QtWidgets.QPushButton('Ok')
        self.btn_ok.clicked.connect(self.btn_ok_trigger)
        self.btn_cancel = QtWidgets.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(self.btn_cancel_trigger)

        self.species_text = QtWidgets.QLineEdit(self.atomic_species)
        self.radii_box = QtWidgets.QDoubleSpinBox()
        self.radii_box.setMinimum(0.00)
        self.radii_box.setMaximum(500.0)
        self.radii_box.setValue(self.parent.dict['atomic_species'][self.atomic_species]['atomic_radii'])
        self.color_edit = GUI_custom_components.RgbaSelector()
        self.color_edit.r_box.setValue(self.parent.dict['atomic_species'][self.atomic_species]['color'][0])
        self.color_edit.g_box.setValue(self.parent.dict['atomic_species'][self.atomic_species]['color'][1])
        self.color_edit.b_box.setValue(self.parent.dict['atomic_species'][self.atomic_species]['color'][2])

        self.set_layout()
        self.exec_()

    def set_layout(self):
        v_layout = QtWidgets.QVBoxLayout()

        grid = QtWidgets.QGridLayout()
        grid.addWidget(QtWidgets.QLabel('Atomic species:'), 0, 0)
        grid.addWidget(QtWidgets.QLabel('Atomic radii (pm):'), 0, 1)
        grid.addWidget(QtWidgets.QLabel('Color (r, g, b):'), 0, 2)

        grid.addWidget(self.species_text, 1, 0)
        grid.addWidget(self.radii_box, 1, 1)
        grid.addWidget(self.color_edit, 1, 2)

        v_layout.addLayout(grid)

        btn_layout = QtWidgets.QHBoxLayout()
        btn_layout.addStretch()
        btn_layout.addWidget(self.btn_ok)
        btn_layout.addWidget(self.btn_cancel)
        btn_layout.addStretch()

        v_layout.addLayout(btn_layout)

        self.setLayout(v_layout)

    def btn_ok_trigger(self):
        for atomic_species in self.parent.dict['atomic_species']:
            if atomic_species == self.species_text.text() and not atomic_species == self.atomic_species:
                string = 'The chosen species name already exist!'
                message = QtWidgets.QMessageBox()
                message.setText(string)
                message.exec_()
                break
        else:
            self.close()
            if not self.species_text.text() == self.atomic_species:
                del self.parent.dict['atomic_species'][self.advanced_species]
            self.parent.dict['atomic_species'][self.species_text.text()] = {
                'atomic_radii': self.radii_box.value(),
                'color': self.color_edit.get_triple()
            }

    def btn_cancel_trigger(self):
        self.close()


class AdvancedSpeciesNew(QtWidgets.QDialog):

    def __init__(self, *args, parent=None):
        super().__init__(*args)

        self.parent = parent

        self.btn_ok = QtWidgets.QPushButton('Ok')
        self.btn_ok.clicked.connect(self.btn_ok_trigger)
        self.btn_cancel = QtWidgets.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(self.btn_cancel_trigger)

        self.species_text = QtWidgets.QLineEdit('New')
        self.symmetry_box = QtWidgets.QSpinBox()
        self.symmetry_box.setMaximum(10)
        self.symmetry_box.setMinimum(1)
        self.cmb_atom = QtWidgets.QComboBox()
        for atomic_species in self.parent.dict['atomic_species']:
            self.cmb_atom.addItem(atomic_species)
        self.color_edit = GUI_custom_components.RgbaSelector()
        self.color_edit.r_box.setValue(100)
        self.color_edit.g_box.setValue(100)
        self.color_edit.b_box.setValue(100)
        self.description_text = QtWidgets.QLineEdit('')

        self.set_layout()
        self.exec_()

    def set_layout(self):
        v_layout = QtWidgets.QVBoxLayout()

        grid = QtWidgets.QGridLayout()
        grid.addWidget(QtWidgets.QLabel('Advanced species:'), 0, 0)
        grid.addWidget(QtWidgets.QLabel('Symmetry number (n):'), 0, 1)
        grid.addWidget(QtWidgets.QLabel('Atomic species:'), 0, 2)
        grid.addWidget(QtWidgets.QLabel('Color (r, g, b):'), 0, 3)
        grid.addWidget(QtWidgets.QLabel('Description:'), 0, 4)

        grid.addWidget(self.species_text, 1, 0)
        grid.addWidget(self.symmetry_box, 1, 1)
        grid.addWidget(self.cmb_atom, 1, 2)
        grid.addWidget(self.color_edit, 1, 3)
        grid.addWidget(self.description_text, 1, 4)

        v_layout.addLayout(grid)

        btn_layout = QtWidgets.QHBoxLayout()
        btn_layout.addStretch()
        btn_layout.addWidget(self.btn_ok)
        btn_layout.addWidget(self.btn_cancel)
        btn_layout.addStretch()

        v_layout.addLayout(btn_layout)

        self.setLayout(v_layout)

    def btn_ok_trigger(self):
        if self.species_text.text() in self.parent.dict['advanced_species']:
            string = 'The chosen species name already exist!'
            message = QtWidgets.QMessageBox()
            message.setText(string)
            message.exec_()
        else:
            self.close()
            self.parent.dict['advanced_species'][self.species_text.text()] = {
                'n': self.symmetry_box.value(),
                'atomic_species': self.cmb_atom.currentText(),
                'color': self.color_edit.get_triple(),
                'description': self.description_text.text()
            }

    def btn_cancel_trigger(self):
        self.close()


class AdvancedSpeciesEdit(QtWidgets.QDialog):

    def __init__(self, *args, parent=None):
        super().__init__(*args)

        self.parent = parent
        self.advanced_species = self.parent.left_lst.currentItem().text()

        self.btn_ok = QtWidgets.QPushButton('Ok')
        self.btn_ok.clicked.connect(self.btn_ok_trigger)
        self.btn_cancel = QtWidgets.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(self.btn_cancel_trigger)

        self.species_text = QtWidgets.QLineEdit(self.advanced_species)
        self.symmetry_box = QtWidgets.QSpinBox()
        self.symmetry_box.setMaximum(10)
        self.symmetry_box.setMinimum(1)
        self.symmetry_box.setValue(self.parent.dict['advanced_species'][self.advanced_species]['n'])
        self.cmb_atom = QtWidgets.QComboBox()
        for atomic_species in self.parent.dict['atomic_species']:
            self.cmb_atom.addItem(atomic_species)
        self.color_edit = GUI_custom_components.RgbaSelector()
        self.color_edit.r_box.setValue(self.parent.dict['advanced_species'][self.advanced_species]['color'][0])
        self.color_edit.g_box.setValue(self.parent.dict['advanced_species'][self.advanced_species]['color'][1])
        self.color_edit.b_box.setValue(self.parent.dict['advanced_species'][self.advanced_species]['color'][2])
        self.description_text = QtWidgets.QLineEdit(self.parent.dict['advanced_species'][self.advanced_species]['description'])

        self.set_layout()
        self.exec_()

    def set_layout(self):
        v_layout = QtWidgets.QVBoxLayout()

        grid = QtWidgets.QGridLayout()
        grid.addWidget(QtWidgets.QLabel('Advanced species:'), 0, 0)
        grid.addWidget(QtWidgets.QLabel('Symmetry number (n):'), 0, 1)
        grid.addWidget(QtWidgets.QLabel('Atomic species:'), 0, 2)
        grid.addWidget(QtWidgets.QLabel('Color (r, g, b):'), 0, 3)
        grid.addWidget(QtWidgets.QLabel('Description:'), 0, 4)

        grid.addWidget(self.species_text, 1, 0)
        grid.addWidget(self.symmetry_box, 1, 1)
        grid.addWidget(self.cmb_atom, 1, 2)
        grid.addWidget(self.color_edit, 1, 3)
        grid.addWidget(self.description_text, 1, 4)

        v_layout.addLayout(grid)

        btn_layout = QtWidgets.QHBoxLayout()
        btn_layout.addStretch()
        btn_layout.addWidget(self.btn_ok)
        btn_layout.addWidget(self.btn_cancel)
        btn_layout.addStretch()

        v_layout.addLayout(btn_layout)

        self.setLayout(v_layout)

    def btn_ok_trigger(self):
        for advanced_species in self.parent.dict['advanced_species']:
            if advanced_species == self.species_text.text() and not advanced_species == self.advanced_species:
                string = 'The chosen species name already exist!'
                message = QtWidgets.QMessageBox()
                message.setText(string)
                message.exec_()
                break
        else:
            self.close()
            if self.species_text.text() == self.advanced_species:
                self.parent.dict['advanced_species'][self.species_text.text()] = {
                    'n': self.symmetry_box.value(),
                    'atomic_species': self.cmb_atom.currentText(),
                    'color': self.color_edit.get_triple(),
                    'description': self.description_text.text()
                }
            else:
                del self.parent.dict['advanced_species'][self.advanced_species]
                self.parent.dict['advanced_species'][self.species_text.text()] = {
                    'n': self.symmetry_box.value(),
                    'atomic_species': self.cmb_atom.currentText(),
                    'color': self.color_edit.get_triple(),
                    'description': self.description_text.text()
                }

    def btn_cancel_trigger(self):
        self.close()


class HeatMapWizard(QtWidgets.QDialog):

    def __init__(self, *args, ui_obj=None):

        super().__init__(*args)

        self.ui_obj = ui_obj

        self.setWindowTitle('Make heat-map')

        self.cmb_attributes = QtWidgets.QComboBox()
        self.cmb_measure = QtWidgets.QComboBox()

        self.txt_title = QtWidgets.QLineEdit('Heat map')

        self.num_box_kernel = QtWidgets.QDoubleSpinBox()
        self.num_box_kernel.setMaximum(2.0)
        self.num_box_kernel.setMinimum(0.0)
        self.num_box_kernel.setValue(1.0)

        self.num_box_step = QtWidgets.QSpinBox()
        self.num_box_step.setMaximum(1000)
        self.num_box_step.setMinimum(1)

        self.cmb_kernel = QtWidgets.QComboBox()

        self.chb_overlay = QtWidgets.QCheckBox('Overlay on image')

        self.btn_cancel = QtWidgets.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(self.btn_cancel_trigger)
        self.btn_make = QtWidgets.QPushButton('Generate')
        self.btn_make.clicked.connect(self.btn_make_trigger)

        self.populate()
        self.set_layout()
        self.exec_()

    def populate(self):

        for attribute in GUI.core.graph_2.Vertex.numerical_attributes:
            self.cmb_attributes.addItem(attribute)

        self.cmb_measure.addItem('Variance')
        self.cmb_measure.addItem('Mean')
        self.cmb_measure.addItem('Sum')

        self.cmb_kernel.addItem('square')
        self.cmb_kernel.addItem('city')

    def set_layout(self):

        grid = QtWidgets.QGridLayout()

        grid.addWidget(QtWidgets.QLabel('Attribute:'), 0, 0)
        grid.addWidget(QtWidgets.QLabel('Measure:'), 1, 0)
        grid.addWidget(QtWidgets.QLabel('Kernel Size:'), 2, 0)
        grid.addWidget(QtWidgets.QLabel('Step_size:'), 3, 0)
        grid.addWidget(QtWidgets.QLabel('Kernel type: '), 4, 0)
        grid.addWidget(QtWidgets.QLabel('Title: '), 5, 0)

        grid.addWidget(self.cmb_attributes, 0, 1)
        grid.addWidget(self.cmb_measure, 1, 1)
        grid.addWidget(self.num_box_kernel, 2, 1)
        grid.addWidget(self.num_box_step, 3, 1)
        grid.addWidget(self.cmb_kernel, 4, 1)
        grid.addWidget(self.txt_title, 5, 1)

        grid.addWidget(QtWidgets.QLabel(' * a (pm)'), 2, 2)
        grid.addWidget(QtWidgets.QLabel('(Pixels)'), 3, 2)

        v_layout = QtWidgets.QVBoxLayout()
        v_layout.addLayout(grid)
        v_layout.addWidget(self.chb_overlay)

        btn_layout = QtWidgets.QHBoxLayout()
        btn_layout.addWidget(self.btn_cancel)
        btn_layout.addWidget(self.btn_make)

        v_layout.addLayout(btn_layout)

        self.setLayout(v_layout)

    def btn_make_trigger(self):

        self.close()
        self.ui_obj.sys_message('Working...')
        self.ui_obj.gs_heat = HeatMap(
            ui_obj=self.ui_obj,
            kernel_size=self.num_box_kernel.value(),
            step_size=self.num_box_step.value(),
            attribute=self.cmb_attributes.currentText(),
            measure_type=self.cmb_measure.currentText(),
            kernel_type=self.cmb_kernel.currentText(),
            legend=self.ui_obj.control_window.widgets['heat_maps']['chb_legend']['widget'].isChecked(),
            title=self.txt_title.text()
        )
        self.ui_obj.control_window.widgets['heat_maps']['sbtn_current_map']['widget'].set_value(self.txt_title.text())
        self.ui_obj.gs_heat.make()
        self.ui_obj.gv_heat.setScene(self.ui_obj.gs_heat)

        self.ui_obj.sys_message('Ready.')

    def btn_cancel_trigger(self):
        self.close()


class PlotProject(QtWidgets.QDialog):

    def __init__(self, *args, ui_obj=None):
        super().__init__(*args)

        self.ui_obj = ui_obj

        self.setWindowTitle('Plotting wizard')

        self.btn_next = QtWidgets.QPushButton('Next')
        self.btn_next.clicked.connect(self.btn_next_trigger)
        self.btn_back = QtWidgets.QPushButton('Back')
        self.btn_back.clicked.connect(self.btn_back_trigger)
        self.btn_cancel = QtWidgets.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(self.btn_cancel_trigger)

        self.btn_layout = QtWidgets.QHBoxLayout()
        self.stack_layout = QtWidgets.QStackedLayout()
        self.top_layout = QtWidgets.QVBoxLayout()

        self.widget_frame_0 = QtWidgets.QWidget()
        self.widget_frame_1 = QtWidgets.QWidget()
        self.widget_frame_2 = QtWidgets.QWidget()
        self.widget_frame_3 = QtWidgets.QWidget()

        # Frame 0
        self.lbl_nominal_data = QtWidgets.QLabel('Select a nominal data type to categorize by: ')
        self.cmb_nominal_data = QtWidgets.QComboBox()

        # Frame 1
        self.list_1 = QtWidgets.QListWidget()
        self.list_2 = QtWidgets.QListWidget()
        self.btn_list_1_up = QtWidgets.QPushButton('Move up')
        self.btn_list_1_up.clicked.connect(self.btn_list_1_up_trigger)
        self.btn_list_1_down = QtWidgets.QPushButton('Move down')
        self.btn_list_1_down.clicked.connect(self.btn_list_1_down_trigger)
        self.btn_list_2_up = QtWidgets.QPushButton('Move up')
        self.btn_list_2_up.clicked.connect(self.btn_list_2_up_trigger)
        self.btn_list_2_down = QtWidgets.QPushButton('Move down')
        self.btn_list_2_down.clicked.connect(self.btn_list_2_down_trigger)
        self.btn_add = QtWidgets.QPushButton('Add')
        self.btn_add.clicked.connect(self.btn_add_item_trigger)
        self.btn_remove = QtWidgets.QPushButton('Remove')
        self.btn_remove.clicked.connect(self.btn_remove_item_trigger)
        self.lbl_included_data = QtWidgets.QLabel('Included attributes:')
        self.lbl_available_data = QtWidgets.QLabel('Available attributes:')

        # Frame 2
        self.lbl_filter = QtWidgets.QLabel('Set exclusionn filter: (Columns with checked properties will not be included)')
        self.chb_edge_columns = QtWidgets.QCheckBox('Exclude edge columns')
        self.chb_matrix_columns = QtWidgets.QCheckBox('Exclude matrix columns')
        self.chb_particle_columns = QtWidgets.QCheckBox('Exclude particle columns')
        self.chb_hidden_columns = QtWidgets.QCheckBox('Exclude columns that are set as hidden in the overlay')
        self.chb_flag_1 = QtWidgets.QCheckBox('Exclude columns where flag 1 is set to True')
        self.chb_flag_2 = QtWidgets.QCheckBox('Exclude columns where flag 2 is set to True')
        self.chb_flag_3 = QtWidgets.QCheckBox('Exclude columns where flag 3 is set to True')
        self.chb_flag_4 = QtWidgets.QCheckBox('Exclude columns where flag 4 is set to True')

        # Frame 3
        self.chb_recalculate_graphs = QtWidgets.QCheckBox('Recalculate graph data before compiling data (might be very slow)')

        self.set_layout()
        self.exec_()

    def frame_0_layout(self):

        self.cmb_nominal_data.addItems(list(GUI.core.graph_2.Vertex.nominal_attributes))

        h_layout = QtWidgets.QHBoxLayout()
        v_layout = QtWidgets.QVBoxLayout()

        v_layout.addStretch()
        v_layout.addWidget(self.lbl_nominal_data)
        v_layout.addWidget(self.cmb_nominal_data)
        v_layout.addStretch()

        h_layout.addStretch()
        h_layout.addLayout(v_layout)
        h_layout.addStretch()

        self.widget_frame_0.setLayout(h_layout)

    def frame_1_layout(self):
        self.list_2.addItems(list(GUI.core.graph_2.Vertex.numerical_attributes))

        h_layout = QtWidgets.QHBoxLayout()
        v_layout_1 = QtWidgets.QVBoxLayout()
        v_layout_2 = QtWidgets.QVBoxLayout()
        v_layout_3 = QtWidgets.QVBoxLayout()
        v_layout_4 = QtWidgets.QVBoxLayout()
        v_layout_5 = QtWidgets.QVBoxLayout()

        v_layout_1.addStretch()
        v_layout_1.addWidget(self.btn_list_1_up)
        v_layout_1.addWidget(self.btn_list_1_down)
        v_layout_1.addStretch()

        v_layout_2.addWidget(self.lbl_included_data)
        v_layout_2.addWidget(self.list_1)

        v_layout_3.addStretch()
        v_layout_3.addWidget(self.btn_add)
        v_layout_3.addWidget(self.btn_remove)
        v_layout_3.addStretch()

        v_layout_4.addWidget(self.lbl_available_data)
        v_layout_4.addWidget(self.list_2)

        v_layout_5.addStretch()
        v_layout_5.addWidget(self.btn_list_2_up)
        v_layout_5.addWidget(self.btn_list_2_down)
        v_layout_5.addStretch()

        h_layout.addLayout(v_layout_1)
        h_layout.addLayout(v_layout_2)
        h_layout.addLayout(v_layout_3)
        h_layout.addLayout(v_layout_4)
        h_layout.addLayout(v_layout_5)

        self.widget_frame_1.setLayout(h_layout)

    def frame_2_layout(self):
        self.chb_edge_columns.setChecked(True)
        self.chb_matrix_columns.setChecked(False)
        self.chb_particle_columns.setChecked(False)
        self.chb_hidden_columns.setChecked(False)
        self.chb_flag_1.setChecked(False)
        self.chb_flag_2.setChecked(False)
        self.chb_flag_3.setChecked(False)
        self.chb_flag_4.setChecked(False)

        h_layout = QtWidgets.QHBoxLayout()
        v_layout = QtWidgets.QVBoxLayout()

        v_layout.addStretch()
        v_layout.addWidget(self.lbl_filter)
        v_layout.addWidget(self.chb_edge_columns)
        v_layout.addWidget(self.chb_matrix_columns)
        v_layout.addWidget(self.chb_particle_columns)
        v_layout.addWidget(self.chb_hidden_columns)
        v_layout.addWidget(self.chb_flag_1)
        v_layout.addWidget(self.chb_flag_2)
        v_layout.addWidget(self.chb_flag_3)
        v_layout.addWidget(self.chb_flag_4)
        v_layout.addStretch()

        h_layout.addStretch()
        h_layout.addLayout(v_layout)
        h_layout.addStretch()

        self.widget_frame_2.setLayout(h_layout)

    def frame_3_layout(self):
        self.chb_recalculate_graphs.setChecked(False)
        v_layout = QtWidgets.QVBoxLayout()
        v_layout.addWidget(self.chb_recalculate_graphs)
        v_layout.addStretch()
        self.widget_frame_3.setLayout(v_layout)

    def set_layout(self):
        self.btn_layout.addStretch()
        self.btn_layout.addWidget(self.btn_cancel)
        self.btn_layout.addWidget(self.btn_back)
        self.btn_layout.addWidget(self.btn_next)
        self.btn_layout.addStretch()

        self.frame_0_layout()
        self.stack_layout.addWidget(self.widget_frame_0)

        self.frame_1_layout()
        self.stack_layout.addWidget(self.widget_frame_1)

        self.frame_2_layout()
        self.stack_layout.addWidget(self.widget_frame_2)

        self.frame_3_layout()
        self.stack_layout.addWidget(self.widget_frame_3)

        self.top_layout.addLayout(self.stack_layout)
        self.top_layout.addLayout(self.btn_layout)

        self.setLayout(self.top_layout)

    def btn_list_2_up_trigger(self):
        if self.list_2.currentItem() is not None:
            text = self.list_2.currentItem().text()
            index = self.list_2.currentRow()
            if not index == 0:
                self.list_2.takeItem(self.list_2.row(self.list_2.currentItem()))
                self.list_2.insertItem(index - 1, text)
                self.list_2.setCurrentRow(index - 1)

    def btn_list_2_down_trigger(self):
        if self.list_2.currentItem() is not None:
            text = self.list_2.currentItem().text()
            index = self.list_2.currentRow()
            if not index == self.list_2.count() - 1:
                self.list_2.takeItem(self.list_2.row(self.list_2.currentItem()))
                self.list_2.insertItem(index + 1, text)
                self.list_2.setCurrentRow(index + 1)

    def btn_list_1_up_trigger(self):
        if self.list_1.currentItem() is not None:
            text = self.list_1.currentItem().text()
            index = self.list_1.currentRow()
            if not index == 0:
                self.list_1.takeItem(self.list_1.row(self.list_1.currentItem()))
                self.list_1.insertItem(index - 1, text)
                self.list_1.setCurrentRow(index - 1)

    def btn_list_1_down_trigger(self):
        if self.list_1.currentItem() is not None:
            text = self.list_1.currentItem().text()
            index = self.list_1.currentRow()
            if not index == self.list_1.count() - 1:
                self.list_1.takeItem(self.list_1.row(self.list_1.currentItem()))
                self.list_1.insertItem(index + 1, text)
                self.list_1.setCurrentRow(index + 1)

    def btn_add_item_trigger(self):
        if self.list_2.currentItem() is not None:
            self.list_1.addItem(self.list_2.currentItem().text())
            self.list_2.takeItem(self.list_2.row(self.list_2.currentItem()))

    def btn_remove_item_trigger(self):
        if self.list_1.currentItem() is not None:
            self.list_2.addItem(self.list_1.currentItem().text())
            self.list_1.takeItem(self.list_1.row(self.list_1.currentItem()))

    def btn_next_trigger(self):
        if self.stack_layout.currentIndex() == len(self.stack_layout) - 1:
            self.close()
            self.ui_obj.sys_message('Working...')
            # Prepare filter
            filter_ = {
                'exclude_edge_columns': self.chb_edge_columns.isChecked(),
                'exclude_matrix_columns': self.chb_matrix_columns.isChecked(),
                'exclude_particle_columns': self.chb_particle_columns.isChecked(),
                'exclude_hidden_columns': self.chb_hidden_columns.isChecked(),
                'exclude_flag_1_columns': self.chb_flag_1.isChecked(),
                'exclude_flag_2_columns': self.chb_flag_2.isChecked(),
                'exclude_flag_3_columns': self.chb_flag_3.isChecked(),
                'exclude_flag_4_columns': self.chb_flag_4.isChecked()
            }
            # Prepare keys?
            attr_keys = []
            for i in range(self.list_1.count()):
                attr_keys.append(self.list_1.item(i).text())
            cat_key = self.cmb_nominal_data.currentText()
            self.ui_obj.project_instance.save('temp/temp_instance')
            # Initiate manager
            manager = data_module.VertexDataManager(
                'temp/temp_instance',
                attr_keys,
                filter_=filter_,
                save_filename='temp/temp_model',
                recalc=self.chb_recalculate_graphs.isChecked(),
                category_key=cat_key
            )
            manager.process_data()
            manager.save()
            self.ui_obj.sys_message('Ready.')

            PlotModels(ui_obj=self.ui_obj, model='temp/temp_model')
        else:
            self.stack_layout.setCurrentIndex(self.stack_layout.currentIndex() + 1)

    def btn_back_trigger(self):
        if self.stack_layout.currentIndex() == 0:
            self.close()
        else:
            self.stack_layout.setCurrentIndex(self.stack_layout.currentIndex() - 1)

    def btn_cancel_trigger(self):
        self.close()


class SubGraphWizard(QtWidgets.QDialog):

    def __init__(self, *args, ui_obj=None):
        super().__init__(*args)

        self.ui_obj = ui_obj

        self.setWindowTitle('Get subgraph')

        self.btn_next = QtWidgets.QPushButton('Make')
        self.btn_next.clicked.connect(self.btn_next_trigger)
        self.btn_cancel = QtWidgets.QPushButton('Cancel')
        self.btn_cancel.clicked.connect(self.btn_cancel_trigger)

        self.cmb_export = QtWidgets.QComboBox()
        self.cmb_export.addItem('Mesh centered subgraph')
        self.cmb_export.addItem('Arc centered subgraph')
        self.cmb_export.addItem('Vertex centered subgraph')
        self.cmb_export.currentIndexChanged.connect(self.cmb_export_trigger)
        self.lbl_export = QtWidgets.QLabel('Select subgraph type: ')

        self.sbox_order = QtWidgets.QSpinBox()
        self.sbox_order.setMinimum(0)
        self.sbox_order.setMaximum(4)
        self.sbox_order.valueChanged.connect(self.cmb_export_trigger)
        self.lbl_order = QtWidgets.QLabel('Select subgraph order: ')

        self.lbl_indices = QtWidgets.QLabel('Select defining vertex indices:')
        self.lbl_i = QtWidgets.QLabel('\ti (The previously selected column): ')
        self.lbl_j = QtWidgets.QLabel('\tj (The currently selected column): ')
        self.sbox_i = QtWidgets.QSpinBox()
        self.sbox_i.setMaximum(self.ui_obj.project_instance.num_columns - 1)
        self.sbox_i.setMinimum(0)
        if not self.ui_obj.previous_selected_column == -1:
            self.sbox_i.setValue(self.ui_obj.previous_selected_column)
        else:
            self.sbox_i.setValue(0)
        self.sbox_i.valueChanged.connect(self.cmb_export_trigger)
        self.sbox_j = QtWidgets.QSpinBox()
        self.sbox_j.setMaximum(self.ui_obj.project_instance.num_columns - 1)
        self.sbox_j.setMinimum(0)
        if not self.ui_obj.selected_column == -1:
            self.sbox_j.setValue(self.ui_obj.selected_column)
        else:
            self.sbox_j.setValue(1)
        self.sbox_j.valueChanged.connect(self.cmb_export_trigger)

        self.lbl_definition = QtWidgets.QLabel('Definition: ')

        top_layout = QtWidgets.QVBoxLayout()

        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.lbl_export)
        layout.addWidget(self.cmb_export)
        top_layout.addLayout(layout)

        layout = QtWidgets.QHBoxLayout()
        layout.addWidget(self.lbl_order)
        layout.addWidget(self.sbox_order)
        top_layout.addLayout(layout)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.lbl_indices)

        layout_1 = QtWidgets.QHBoxLayout()
        layout_1.addWidget(self.lbl_i)
        layout_1.addWidget(self.sbox_i)
        layout.addLayout(layout_1)

        layout_2 = QtWidgets.QHBoxLayout()
        layout_2.addWidget(self.lbl_j)
        layout_2.addWidget(self.sbox_j)
        layout.addLayout(layout_2)

        top_layout.addLayout(layout)

        top_layout.addWidget(self.lbl_definition)

        btn_layout = QtWidgets.QHBoxLayout()
        btn_layout.addStretch()
        btn_layout.addWidget(self.btn_cancel)
        btn_layout.addWidget(self.btn_next)
        btn_layout.addStretch()
        top_layout.addLayout(btn_layout)

        self.cmb_export_trigger()

        self.setLayout(top_layout)
        self.exec_()

    def cmb_export_trigger(self):
        if self.cmb_export.currentText() == 'Vertex centered subgraph':
            self.lbl_i.hide()
            self.sbox_i.hide()
        else:
            self.lbl_i.show()
            self.sbox_i.show()

        if self.cmb_export.currentText() == 'Vertex centered subgraph':
            self.lbl_definition.setText('Definition: H_{}^{}({})'.format(
                'vertex',
                self.sbox_order.value(),
                self.sbox_j.value()
            ))
        elif self.cmb_export.currentText() == 'Mesh centered subgraph':
            self.lbl_definition.setText('Definition: H_{}^{}({}, {})'.format(
                'mesh',
                self.sbox_order.value(),
                self.sbox_i.value(),
                self.sbox_j.value()
            ))
        else:
            self.lbl_definition.setText('Definition: H_{}^{}({}, {})'.format(
                'arc',
                self.sbox_order.value(),
                self.sbox_i.value(),
                self.sbox_j.value()
            ))

    def btn_next_trigger(self):
        self.close()

    def btn_cancel_trigger(self):
        self.close()

