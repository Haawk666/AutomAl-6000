# By Haakon Tvedt @ NTNU
# Contributors:
"""Module container for low-level custom GUI-elements"""

# Program imports:
import GUI_settings
import GUI
# External imports:
from PyQt5 import QtWidgets, QtGui, QtCore
import numpy as np
import logging
# Instantiate logger:
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# ----------
# Graphic elements:
# ----------


class InteractiveColumn(QtWidgets.QGraphicsEllipseItem):

    """A general interactive graphical element that is meant to represent atomic columns in the GUI.

    Inherits PyQt5.QtWidgets.QGraphicsEllipseItem()
    """

    def __init__(self, ui_obj=None, vertex=None, scale_factor=1, r=5, movable=True, selectable=True):

        """Initialize with optional reference to a MainUI object.

        parameters
        ----------
        ui_obj : GUI.MainUI, optional
            Reference MainUI object.
        i : int, optional
            An integer that is a reference to the index of the relative vertex in
            ui_obj.project_instance.graph.vertices[i]. If i == -1, it is interpreted as no reference.
        """

        self.ui_obj = ui_obj
        if vertex is not None:
            self.r = r
            self.i = vertex.i
            self.center_coor = vertex.im_pos()
        else:
            self.r = r
            self.i = -2
            self.center_coor = (0, 0)
        self.scale_factor = scale_factor
        self.vertex = vertex

        super().__init__(0, 0, 2 * self.r, 2 * self.r)

        if movable:
            offset = self.boundingRect().center()
            self.moveBy(self.scale_factor * self.center_coor[0] - offset.x(), self.scale_factor * self.center_coor[1] - offset.y())
        if selectable:
            self.setFlag(QtWidgets.QGraphicsItem.ItemIsSelectable)

    def get_vertex_pos(self):
        offset = self.boundingRect().center()
        pos = (self.x() + offset.x(), self.y() + offset.y())
        return pos

    def mouseReleaseEvent(self, event: 'QtWidgets.QGraphicsEllipseItem.mouseReleaseEvent'):
        """Pass a mouse release event on to the ui_obj reference object"""
        self.ui_obj.column_selected(self.i)

    def mousePressEvent(self, *args, **kwargs):
        pass


class InteractivePosColumn(InteractiveColumn):

    def __init__(self, ui_obj=None, vertex=None, scale_factor=1, movable=True, r=16):
        """Initialize a positional interactive column.

        Inherits GUI_custom_components.InteractiveColumn. Is used to highlight atomic positions."""
        super().__init__(ui_obj=ui_obj, vertex=vertex, scale_factor=scale_factor, movable=movable, r=r)

        self.selected_pen = GUI_settings.pen_selected_2
        self.selected_brush = GUI_settings.brush_selected_2
        self.unselected_pen = GUI_settings.pen_atom_pos
        self.unselected_brush = GUI_settings.brush_atom_pos
        self.hidden_pen = GUI_settings.pen_atom_pos_hidden
        self.hidden_brush = GUI_settings.brush_atom_pos_hidden
        self.set_style()

    def set_style(self):
        """Set the appearance of the shape"""
        if self.i == self.ui_obj.selected_column:
            self.setPen(self.selected_pen)
            self.setBrush(self.selected_brush)
        else:
            if self.vertex.show_in_overlay:
                self.setPen(self.hidden_pen)
                self.setBrush(self.hidden_brush)
            else:
                self.setPen(self.unselected_pen)
                self.setBrush(self.hidden_brush)


class InteractiveOverlayColumn(InteractiveColumn):

    def __init__(
            self, ui_obj=None, vertex=None, scale_factor=1, movable=True, selectable=True, r=16,
            pen=QtGui.QPen(QtCore.Qt.red), brush=QtGui.QBrush(QtCore.Qt.red)
    ):
        """Initialize a positional interactive column.

        Inherits GUI_custom_components.InteractiveColumn. Is used to highlight atomic positions."""
        super().__init__(ui_obj=ui_obj, vertex=vertex, scale_factor=scale_factor, movable=movable, selectable=selectable, r=r)

        self.pen = pen
        self.brush = brush
        self.selected_pen = QtGui.QPen(QtCore.Qt.darkCyan)
        self.selected_pen.setWidth(6)
        self.selected_brush = QtGui.QBrush(QtCore.Qt.darkCyan)

        self.set_style()

    def set_style(self):
        """Set the appearance of the shape"""
        if self.i == self.ui_obj.selected_column:
            self.setPen(self.selected_pen)
            self.setBrush(self.selected_brush)
        else:
            self.setPen(self.pen)
            self.setBrush(self.brush)
        if self.vertex is not None:
            if not self.vertex.show_in_overlay:
                self.hide()
            else:
                self.show()
        else:
            self.show()


class InteractiveGraphColumn(InteractiveColumn):

    def __init__(self, ui_obj=None, vertex=None, scale_factor=1, r=16):
        """Initialize a positional interactive column.

        Inherits GUI_custom_components.InteractiveColumn. Is used to highlight atomic positions."""
        super().__init__(ui_obj=ui_obj, vertex=vertex, scale_factor=scale_factor, movable=True, selectable=True, r=r)

        self.selected_pen = GUI_settings.pen_selected_1
        self.unselected_pen = GUI_settings.pen_graph
        self.level_0_brush = GUI_settings.brush_graph_0
        self.level_1_brush = GUI_settings.brush_graph_1

        self.set_style()

    def set_style(self):
        """Set the appearance of the shape"""
        if self.i == self.ui_obj.selected_column:
            self.setPen(self.selected_pen)
        else:
            self.setPen(self.unselected_pen)
        if self.vertex.zeta == 0:
            self.setBrush(self.level_0_brush)
        else:
            self.setBrush(self.level_1_brush)


class Arrow(QtWidgets.QGraphicsItemGroup):

    def __init__(self, *args, i=0, j=1, p1=(0, 0), p2=(1, 1), r=1, scale_factor=1, dual_arc=False, co_planar=False):
        super().__init__(*args)

        self.inconsistent_pen = GUI_settings.pen_inconsistent_edge
        self.dislocation_pen = GUI_settings.pen_dislocation_edge
        self.normal_pen = GUI_settings.pen_edge

        self.dual_arc = dual_arc
        self.co_planar = co_planar
        self.scale_factor = scale_factor
        self.r = r
        self.p1 = p1
        self.p2 = p2
        self.arrow = None, None
        self.i = i
        self.j = j

        self.make_arrow_obj()
        self.set_style()

    def set_style(self):
        if self.dual_arc and not self.co_planar:
            self.childItems()[0].setPen(self.normal_pen)
            self.childItems()[0].show()
            self.childItems()[1].hide()
        elif self.dual_arc and self.co_planar:
            self.childItems()[0].setPen(self.dislocation_pen)
            self.childItems()[0].show()
            self.childItems()[1].setPen(self.dislocation_pen)
            self.childItems()[1].show()
        else:
            self.childItems()[0].setPen(self.inconsistent_pen)
            self.childItems()[0].show()
            self.childItems()[1].setPen(self.inconsistent_pen)
            self.childItems()[1].show()

    def make_arrow_obj(self):

        r_2 = QtCore.QPointF(self.scale_factor * self.p2[0], self.scale_factor * self.p2[1])
        r_1 = QtCore.QPointF(self.scale_factor * self.p1[0], self.scale_factor * self.p1[1])

        r_vec = r_2 - r_1
        r_mag = np.sqrt((r_2.x() - r_1.x()) ** 2 + (r_2.y() - r_1.y()) ** 2)
        factor = self.r / (r_mag * 2)

        k_2 = r_1 + (1 - factor) * r_vec

        theta = np.pi / 4

        l_3 = - factor * QtCore.QPointF(r_vec.x() * np.cos(theta) + r_vec.y() * np.sin(theta), - r_vec.x() * np.sin(theta) + r_vec.y() * np.cos(theta))
        l_3 = k_2 + l_3
        l_4 = - factor * QtCore.QPointF(r_vec.x() * np.cos(-theta) + r_vec.y() * np.sin(-theta), - r_vec.x() * np.sin(-theta) + r_vec.y() * np.cos(-theta))
        l_4 = k_2 + l_4

        tri_2 = (k_2, l_3, l_4)

        poly_2 = QtGui.QPolygonF(tri_2)

        line = QtWidgets.QGraphicsLineItem(self.scale_factor * self.p1[0],
                                           self.scale_factor * self.p1[1],
                                           self.scale_factor * self.p2[0],
                                           self.scale_factor * self.p2[1])
        head_2 = QtWidgets.QGraphicsPolygonItem(poly_2)

        self.addToGroup(line)
        self.addToGroup(head_2)
        self.setZValue(-1)


class DistanceArrow(QtWidgets.QGraphicsItemGroup):

    def __init__(self, *args, color=None, i=0, j=1, p1=(0, 0), p2=(1, 1), r=1, scale_factor=1):
        super().__init__(*args)

        self.inconsistent_pen = GUI_settings.pen_inconsistent_edge
        self.dislocation_pen = GUI_settings.pen_dislocation_edge
        self.normal_pen = GUI_settings.pen_edge

        if color is None:
            self.color = QtGui.QColor(255, 255, 255, 255)
        else:
            self.color = QtGui.QColor(color[0], color[1], color[2], color[3])
        self.pen = QtGui.QPen(self.color)
        self.pen.setWidth(3)

        self.scale_factor = scale_factor
        self.r = r
        self.p1 = p1
        self.p2 = p2
        self.i = i
        self.j = j

        self.make_arrow_obj()
        self.set_style()

    def set_style(self):

            self.childItems()[0].setPen(self.pen)
            self.childItems()[0].show()
            self.show()

    def make_arrow_obj(self):

        line = QtWidgets.QGraphicsLineItem(self.scale_factor * self.p1[0],
                                           self.scale_factor * self.p1[1],
                                           self.scale_factor * self.p2[0],
                                           self.scale_factor * self.p2[1])

        self.addToGroup(line)
        self.setZValue(-1)


class ScaleBar(QtWidgets.QGraphicsItemGroup):

    def __init__(self, *args, length=2, scale=5, r=10, height=512, unit='nm'):
        """Initialize a positional interactive column.

        Inherits GUI_custom_components.InteractiveColumn. Is used to highlight atomic positions."""
        super().__init__(*args)

        self.length = length
        self.height = height
        self.pico_scale = scale
        self.aa_scale = scale / 100
        self.nano_scale = scale / 1000
        self.micro_scale = scale / 1000000
        self.unit = unit
        if unit == 'pm':
            self.relative_scale = self.pico_scale
        elif unit == 'Å':
            self.relative_scale = self.aa_scale
        elif unit == 'nm':
            self.relative_scale = self.nano_scale
        elif unit == '@\mu@m':
            self.relative_scale = self.micro_scale
        else:
            self.relative_scale = self.nano_scale
            self.unit = 'nm'
        self.r = r
        self.color_brush = GUI_settings.brush_white
        self.make()

    def make(self):
        p1 = 0, 0
        p2 = self.length / self.relative_scale - GUI_settings.pen_scalebar.width(), 0
        line = QtWidgets.QGraphicsLineItem(p1[0], p1[1], p2[0], p2[1])
        line.setPen(GUI_settings.pen_scalebar)

        text = QtWidgets.QGraphicsSimpleTextItem()
        text.setText('{} {}'.format(self.length, self.unit))
        text.setFont(GUI_settings.font_scalebar)
        text.setPen(GUI_settings.white_pen)
        text.setBrush(GUI_settings.brush_white)
        rect = text.boundingRect()
        text.setX(p2[0] / 2 - rect.width() / 2)
        text.setY(- rect.height() - GUI_settings.pen_scalebar.width())

        self.addToGroup(line)
        self.addToGroup(text)

        self.moveBy(self.r * 6, self.height - self.r * 6)
        self.setFlag(QtWidgets.QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QtWidgets.QGraphicsItem.ItemIsSelectable, True)


class Legend(QtWidgets.QGraphicsItemGroup):

    def __init__(self, *args, ui_obj=None):
        super().__init__(*args)

        self.ui_obj = ui_obj
        self.transparent = False
        if self.ui_obj.overlay_settings['legend_background'] == 'Black':
            self.backgroundcolor = QtGui.QColor(0, 0, 0)
            self.textcolor = QtGui.QColor(255, 255, 255)
        elif self.ui_obj.overlay_settings['legend_background'] == 'White':
            self.backgroundcolor = QtGui.QColor(255, 255, 255)
            self.textcolor = QtGui.QColor(0, 0, 0)
        elif self.ui_obj.overlay_settings['legend_background'] == 'Gray':
            self.backgroundcolor = QtGui.QColor(120, 120, 120)
            self.textcolor = QtGui.QColor(0, 0, 0)
        elif self.ui_obj.overlay_settings['legend_background'] == 'Transparent':
            self.transparent = True
            self.backgroundcolor = QtGui.QColor(0, 0, 0)
            self.textcolor = QtGui.QColor(255, 255, 255)
        else:
            self.backgroundcolor = QtGui.QColor(0, 0, 0)
            self.textcolor = QtGui.QColor(255, 255, 255)
        self.font = QtGui.QFont('Helvetica [Cronyx]', 22)
        self.pen = QtGui.QPen(self.textcolor)
        self.brush = QtGui.QBrush(self.textcolor)
        if self.ui_obj.project_instance is not None:
            self.species_dict = self.ui_obj.project_instance.species_dict
            self.scale = self.ui_obj.project_instance.scale
        else:
            self.species_dict = self.ui_obj.core.Project.default_species_dict
            self.scale = 1
        self.pen.setWidth(1)
        self.make()
        if not self.ui_obj.control_window.widgets['overlay']['chb_legend']['widget'].isChecked():
            self.hide()
        else:
            self.show()

    def make(self):

        h_margin = 0
        v_margin = 0
        separation_margin = 10
        head_line_width = 3

        categories = {}
        max_r = 0
        max_text_width = 0

        if self.ui_obj.overlay_settings['display_mode'] == 'atomic_species':
            for key, item in self.species_dict['atomic_species'].items():

                r = (item['atomic_radii'] - 15) / self.scale

                label = QtWidgets.QGraphicsSimpleTextItem()
                label.setText(key)
                label.setFont(self.font)
                label.setPen(self.pen)
                label.setBrush(self.brush)

                q_color = QtGui.QColor(
                    item['color'][0],
                    item['color'][1],
                    item['color'][2]
                )
                k_color = QtGui.QColor(0, 0, 0)
                pen = QtGui.QPen(q_color)
                pen.setWidth(int(30 / self.scale))

                object_0 = InteractiveOverlayColumn(
                    ui_obj=self.ui_obj,
                    selectable=False,
                    r=self.ui_obj.overlay_settings['overlay_radii'] * r,
                    pen=pen,
                    brush=QtGui.QBrush(q_color)
                )

                object_1 = InteractiveOverlayColumn(
                    ui_obj=self.ui_obj,
                    selectable=False,
                    r=self.ui_obj.overlay_settings['overlay_radii'] * r,
                    pen=pen,
                    brush=QtGui.QBrush(k_color)
                )

                categories[key] = {
                    'label': label,
                    'item_0': object_0,
                    'item_1': object_1,
                }

                if r > max_r:
                    max_r = r
                text_width = label.boundingRect().width()
                if text_width > max_text_width:
                    max_text_width = text_width
        else:
            for key, item in self.species_dict['advanced_species'].items():
                r = (self.species_dict['atomic_species'][item['atomic_species']]['atomic_radii'] - 15) / self.scale

                label = QtWidgets.QGraphicsSimpleTextItem()
                label.setText(key)
                label.setFont(self.font)
                label.setPen(self.pen)
                label.setBrush(self.brush)

                q_color = QtGui.QColor(
                    item['color'][0],
                    item['color'][1],
                    item['color'][2]
                )
                k_color = QtGui.QColor(0, 0, 0)
                pen = QtGui.QPen(q_color)
                pen.setWidth(int(30 / self.scale))

                object_0 = InteractiveOverlayColumn(
                    ui_obj=self.ui_obj,
                    selectable=False,
                    r=self.ui_obj.overlay_settings['overlay_radii'] * r,
                    pen=pen,
                    brush=QtGui.QBrush(q_color)
                )

                object_1 = InteractiveOverlayColumn(
                    ui_obj=self.ui_obj,
                    selectable=False,
                    r=self.ui_obj.overlay_settings['overlay_radii'] * r,
                    pen=pen,
                    brush=QtGui.QBrush(k_color)
                )

                categories[key] = {
                    'label': label,
                    'item_0': object_0,
                    'item_1': object_1
                }

                if r > max_r:
                    max_r = r
                text_width = label.boundingRect().width()
                if text_width > max_text_width:
                    max_text_width = text_width

        max_r = 2 * max_r

        zeta_text = QtWidgets.QGraphicsSimpleTextItem()
        zeta_text.setText('zeta:')
        zeta_text.setFont(self.font)
        zeta_text.setPen(self.pen)
        zeta_text.setBrush(self.brush)
        zeta_text.setX(h_margin)
        zeta_text.setY(v_margin)
        self.addToGroup(zeta_text)
        text_width = zeta_text.boundingRect().width()
        if text_width > max_text_width:
            max_text_width = text_width

        zeta_text_0 = QtWidgets.QGraphicsSimpleTextItem()
        zeta_text_0.setText('0')
        zeta_text_0.setFont(self.font)
        zeta_text_0.setPen(self.pen)
        zeta_text_0.setBrush(self.brush)
        zeta_text_0.setX(separation_margin + h_margin + max_text_width)
        zeta_text_0.setY(v_margin)
        self.addToGroup(zeta_text_0)

        zeta_text_1 = QtWidgets.QGraphicsSimpleTextItem()
        zeta_text_1.setText('1')
        zeta_text_1.setFont(self.font)
        zeta_text_1.setPen(self.pen)
        zeta_text_1.setBrush(self.brush)
        zeta_text_1.setX(separation_margin + 2 * h_margin + max_text_width + max_r)
        zeta_text_1.setY(v_margin)
        self.addToGroup(zeta_text_1)

        text_height = zeta_text.boundingRect().height()

        line = QtWidgets.QGraphicsLineItem(
            0,
            v_margin + text_height + head_line_width / 2,
            2 * h_margin + max_text_width + separation_margin + 2 * max_r,
            v_margin + text_height + head_line_width / 2
        )
        line_pen = QtGui.QPen(self.textcolor)
        line_pen.setWidth(head_line_width)
        line.setPen(line_pen)
        self.addToGroup(line)

        counter = 0
        for item in categories.values():

            text = item['label']
            text.setX(h_margin)
            text.setY(2 * v_margin + text_height + counter * (max_r + v_margin) + head_line_width)
            self.addToGroup(text)

            object_0 = item['item_0']
            r = 2 * object_0.rect().height()
            object_0.setX(separation_margin + h_margin + max_text_width)
            object_0.setY(2 * v_margin + text_height + counter * (max_r + v_margin) + (max_r - r) / 2 + head_line_width)
            self.addToGroup(object_0)

            object_1 = item['item_1']
            object_1.setX(separation_margin + 2 * h_margin + max_text_width + max_r)
            object_1.setY(2 * v_margin + text_height + counter * (max_r + v_margin) + (max_r - r) / 2 + head_line_width)
            self.addToGroup(object_1)

            counter += 1

        background = QtWidgets.QGraphicsRectItem(self.boundingRect())
        background.setBrush(QtGui.QBrush(self.backgroundcolor))
        background.setZValue(-1)
        if not self.transparent:
            self.addToGroup(background)

        self.moveBy(30, 30)
        self.setFlag(QtWidgets.QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QtWidgets.QGraphicsItem.ItemIsSelectable, True)
        self.setZValue(100)


class HeatLegend(QtWidgets.QGraphicsItemGroup):

    def __init__(self, *args, ui_obj=None, min=0, max=1, margin=20, width=40):
        super().__init__(*args)

        self.ui_obj = ui_obj

        self.min = min
        self.max = max

        box = QtWidgets.QGraphicsRectItem(
            self.ui_obj.project_instance.im_width + margin, 0, width, self.ui_obj.project_instance.im_height
        )
        gradient_ = QtGui.QLinearGradient()
        gradient_.setStart(0.0, self.ui_obj.project_instance.im_height)
        gradient_.setFinalStop(0.0, 0.0)
        box.setBrush(gradient_)

        if GUI_settings.theme == 'dark':
            textcolor = QtGui.QColor(255, 255, 255)
        else:
            textcolor = QtGui.QColor(0, 0, 0)
        pen = QtGui.QPen(textcolor)
        brush = QtGui.QBrush(textcolor)
        font = QtGui.QFont('Helvetica [Cronyx]', 18)

        text_min = QtWidgets.QGraphicsSimpleTextItem()
        text_min.setText('{:.3f}'.format(self.min))
        text_min.setPen(pen)
        text_min.setFont(font)
        text_min.setBrush(brush)
        text_height = text_min.boundingRect().height()
        text_min.setX(self.ui_obj.project_instance.im_width + margin + width + 10)
        text_min.setY(self.ui_obj.project_instance.im_height - text_height / 2)

        text_max = QtWidgets.QGraphicsSimpleTextItem()
        text_max.setText('{:.3f}'.format(self.max))
        text_max.setPen(pen)
        text_max.setFont(font)
        text_max.setBrush(brush)
        text_height = text_max.boundingRect().height()
        text_max.setX(self.ui_obj.project_instance.im_width + margin + width + 10)
        text_max.setY(0 - text_height / 2)

        text_mid = QtWidgets.QGraphicsSimpleTextItem()
        text_mid.setText('{:.3f}'.format(self.min + (self.max - self.min) / 2))
        text_mid.setPen(pen)
        text_mid.setFont(font)
        text_mid.setBrush(brush)
        text_height = text_mid.boundingRect().height()
        text_mid.setX(self.ui_obj.project_instance.im_width + margin + width + 10)
        text_mid.setY(self.ui_obj.project_instance.im_height / 2 - text_height / 2)

        self.addToGroup(box)
        self.addToGroup(text_min)
        self.addToGroup(text_mid)
        self.addToGroup(text_max)

        self.setFlag(QtWidgets.QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QtWidgets.QGraphicsItem.ItemIsSelectable, True)


class MeshDetail(QtWidgets.QGraphicsItemGroup):

    def __init__(self, *args, mesh=None):
        super().__init__(*args)
        self.mesh = mesh
        if self.mesh is not None:
            self.make()
        self.setFlag(QtWidgets.QGraphicsItem.ItemIsMovable, True)
        self.setFlag(QtWidgets.QGraphicsItem.ItemIsSelectable, True)

    def make(self):
        num_corners = self.mesh.order
        text = QtWidgets.QGraphicsSimpleTextItem()
        text.setText('{}'.format(num_corners))
        text.setFont(GUI_settings.font_mesh_details)
        if num_corners == 4:
            text.setPen(GUI_settings.pen_boarder)
        else:
            text.setPen(GUI_settings.pen_skinny_red)
        text.setBrush(GUI_settings.brush_black)
        rect = text.boundingRect()
        text.setX(2 * self.mesh.cm[0] - rect.width() / 2)
        text.setY(2 * self.mesh.cm[1] - rect.height() / 2)

        self.addToGroup(text)

        self.setFlag(QtWidgets.QGraphicsItem.ItemIsMovable, False)
        self.setFlag(QtWidgets.QGraphicsItem.ItemIsSelectable, True)

    def mouseReleaseEvent(self, event: 'QtWidgets.QGraphicsItemGroup.mouseReleaseEvent'):
        print('{}'.format(self.mesh.mesh_index))
        # super().__init__()


class RgbaSelector(QtWidgets.QGroupBox):

    def __init__(self, *args, r=100, g=100, b=100, a=100):
        super().__init__(*args)

        self.r_box = QtWidgets.QSpinBox()
        self.r_box.setMaximum(255)
        self.r_box.setMinimum(0)
        if 0 <= r <= 255:
            self.r_box.setValue(r)
        else:
            self.r_box.setValue(100)

        self.g_box = QtWidgets.QSpinBox()
        self.g_box.setMaximum(255)
        self.g_box.setMinimum(0)
        if 0 <= g <= 255:
            self.g_box.setValue(g)
        else:
            self.g_box.setValue(100)

        self.b_box = QtWidgets.QSpinBox()
        self.b_box.setMaximum(255)
        self.b_box.setMinimum(0)
        if 0 <= b <= 255:
            self.b_box.setValue(b)
        else:
            self.b_box.setValue(100)

        self.a_box = QtWidgets.QSpinBox()
        self.a_box.setMaximum(255)
        self.a_box.setMinimum(0)
        if 0 <= a <= 255:
            self.a_box.setValue(a)
        else:
            self.a_box.setValue(100)

        self.lbl_r = QtWidgets.QLabel('R:')
        self.lbl_g = QtWidgets.QLabel('G:')
        self.lbl_b = QtWidgets.QLabel('B:')
        self.lbl_a = QtWidgets.QLabel('A:')

        layout = QtWidgets.QGridLayout()

        layout.addWidget(self.lbl_r, 0, 0)
        layout.addWidget(self.lbl_g, 0, 1)
        layout.addWidget(self.lbl_b, 0, 2)
        layout.addWidget(self.lbl_a, 0, 3)
        layout.addWidget(self.r_box, 1, 0)
        layout.addWidget(self.g_box, 1, 1)
        layout.addWidget(self.b_box, 1, 2)
        layout.addWidget(self.a_box, 1, 3)

        self.setLayout(layout)

        self.setTitle('')

    def get_triple(self):
        triple = (
            self.r_box.value(),
            self.g_box.value(),
            self.b_box.value()
        )
        return triple

    def set_triple(self, triple):
        self.r_box.setValue(triple[0])
        self.g_box.setValue(triple[1])
        self.b_box.setValue(triple[2])


class CenterBox(QtWidgets.QGroupBox):

    def __init__(self, *args, widget=None):
        super().__init__(*args)

        self.widget = widget

        h_layout = QtWidgets.QHBoxLayout()
        h_layout.addStretch()
        h_layout.addWidget(self.widget)
        h_layout.addStretch()

        v_layout = QtWidgets.QVBoxLayout()
        v_layout.addStretch()
        v_layout.addLayout(h_layout)
        v_layout.addStretch()

        self.setLayout(v_layout)

        self.setStyleSheet("border:0;")

        self.widget.setStyleSheet('')


class HCenterBox(QtWidgets.QGroupBox):

    def __init__(self, *args, widget=None):
        super().__init__(*args)

        self.widget = widget

        h_layout = QtWidgets.QHBoxLayout()
        h_layout.addStretch()
        h_layout.addWidget(self.widget)
        h_layout.addStretch()

        self.setLayout(h_layout)

        self.setFlat(True)
        self.setStyleSheet("border:0;")


class VCenterBox(QtWidgets.QGroupBox):

    def __init__(self, *args, widget=None):
        super().__init__(*args)

        self.widget = widget

        v_layout = QtWidgets.QVBoxLayout()
        v_layout.addStretch()
        v_layout.addWidget(self.widget)
        v_layout.addStretch()

        self.setLayout(v_layout)

        self.setFlat(True)
        self.setStyleSheet("border:0;")


# ----------
# Convenience re-implementations:
# ----------


class SmallButton(QtWidgets.QPushButton):

    def __init__(self, *args, trigger_func=None):
        super().__init__(*args)

        self.trigger_func = trigger_func
        self.clicked.connect(trigger_func)
        self.setMaximumHeight(15)
        self.setMaximumWidth(50)
        self.setFont(GUI_settings.font_tiny)


class MediumButton(QtWidgets.QPushButton):

    def __init__(self, *args, trigger_func=None):
        super().__init__(*args)

        self.clicked.connect(trigger_func)
        self.setMaximumHeight(20)
        self.setMaximumWidth(200)
        self.setFont(GUI_settings.font_tiny)


class SetButton(QtWidgets.QPushButton):

    def __init__(self, obj, trigger_func=None):
        super().__init__('Set', obj)

        self.trigger_func = trigger_func
        self.clicked.connect(trigger_func)
        self.setMaximumHeight(15)
        self.setMaximumWidth(30)
        self.setFont(GUI_settings.font_tiny)


class SetButtonLayout(QtWidgets.QHBoxLayout):

    def __init__(self, *args, obj=None, trigger_func=None, text='', value=''):
        super().__init__(*args)
        self.text = text
        self.value = value
        self.label = QtWidgets.QLabel('{}: {}'.format(text, value))
        self.btn = SetButton(obj, trigger_func)

        self.addWidget(self.btn)
        self.addWidget(self.label)
        self.addStretch()

    def set_value(self, value):
        self.value = value
        self.label.setText('{}: {}'.format(self.text, self.value))

    def setToolTip(self, tooltip):
        self.label.setToolTip(tooltip)
        self.btn.setToolTip(tooltip)

    def setDisabled(self, state):
        self.btn.setDisabled(state)


class MoveControls(QtWidgets.QHBoxLayout):

    def __init__(self, *args, text='', trigger_chb=None, trigger_1=None, trigger_2=None):
        super().__init__(*args)

        self.chb = CheckBox(label=text, trigger_func=trigger_chb)
        self.btn_1 = SmallButton('Accept', trigger_func=trigger_1)
        self.btn_1.setDisabled(True)
        self.btn_2 = SmallButton('Cancel', trigger_func=trigger_2)
        self.btn_2.setDisabled(True)

        self.addWidget(self.chb)
        self.addWidget(self.btn_1)
        self.addWidget(self.btn_2)
        self.addStretch()

    def setToolTip(self, tooltip):
        self.chb.setToolTip(tooltip)
        self.btn_1.setToolTip(tooltip)
        self.btn_2.setToolTip(tooltip)

    def isChecked(self):
        return self.chb.isChecked()

    def setChecked(self, state):
        self.chb.setChecked(state)

    def setDisabled(self, on):
        self.btn_1.setDisabled(on)
        self.btn_2.setDisabled(on)

    def blockSignals(self, state):
        self.chb.blockSignals(state)


class CheckBox(QtWidgets.QCheckBox):

    def __init__(self, label='', trigger_func=None, default_state=False):
        super().__init__(label)
        self.default_state = default_state
        self.setChecked(self.default_state)
        if trigger_func is not None:
            self.toggled.connect(trigger_func)

    def set_state_no_trigger(self, state):
        self.blockSignals(True)
        self.setChecked(state)
        self.blockSignals(False)


class Label(QtWidgets.QLabel):

    def __init__(self, text='', value=''):
        self.text = text
        self.value = value
        super().__init__('{}: {}'.format(self.text, self.value))

    def set_value(self, value):
        self.value = value
        self.setText('{}: {}'.format(self.text, self.value))


class StudioButton(QtWidgets.QWidget):

    def __init__(self, *args, ui_obj=None, nominal_attribute='', category=''):
        super().__init__(*args)

        self.ui_obj = ui_obj
        self.nominal_attribute = nominal_attribute
        self.category = category

        self.lbl_category = QtWidgets.QLabel('{} = {}'.format(self.nominal_attribute, self.category))

        self.btn_hide = QtWidgets.QPushButton('Hide')
        self.btn_hide.setCheckable(True)
        self.btn_hide.setChecked(False)
        self.btn_hide.toggled.connect(self.hide_func)

        self.btn_show = QtWidgets.QPushButton('Show')
        self.btn_show.setCheckable(True)
        self.btn_show.setChecked(False)
        self.btn_show.toggled.connect(self.show_func)

        btn_layout = QtWidgets.QVBoxLayout()
        btn_layout.addWidget(self.btn_hide)
        btn_layout.addWidget(self.btn_show)

        layout = QtWidgets.QHBoxLayout()
        layout.addLayout(btn_layout)
        layout.addWidget(self.lbl_category)
        layout.addStretch()

    def hide_func(self):
        if self.ui_obj.project_instance is not None:
            for vertex in self.ui_obj.project_instance.graph:
                if eval('vertex.{}'.format(self.nominal_attribute)) == self.category:
                    vertex.show_in_overlay = False

    def show_func(self):
        if self.ui_obj.project_instance is not None:
            for vertex in self.ui_obj.project_instance.graph:
                if eval('vertex.{}'.format(self.nominal_attribute)) == self.category:
                    vertex.show_in_overlay = True


class StudioView(QtWidgets.QWidget):

    def __init__(self, *args, ui_obj=None):
        super().__init__(*args)

        self.ui_obj = ui_obj

        if self.ui_obj.project_instance is None:
            self.species_dict = GUI.core.Project.default_species_dict
        else:
            self.species_dict = self.ui_obj.project_instance.species_dict

        self.channels_list = []

        self.refresh_layout()

    def populate(self):

        self.channels_list = []
        if self.ui_obj.overlay_settings['display_mode'] == 'advanced_species':
            for category in GUI.core.graph_2.AtomicGraph.get_advanced_species_list(None, species_dict=self.species_dict):
                self.channels_list.append(
                    {
                        'nominal_attribute': 'advanced_species',
                        'category': category,
                        'sbtn': StudioButton(
                            ui_obj=self.ui_obj,
                            nominal_attribute='advanced_species',
                            category=category
                        )
                    }
                )
        else:
            for category in GUI.core.graph_2.AtomicGraph.get_atomic_species_list(None, species_dict=self.species_dict):
                self.channels_list.append(
                    {
                        'nominal_attribute': 'atomic_species',
                        'category': category,
                        'sbtn': StudioButton(
                            ui_obj=self.ui_obj,
                            nominal_attribute='atomic_species',
                            category=category
                        )
                    }
                )
        for nominal_attribute in GUI.core.graph_2.Vertex.nominal_attributes:
            if nominal_attribute == 'advanced_species' or \
                    nominal_attribute == 'atomic_species' or \
                    nominal_attribute == 'i' or \
                    nominal_attribute == 'flag_1' or \
                    nominal_attribute == 'flag_2' or \
                    nominal_attribute == 'flag_3' or \
                    nominal_attribute == 'flag_4' or \
                    nominal_attribute == 'flag_5' or \
                    nominal_attribute == 'flag_6' or \
                    nominal_attribute == 'flag_7' or \
                    nominal_attribute == 'flag_8' or \
                    nominal_attribute == 'flag_9':
                pass
            else:
                range_list = set()
                if self.ui_obj.project_instance is not None:
                    for vertex in self.ui_obj.project_instance.graph.vertices:
                        range_list.add(eval('vertex.{}'.format(nominal_attribute)))
                    range_list = list(range_list)
                    range_list.sort()
                for value in range_list:
                    self.channels_list.append(
                        {
                            'nominal_attribute': nominal_attribute,
                            'category': value,
                            'sbtn': StudioButton(
                                ui_obj=self.ui_obj,
                                nominal_attribute=nominal_attribute,
                                category=value
                            )
                        }
                    )

    def refresh_layout(self):
        self.populate()
        layout = QtWidgets.QVBoxLayout()
        for channel in self.channels_list:
            layout.addWidget(channel['sbtn'])
        self.setLayout(layout)


class GroupBox(QtWidgets.QGroupBox):

    def __init__(self, title):
        super().__init__(title)

        self.setStyleSheet('QGroupBox { font-weight: bold; } ')
        self.visible = True

    def set_state(self, state):
        self.visible = state
        if self.visible:
            for widget in self.children():
                if widget is not None and widget.isWidgetType():
                    widget.show()
        else:
            for widget in self.children():
                if widget is not None and widget.isWidgetType():
                    widget.hide()

    def set_visible(self):
        self.set_state(True)

    def set_hidden(self):
        self.set_state(False)

    def toggle(self):
        self.set_state(not self.visible)

    def mouseDoubleClickEvent(self, *args):
        self.toggle()


