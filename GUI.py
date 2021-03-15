# By Haakon Tvedt @ NTNU
# Contributors:
"""Module containing the Main Window class. This is the top-level GUI and contains the *business logic* of the
GUI."""

# Program imports:
import core
import GUI_elements
import utils
import data_module
# External imports:
from PyQt5 import QtWidgets, QtGui, QtCore
import numpy as np
import logging
# Instantiate logger
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
logger.name = 'GUI'


class MainUI(QtWidgets.QMainWindow):
    """Main GUI. Inherits PyQt5.QtWidgets.QMainWindow."""

    def __init__(self, *args, settings_file=None):
        super().__init__(*args)

        self.config = settings_file

        # Initialize in an 'empty state'
        self.project_instance = None
        self.savefile = None
        self.control_instance = None
        self.selected_column = -1
        self.previous_selected_column = -1
        self.selection_history = []
        self.lock_views = False
        self.overlay_settings = {
            'display_mode': 'atomic_species',  # 'atomic_species' or 'advanced_species'
            'overlay_radii': 0.5,  # Overlay radii * atomic radii
            'scalebar_length': 2,
            'scalebar_unit': 'nm',  # 'pm', 'nm', '$\mu$m' or 'Å'
            'scalebar_color': (255, 255, 255),
            'legend_background': 'Black'  # 'Black', 'White', 'Gray' or 'Transparent'
        }

        # Create menu bar
        self.menu = GUI_elements.MenuBar(self.menuBar(), self)

        # Create Control window
        self.control_window = None
        self.control_window_dock = None
        self.redraw_control_window()

        # Create terminal window
        self.terminal_window = GUI_elements.Terminal(obj=self)
        self.terminal_window.handler.set_mode(False)
        logger.addHandler(self.terminal_window.handler)

        self.terminal_window_scroll = QtWidgets.QScrollArea()
        self.terminal_window_scroll.setWidget(self.terminal_window)
        self.terminal_window_scroll.setWidgetResizable(True)
        self.terminal_window_scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)

        self.terminal_window_dock = QtWidgets.QDockWidget()
        self.terminal_window_dock.setWidget(self.terminal_window_scroll)
        self.terminal_window_dock.setWindowTitle('Terminal window')
        self.terminal_window_dock.setMinimumWidth(300)

        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, self.terminal_window_dock)

        # Tab contents:
        self.no_graphic = QtGui.QPixmap('Images\\no_image.png')
        self.graphic = QtGui.QPixmap('Images\\no_image.png')

        # gs = QGraphicsScene
        self.gs_raw_image = GUI_elements.RawImage(ui_obj=self, background=self.no_graphic)
        self.gs_atomic_positions = GUI_elements.AtomicPositions(ui_obj=self, background=self.no_graphic)
        self.gs_overlay_composition = GUI_elements.OverlayComposition(ui_obj=self, background=self.no_graphic)
        self.gs_atomic_graph = GUI_elements.AtomicGraph(ui_obj=self, background=self.no_graphic)
        self.gs_zeta_graph = GUI_elements.AtomicGraph(ui_obj=self, background=self.no_graphic, mode='zeta')
        self.gs_anti_graph = GUI_elements.AntiGraph(ui_obj=self, background=self.no_graphic)
        self.gs_heat = GUI_elements.HeatMap(ui_obj=self, background=self.no_graphic)
        self.gs_atomic_sub_graph = GUI_elements.AtomicSubGraph(ui_obj=self, background=self.no_graphic)
        self.gs_search_matrix = GUI_elements.RawImage(ui_obj=self, background=self.no_graphic)
        self.gs_fft = GUI_elements.RawImage(ui_obj=self, background=self.no_graphic)

        # gv = QGraphicsView
        self.gv_raw_image = GUI_elements.ZoomGraphicsView(self.gs_raw_image, ui_obj=self, trigger_func=self.key_press_trigger, tab_index=0, scale=1)
        self.gv_atomic_positions = GUI_elements.ZoomGraphicsView(self.gs_atomic_positions, ui_obj=self, trigger_func=self.key_press_trigger, tab_index=1, scale=1)
        self.gv_overlay_composition = GUI_elements.ZoomGraphicsView(self.gs_overlay_composition, ui_obj=self, trigger_func=self.key_press_trigger, tab_index=2, scale=1)
        self.gv_atomic_graph = GUI_elements.ZoomGraphicsView(self.gs_atomic_graph, ui_obj=self, trigger_func=self.key_press_trigger, tab_index=3, scale=2)
        self.gv_zeta_graph = GUI_elements.ZoomGraphicsView(self.gs_zeta_graph, ui_obj=self, trigger_func=self.key_press_trigger, tab_index=4, scale=2)
        self.gv_anti_graph = GUI_elements.ZoomGraphicsView(self.gs_anti_graph, ui_obj=self, trigger_func=self.key_press_trigger, tab_index=5, scale=2)
        self.gv_heat = GUI_elements.ZoomGraphicsView(self.gs_heat, ui_obj=self, trigger_func=self.key_press_trigger, tab_index=6, scale=1)
        self.gv_atomic_sub_graph = GUI_elements.ZoomGraphicsView(self.gs_atomic_sub_graph, ui_obj=self, trigger_func=self.key_press_trigger, tab_index=7, scale=2)
        self.gv_search_matrix = GUI_elements.ZoomGraphicsView(self.gs_search_matrix, ui_obj=self, trigger_func=self.key_press_trigger, tab_index=8, scale=1)
        self.gv_fft = GUI_elements.ZoomGraphicsView(self.gs_fft, ui_obj=self, trigger_func=self.key_press_trigger, tab_index=9, scale=1)

        self.gv_list = [
            self.gv_raw_image,
            self.gv_atomic_positions,
            self.gv_overlay_composition,
            self.gv_atomic_graph,
            self.gv_zeta_graph,
            self.gv_anti_graph,
            self.gv_heat,
            self.gv_atomic_sub_graph,
            self.gv_search_matrix,
            self.gv_fft
        ]

        # Set up tabs for central widget
        self.tabs = QtWidgets.QTabWidget()

        self.tab_raw_image = self.tabs.addTab(self.gv_raw_image, 'Raw image')
        self.tab_atomic_positions = self.tabs.addTab(self.gv_atomic_positions, 'Atomic positions')
        self.tab_overlay_composition = self.tabs.addTab(self.gv_overlay_composition, 'Overlay composition')
        self.tab_atomic_graph = self.tabs.addTab(self.gv_atomic_graph, 'Atomic graph')
        self.tab_zeta_graph = self.tabs.addTab(self.gv_zeta_graph, 'Zeta graph')
        self.tab_anti_graph = self.tabs.addTab(self.gv_anti_graph, 'Anti-graph')
        self.tab_heat = self.tabs.addTab(self.gv_heat, 'Heat-map')
        self.tab_atomic_sub_graph = self.tabs.addTab(self.gv_atomic_sub_graph, 'Atomic sub-graph')
        self.tab_search_matrix = self.tabs.addTab(self.gv_search_matrix, 'Search matrix')
        self.tab_fft = self.tabs.addTab(self.gv_fft, 'FFT image')

        self.setCentralWidget(self.tabs)

        # Generate elements
        self.setWindowTitle('AutomAl 6000')
        self.resize(1500, 900)
        self.move(50, 30)
        self.statusBar().showMessage('Ready')

        self.debug_mode = False

        # Display
        self.show()

        # Intro
        logger.info('Welcome to AutomAl 6000 by Haakon Tvedt')
        logger.info('version: {}.{}.{}'.format(core.Project.version[0], core.Project.version[1], core.Project.version[2]))
        logger.info('Build: Alpha 2.0')
        logger.info('See http://automal.org for help\n------------------------')

    # ----------
    # Business methods:
    # ----------

    def set_advanced_species(self, h):
        """Set advanced species of selected column"""
        if self.project_instance is not None and not self.selected_column == -1:
            # Update relevant graphics:
            self.project_instance.graph.set_advanced_species(self.selected_column, h)
            self.gs_overlay_composition.re_draw_vertex(self.selected_column)
            self.gs_atomic_graph.redraw_neighbourhood(self.selected_column)
            self.gs_zeta_graph.redraw_neighbourhood(self.selected_column)
            # Update control window info:
            self.control_window.select_column()

    def set_atomic_species(self, h):
        """Set advanced species of selected column"""
        if self.project_instance is not None and not self.selected_column == -1:
            # Update relevant graphics:
            self.project_instance.graph.set_atomic_species(self.selected_column, h)
            self.gs_overlay_composition.re_draw_vertex(self.selected_column)
            self.gs_atomic_graph.redraw_neighbourhood(self.selected_column)
            self.gs_zeta_graph.redraw_neighbourhood(self.selected_column)
            # Update control window info:
            self.control_window.select_column()

    def set_zeta(self, zeta):
        """Set level of selected column"""
        if self.project_instance is not None and not self.selected_column == -1:
            # Update relevant graphics:
            self.project_instance.graph.set_zeta(self.selected_column, zeta)
            self.gs_overlay_composition.re_draw_vertex(self.selected_column)
            self.gs_atomic_graph.interactive_vertex_objects[self.selected_column].set_style()
            self.gs_atomic_graph.redraw_neighbourhood(self.selected_column)
            self.gs_zeta_graph.interactive_vertex_objects[self.selected_column].set_style()
            self.gs_zeta_graph.redraw_neighbourhood(self.selected_column)
            # Update control window info:
            self.control_window.widgets['selected_column']['sbtn_zeta']['widget'].set_value(zeta)

    # ----------
    # Self state methods:
    # ----------

    def update_display(self):
        self.update_central_widget()
        self.redraw_control_window()
        self.sys_message('Ready.')

    def update_central_widget(self):
        if self.project_instance is not None:
            utils.im_out_static(self.project_instance.im_mat.astype(np.float64), 'Images\Outputs\Buffers\\raw_image.png')
            utils.im_out_static(self.project_instance.search_mat.astype(np.float64), 'Images\Outputs\Buffers\search_image.png')
            utils.im_out_static(self.project_instance.fft_im_mat.astype(np.float64), 'Images\Outputs\Buffers\FFT.png')
            void = False
        else:
            void = True
        self.update_raw_image(void=void)
        self.update_column_positions(void=void)
        self.update_overlay(void=void)
        self.update_graph()
        self.update_zeta_graph()

        self.update_search_matrix(void=void)
        self.update_fft(void=void)

    def update_raw_image(self, void=False):
        if void:
            graphic_ = self.no_graphic
        else:
            graphic_ = QtGui.QPixmap('Images\Outputs\Buffers\\raw_image.png')
        self.gs_raw_image = GUI_elements.RawImage(ui_obj=self, background=graphic_)
        self.gv_raw_image.setScene(self.gs_raw_image)

    def update_column_positions(self, void=False):
        if void:
            graphic_ = self.no_graphic
        else:
            graphic_ = QtGui.QPixmap('Images\Outputs\Buffers\\raw_image.png')
        self.gs_atomic_positions = GUI_elements.AtomicPositions(ui_obj=self, background=graphic_)
        self.gs_atomic_positions.re_draw()
        self.gv_atomic_positions.setScene(self.gs_atomic_positions)

    def update_overlay(self, void=False):
        if void:
            graphic_ = self.no_graphic
        else:
            graphic_ = QtGui.QPixmap('Images\Outputs\Buffers\\raw_image.png')
        self.gs_overlay_composition = GUI_elements.OverlayComposition(ui_obj=self, background=graphic_)
        self.gs_overlay_composition.re_draw()
        self.gv_overlay_composition.setScene(self.gs_overlay_composition)

    def update_graph(self):
        self.gs_atomic_graph = GUI_elements.AtomicGraph(ui_obj=self, scale_factor=self.gv_atomic_graph.scaling_factor)
        self.gv_atomic_graph.setScene(self.gs_atomic_graph)

    def update_zeta_graph(self):
        self.gs_zeta_graph = GUI_elements.AtomicGraph(ui_obj=self, scale_factor=self.gv_zeta_graph.scaling_factor, mode='zeta')
        self.gv_zeta_graph.setScene(self.gs_zeta_graph)

    def update_sub_graph(self):
        self.gs_atomic_sub_graph = GUI_elements.AtomicSubGraph(ui_obj=self, scale_factor=self.gv_atomic_sub_graph.scaling_factor)
        self.gv_atomic_sub_graph.setScene(self.gs_atomic_sub_graph)

    def update_anti_graph(self):
        if self.project_instance is not None:
            graph = self.project_instance.graph
        else:
            graph = None
        self.gs_anti_graph = GUI_elements.AntiGraph(
            ui_obj=self,
            scale_factor=self.gv_anti_graph.scaling_factor,
            graph=graph)
        self.gv_anti_graph.setScene(self.gs_anti_graph)

    def update_search_matrix(self, void=False):
        if void:
            graphic_ = self.no_graphic
        else:
            graphic_ = QtGui.QPixmap('Images\Outputs\Buffers\search_image.png')
        scene = GUI_elements.StaticImage(ui_obj=self, background=graphic_)
        self.gv_search_matrix.setScene(scene)

    def update_fft(self, void=False):
        if void:
            graphic_ = self.no_graphic
        else:
            graphic_ = QtGui.QPixmap('Images\Outputs\Buffers\FFT.png')
        scene = GUI_elements.StaticImage(ui_obj=self, background=graphic_)
        self.gv_fft.setScene(scene)

    def update_heat(self, kernel_size=1, step_size=1, attribute='normalized_peak_gamma', measure_type='variance', kernel_type='square'):
        if self.project_instance is not None:
            self.gs_heat = GUI_elements.HeatMap(
                ui_obj=self,
                background=QtGui.QPixmap('Images\Outputs\Buffers\\raw_image.png'),
                kernel_size=kernel_size,
                step_size=step_size,
                attribute=attribute,
                measure_type=measure_type,
                kernel_type=kernel_type
            )
            self.gs_heat.make()
            self.gv_heat.setScene(self.gs_heat)

    def update_control_window_test(self):
        self.control_window.update_display()

    def redraw_control_window(self):
        if self.control_window_dock is None:
            pos = 0
        else:
            pos = self.control_window_dock.widget().verticalScrollBar().value()
        self.removeDockWidget(self.control_window_dock)
        # Create Control window
        new_control = GUI_elements.ControlWindow(obj=self)
        if self.control_window is not None:
            for group in new_control.groups:
                new_control.groups[group].set_state(self.control_window.groups[group].visible)
                for widget in new_control.widgets[group]:
                    if widget in self.control_window.widgets[group]:
                        if new_control.widgets[group][widget]['widget_type'] == 'chb':
                            new_control.widgets[group][widget]['widget'].set_state_no_trigger(self.control_window.widgets[group][widget]['widget'].isChecked())

        self.control_window = new_control

        control_window_scroll = QtWidgets.QScrollArea()
        control_window_scroll.setWidget(self.control_window)
        control_window_scroll.setWidgetResizable(True)
        control_window_scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        control_window_scroll.verticalScrollBar().setValue(pos)

        self.control_window_dock = QtWidgets.QDockWidget()
        self.control_window_dock.setWidget(control_window_scroll)
        self.control_window_dock.setWindowTitle('Control window')
        self.control_window_dock.setMinimumWidth(300)

        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, self.control_window_dock)
        self.control_window.update_display()

    def column_selected(self, i):
        if self.control_window.widgets['selected_column']['mchb_move']['widget'].isChecked():
            pass
        else:
            self.previous_selected_column, self.selected_column = self.selected_column, i
            self.control_window.select_column()
            j = self.previous_selected_column
            if not j == -1:
                self.gs_atomic_positions.interactive_position_objects[j].set_style()
                self.gs_overlay_composition.interactive_overlay_objects[j].set_style()
                self.gs_atomic_graph.interactive_vertex_objects[j].set_style()
                self.gs_zeta_graph.interactive_vertex_objects[j].set_style()
            if not i == -1:
                self.gs_atomic_positions.interactive_position_objects[i].set_style()
                self.gs_overlay_composition.interactive_overlay_objects[i].set_style()
                self.gs_atomic_graph.interactive_vertex_objects[i].set_style()
                self.gs_zeta_graph.interactive_vertex_objects[i].set_style()
            if self.control_window.widgets['atomic_graph']['chb_permute_mode']['widget'].isChecked():
                if len(self.selection_history) == 2:
                    if self.tabs.currentIndex() == 3:
                        self.gs_atomic_graph.perturb_edge(self.selection_history[0], self.selection_history[1], self.selected_column)
                    elif self.tabs.currentIndex() == 4:
                        self.gs_zeta_graph.perturb_edge(self.selection_history[0], self.selection_history[1], self.selected_column)
                    self.selection_history = []
                else:
                    self.selection_history.append(self.selected_column)
            if self.control_window.widgets['atomic_graph']['chb_enable_ruler']['widget'].isChecked():
                projected_distance = self.project_instance.graph.get_projected_separation(self.selected_column, self.previous_selected_column)
                spatial_distance = self.project_instance.graph.get_separation(self.selected_column, self.previous_selected_column)
                expected_hard_sphere_distance = self.project_instance.graph.get_hard_sphere_separation(self.selected_column, self.previous_selected_column)
                string = 'Distance between vertex {} and {}\n' \
                         '    Projected distance: {} pm\n' \
                         '    Spatial distance: {} pm\n' \
                         '    Expected hard-sphere distance: {} pm\n' \
                         '    Deviation from hard sphere: {} pm\n'.format(self.previous_selected_column,
                                                                          self.selected_column, projected_distance,
                                                                          spatial_distance,
                                                                          expected_hard_sphere_distance,
                                                                          expected_hard_sphere_distance - spatial_distance)
                logger.info(string)

    def sys_message(self, msg):
        self.statusBar().showMessage(msg)
        QtWidgets.QApplication.processEvents()

    # ----------
    # Keyboard press methods:
    # ----------

    def keyPressEvent(self, event):
        """Handles key-presses when central widget has focus. Used to switch between tabs"""
        self.key_press_trigger(event.key())

    def key_press_trigger(self, key):
        """Process key-press events from graphic elements"""
        if self.project_instance is not None:
            if key == QtCore.Qt.Key_1:
                if not self.selected_column == -1:
                    self.project_instance.graph.vertices[self.selected_column].is_set_by_user = True
                    self.set_atomic_species('Si')
            elif key == QtCore.Qt.Key_2:
                if not self.selected_column == -1:
                    self.project_instance.graph.vertices[self.selected_column].is_set_by_user = True
                    self.set_atomic_species('Cu')
            elif key == QtCore.Qt.Key_3:
                if not self.selected_column == -1:
                    self.project_instance.graph.vertices[self.selected_column].is_set_by_user = True
                    self.set_atomic_species('Al')
            elif key == QtCore.Qt.Key_4:
                if not self.selected_column == -1:
                    self.project_instance.graph.vertices[self.selected_column].is_set_by_user = True
                    self.set_atomic_species('Mg')
            elif key == QtCore.Qt.Key_Plus:
                if not self.selected_column == -1:
                    self.set_zeta(self.project_instance.graph.vertices[self.selected_column].anti_zeta())
            elif key == QtCore.Qt.Key_W:
                if self.control_window.widgets['selected_column']['mchb_move']['widget'].isChecked():
                    self.gs_atomic_positions.interactive_position_objects[self.selected_column].moveBy(0.0, -1.0)
                else:
                    self.gv_list[self.tabs.currentIndex()].translate_w(10)
            elif key == QtCore.Qt.Key_S:
                if self.control_window.widgets['selected_column']['mchb_move']['widget'].isChecked():
                    self.gs_atomic_positions.interactive_position_objects[self.selected_column].moveBy(0.0, 1.0)
                else:
                    self.gv_list[self.tabs.currentIndex()].translate_s(10)
            elif key == QtCore.Qt.Key_A:
                if self.control_window.widgets['selected_column']['mchb_move']['widget'].isChecked():
                    self.gs_atomic_positions.interactive_position_objects[self.selected_column].moveBy(-1.0, 0.0)
                else:
                    self.gv_list[self.tabs.currentIndex()].translate_a(10)
            elif key == QtCore.Qt.Key_D:
                if self.control_window.widgets['selected_column']['mchb_move']['widget'].isChecked():
                    self.gs_atomic_positions.interactive_position_objects[self.selected_column].moveBy(1.0, 0.0)
                else:
                    self.gv_list[self.tabs.currentIndex()].translate_d(10)
            elif key == QtCore.Qt.Key_E:
                self.control_window.btn_align_views_trigger()
            elif key == QtCore.Qt.Key_Q:
                self.control_window.btn_snap_trigger()
            elif key == QtCore.Qt.Key_Z:
                if self.tabs.currentIndex() == 0:
                    self.tabs.setCurrentIndex(len(self.tabs) - 1)
                else:
                    self.tabs.setCurrentIndex(self.tabs.currentIndex() - 1)
            elif key == QtCore.Qt.Key_X:
                if self.tabs.currentIndex() == len(self.tabs) - 1:
                    self.tabs.setCurrentIndex(0)
                else:
                    self.tabs.setCurrentIndex(self.tabs.currentIndex() + 1)
            elif key == QtCore.Qt.Key_F1:
                if not self.selected_column == -1:
                    self.project_instance.graph.vertices[self.selected_column].flag_1 = not self.project_instance.graph.vertices[self.selected_column].flag_1
                    logger.info('vertex {}, flag 1 set to {}'.format(self.selected_column, self.project_instance.graph.vertices[self.selected_column].flag_1))
            elif key == QtCore.Qt.Key_F2:
                if not self.selected_column == -1:
                    self.project_instance.graph.vertices[self.selected_column].flag_2 = not self.project_instance.graph.vertices[self.selected_column].flag_2
                    logger.info('vertex {}, flag 2 set to {}'.format(self.selected_column, self.project_instance.graph.vertices[self.selected_column].flag_2))
            elif key == QtCore.Qt.Key_F3:
                if not self.selected_column == -1:
                    self.project_instance.graph.vertices[self.selected_column].flag_3 = not self.project_instance.graph.vertices[self.selected_column].flag_3
                    logger.info('vertex {}, flag 3 set to {}'.format(self.selected_column, self.project_instance.graph.vertices[self.selected_column].flag_3))
            elif key == QtCore.Qt.Key_F4:
                if not self.selected_column == -1:
                    self.project_instance.graph.vertices[self.selected_column].flag_4 = not self.project_instance.graph.vertices[self.selected_column].flag_4
                    logger.info('vertex {}, flag 4 set to {}'.format(self.selected_column, self.project_instance.graph.vertices[self.selected_column].flag_4))
            elif key == QtCore.Qt.Key_Space:
                if self.tabs.currentIndex() == 3:
                    self.control_window.widgets['atomic_graph']['chb_permute_mode']['widget'].toggle()
                    self.selection_history = []
                else:
                    pass
            elif key == QtCore.Qt.Key_R:
                if not self.selected_column == -1:
                    self.btn_show_stats_trigger()
            elif key == QtCore.Qt.Key_P:
                if self.control_window.widgets['selected_column']['mchb_move']['widget'].isChecked():
                    self.control_window.mchb_move_trigger(False)
                else:
                    self.control_window.mchb_move_trigger(True)
            elif key == QtCore.Qt.Key_Enter:
                if self.control_window.widgets['selected_column']['mchb_move']['widget'].isChecked():
                    self.control_window.btn_set_position_trigger()

    # ----------
    # PNG export:
    # ----------

    def export_raw_image_trigger(self):
        if self.project_instance is not None:
            filename = QtWidgets.QFileDialog.getSaveFileName(self, "Save image", '',
                                                             "PNG (*.png);;BMP Files (*.bmp);;JPEG (*.JPEG)")
            if filename[0]:
                self.update_raw_image()
                rect_f = self.gs_raw_image.sceneRect()
                img = QtGui.QImage(rect_f.size().toSize(), QtGui.QImage.Format_ARGB32)
                img.fill(QtCore.Qt.white)
                p = QtGui.QPainter(img)
                self.gs_raw_image.render(p, target=QtCore.QRectF(img.rect()), source=rect_f)
                p.end()
                saved = img.save(filename[0])
                if saved:
                    logger.info('Successfully exported raw image to file!')
                else:
                    logger.error('Could not export image!')

    def export_column_position_image_trigger(self):
        if self.project_instance is not None:
            filename = QtWidgets.QFileDialog.getSaveFileName(self, "Save image", '',
                                                             "PNG (*.png);;BMP Files (*.bmp);;JPEG (*.JPEG)")
            if filename[0]:
                self.update_column_positions()
                rect_f = self.gs_atomic_positions.sceneRect()
                img = QtGui.QImage(rect_f.size().toSize(), QtGui.QImage.Format_ARGB32)
                img.fill(QtCore.Qt.white)
                p = QtGui.QPainter(img)
                self.gs_atomic_positions.render(p, target=QtCore.QRectF(img.rect()), source=rect_f)
                p.end()
                saved = img.save(filename[0])
                if saved:
                    logger.info('Successfully exported column positions image to file!')
                else:
                    logger.error('Could not export image!')

    def export_overlay_image_trigger(self):
        if self.project_instance is not None:
            filename = QtWidgets.QFileDialog.getSaveFileName(self, "Save image", '',
                                                             "PNG (*.png);;BMP Files (*.bmp);;JPEG (*.JPEG)")
            if filename[0]:
                rect_f = self.gs_overlay_composition.sceneRect()
                img = QtGui.QImage(rect_f.size().toSize(), QtGui.QImage.Format_ARGB32)
                img.fill(QtCore.Qt.white)
                p = QtGui.QPainter(img)
                self.gs_overlay_composition.render(p, target=QtCore.QRectF(img.rect()), source=rect_f)
                p.end()
                saved = img.save(filename[0])
                if saved:
                    logger.info('Successfully exported overlay image to file!')
                else:
                    logger.error('Could not export image!')

    def export_atomic_graph_trigger(self):
        if self.project_instance is not None:
            filename = QtWidgets.QFileDialog.getSaveFileName(self, "Save image", '',
                                                             "PNG (*.png);;BMP Files (*.bmp);;JPEG (*.JPEG)")
            if filename[0]:
                self.update_graph()
                rect_f = self.gs_atomic_graph.sceneRect()
                img = QtGui.QImage(rect_f.size().toSize(), QtGui.QImage.Format_ARGB32)
                img.fill(QtCore.Qt.white)
                p = QtGui.QPainter(img)
                self.gs_atomic_graph.render(p, target=QtCore.QRectF(img.rect()), source=rect_f)
                p.end()
                saved = img.save(filename[0])
                if saved:
                    logger.info('Successfully exported graph image to file!')
                else:
                    logger.error('Could not export image!')

    def export_zeta_graph_trigger(self):
        if self.project_instance is not None:
            filename = QtWidgets.QFileDialog.getSaveFileName(self, "Save image", '',
                                                             "PNG (*.png);;BMP Files (*.bmp);;JPEG (*.JPEG)")
            if filename[0]:
                self.update_graph()
                rect_f = self.gs_zeta_graph.sceneRect()
                img = QtGui.QImage(rect_f.size().toSize(), QtGui.QImage.Format_ARGB32)
                img.fill(QtCore.Qt.white)
                p = QtGui.QPainter(img)
                self.gs_zeta_graph.render(p, target=QtCore.QRectF(img.rect()), source=rect_f)
                p.end()
                saved = img.save(filename[0])
                if saved:
                    logger.info('Successfully exported graph image to file!')
                else:
                    logger.error('Could not export image!')

    def export_heat_trigger(self):
        if self.project_instance is not None:
            filename = QtWidgets.QFileDialog.getSaveFileName(self, "Save image", '',
                                                             "PNG (*.png);;BMP Files (*.bmp);;JPEG (*.JPEG)")
            if filename[0]:
                rect_f = self.gs_heat.sceneRect()
                img = QtGui.QImage(rect_f.size().toSize(), QtGui.QImage.Format_ARGB32)
                img.fill(QtCore.Qt.white)
                p = QtGui.QPainter(img)
                self.gs_heat.render(p, target=QtCore.QRectF(img.rect()), source=rect_f)
                p.end()
                if img.save(filename[0]):
                    logger.info('Successfully exported heat map to file!')
                else:
                    logger.error('Could not export image!')

    # ----------
    # Menu triggers:
    # ----------

    def import_trigger(self, file_type, scale=1.0):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Select dm3', '')
        if filename[0]:
            self.sys_message('Working...')
            self.project_instance = core.Project.import_from_file(filename[0], file_type, scale=scale)
            self.savefile = None
            self.control_instance = None
            self.update_display()
            self.control_window.btn_species_dict_trigger()
        else:
            self.sys_message('Ready')

    def menu_new_trigger(self):
        GUI_elements.ImportWizard(ui_obj=self)

    def menu_open_trigger(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', '')
        if filename[0]:
            logger.info('Opening file {}'.format(filename[0]))
            self.sys_message('Working...')
            self.project_instance = core.Project.load(filename[0])
            if self.project_instance is not None:
                self.project_instance.debug_mode = False
                self.terminal_window.handler.set_mode(False)
                self.control_instance = None
                self.savefile = filename[0]
                self.update_display()
                self.sys_message('Ready')
            else:
                logger.info('File was not loaded. Something must have gone wrong!')
                self.sys_message('Ready')
        else:
            self.sys_message('Ready')

    def menu_save_trigger(self):
        if self.savefile is None:
            filename = QtWidgets.QFileDialog.getSaveFileName(self, 'Save file', '')
        else:
            filename = QtWidgets.QFileDialog.getSaveFileName(self, 'Save file', self.savefile)

        if filename[0]:
            self.sys_message('Working...')
            self.project_instance.save(filename[0])
            self.savefile = filename[0]
            self.update_display()
            self.sys_message('Ready.')
        else:
            self.sys_message('Ready.')

    def menu_close_trigger(self):
        self.sys_message('Working...')
        self.control_window.btn_cancel_move_trigger()
        self.column_selected(-1)
        self.project_instance = None
        self.control_instance = None
        self.update_display()
        self.sys_message('Ready.')

    def menu_exit_trigger(self):
        self.close()

    def menu_update_display(self):
        self.sys_message('Working...')
        self.update_display()
        self.sys_message('Ready.')

    def menu_toggle_debug_mode_trigger(self, state):
        if state:
            self.control_window.debug_box.set_visible()
        else:
            self.control_window.debug_box.set_hidden()

    def menu_add_mark_trigger(self):
        logger.info('Column positions: {}'.format(self.gv_atomic_positions.size()))
        logger.info('Graph: {}'.format(self.gv_atomic_positions.size()))
        logger.info('-------------')

    def menu_clear_flags_trigger(self):
        if self.project_instance is not None:
            self.project_instance.graph.reset_all_flags()

    def menu_set_control_file_trigger(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self, 'Open control file', '')
        if filename[0]:
            self.control_instance = core.Project.load(filename[0])
            if self.control_instance is not None:
                self.project_instance.debug_mode = False
            else:
                logger.info('Control-file was not loaded. Something must have gone wrong!')

    def menu_run_benchmark_trigger(self):
        pass

    def menu_display_deviations_trigger(self):
        if self.project_instance is not None and self.control_instance is not None and self.project_instance.num_columns > 0:
            num_errors = 0
            num_precipitate_errors = 0
            num_columns = 0
            num_precipitate_columns = 0

            for vertex in self.project_instance.graph.vertices:
                if not vertex.is_edge_column:
                    if not vertex.atomic_species == self.control_instance.graph.vertices[vertex.i].atomic_species:
                        num_errors += 1
                        if vertex.is_in_precipitate:
                            num_precipitate_errors += 1
                            num_precipitate_columns += 1
                    else:
                        if vertex.is_in_precipitate:
                            num_precipitate_columns += 1
                    num_columns += 1

            if not num_columns == 0:
                percent_error = 100 * (num_errors / num_columns)
            else:
                percent_error = 0

            if not num_precipitate_columns == 0:
                percent_precipitate_error = 100 * (num_precipitate_errors / num_precipitate_columns)
            else:
                percent_precipitate_error = 0

            string = 'Control comparison:\n'
            string += '    precipitate errors (#): {}\n'.format(num_precipitate_errors)
            string += '    precipitate error (%): {}\n\n'.format(percent_precipitate_error)
            string += '    total errors (#): {}\n'.format(num_errors)
            string += '    total error (%): {}\n'.format(percent_error)

            logger.info(string)

    def menu_test_consistency_trigger(self):
        pass

    def menu_ad_hoc_trigger(self):
        pass

    def menu_toggle_tooltips_trigger(self, state):
        self.control_window.mode_tooltip(state)
        self.terminal_window.mode_tooltip(state)
        self.config.set('tooltips', 'tooltips', str(state))
        with open('config.ini', 'w') as configfile:
            self.config.write(configfile)

    def menu_set_theme_trigger(self):
        items = ('dark', 'classic')
        item, ok_pressed = QtWidgets.QInputDialog.getItem(self, "Select theme", "Theme:", items, 0, False)
        if ok_pressed and item:
            self.config.set('theme', 'theme', item)
            with open('config.ini', 'w') as configfile:
                self.config.write(configfile)
            message = QtWidgets.QMessageBox()
            message.setText('Save your work and restart the program for the changes to take effect!')
            message.exec_()

    def menu_there_is_no_help_trigger(self):
        message = QtWidgets.QMessageBox()
        message.setText('Not implemented yet')
        message.exec_()

    # ----------
    # Other button triggers:
    # ----------

    def btn_save_log_trigger(self):
        filename = QtWidgets.QFileDialog.getSaveFileName(self, 'Save log-file', '')
        if filename[0]:
            self.sys_message('Working...')
            logger.info('Saving log file...')
            string = self.terminal_window.handler.widget.toPlainText()
            with open(filename[0], 'w') as f:
                for line in iter(string.splitlines()):
                    f.write('{}\n'.format(line))
            f.close()
            self.sys_message('Ready.')
            logger.info('Saved log to {}'.format(filename[0]))

    def btn_clear_log_trigger(self):
        self.terminal_window.handler.widget.clear()

    def btn_eval_trigger(self, string):
        try:
            eval(string)
        except:
            logger.info('The input\n    {}\n    threw an exception!'.format(string))
        else:
            logger.info('User input: {}'.format(string))

    # ----------
    # Checkbox triggers (DEPRECATED):
    # ----------

    def chb_placeholder_trigger(self):
        pass


