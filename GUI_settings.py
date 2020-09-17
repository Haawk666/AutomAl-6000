# By Haakon Tvedt @ NTNU
"""Module container for style settings of the GUI"""

# External imports:
from PyQt5 import QtGui, QtCore

theme = 'dark'
tooltips = True

dark_palette = QtGui.QPalette()
dark_palette.setColor(QtGui.QPalette.Window, QtGui.QColor(53, 53, 53))
dark_palette.setColor(QtGui.QPalette.WindowText, QtGui.QColor(200, 200, 200))
dark_palette.setColor(QtGui.QPalette.Base, QtGui.QColor(25, 25, 25))
dark_palette.setColor(QtGui.QPalette.AlternateBase, QtGui.QColor(53, 53, 53))
dark_palette.setColor(QtGui.QPalette.ToolTipBase, QtCore.Qt.white)
dark_palette.setColor(QtGui.QPalette.ToolTipText, QtCore.Qt.white)
dark_palette.setColor(QtGui.QPalette.Text, QtGui.QColor(200, 200, 200))
dark_palette.setColor(QtGui.QPalette.Button, QtGui.QColor(53, 53, 53))
dark_palette.setColor(QtGui.QPalette.ButtonText, QtGui.QColor(200, 200, 200))
dark_palette.setColor(QtGui.QPalette.BrightText, QtCore.Qt.red)
dark_palette.setColor(QtGui.QPalette.Link, QtGui.QColor(42, 130, 218))
dark_palette.setColor(QtGui.QPalette.Highlight, QtGui.QColor(42, 130, 218))
dark_palette.setColor(QtGui.QPalette.HighlightedText, QtCore.Qt.black)

font_tiny = QtGui.QFont()
font_tiny.setPixelSize(9)

font_scalebar = QtGui.QFont('Helvetica [Cronyx]', 22)
font_mesh_details = QtGui.QFont('Helvetica [Cronyx]', 12)

pen_scalebar = QtGui.QPen(QtCore.Qt.white)
pen_scalebar.setWidth(12)

pen_skinny_red = QtGui.QPen(QtCore.Qt.red)
pen_skinny_red.setWidth(2)

white_font_tiny = QtGui.QFont()
white_font_tiny.setPixelSize(9)

white_pen = QtGui.QPen(QtCore.Qt.white)
white_pen.setWidth(1)

pen_boarder = QtGui.QPen(QtCore.Qt.black)
pen_boarder.setWidth(1)

brush_black = QtGui.QBrush(QtCore.Qt.black)
brush_white = QtGui.QBrush(QtCore.Qt.white)
background_brush = QtGui.QBrush(QtGui.QColor(180, 180, 180))

pen_al = QtGui.QPen(QtCore.Qt.green)
pen_al.setWidth(5)
brush_al = QtGui.QBrush(QtCore.Qt.green)

pen_mg = QtGui.QPen(QtGui.QColor(143, 0, 255))
pen_mg.setWidth(5)
brush_mg = QtGui.QBrush(QtGui.QColor(143, 0, 255))

pen_si = QtGui.QPen(QtCore.Qt.red)
pen_si.setWidth(5)
brush_si = QtGui.QBrush(QtCore.Qt.red)

pen_cu = QtGui.QPen(QtCore.Qt.yellow)
pen_cu.setWidth(5)
brush_cu = QtGui.QBrush(QtCore.Qt.yellow)

pen_zn = QtGui.QPen(QtGui.QColor(100, 100, 100))
pen_zn.setWidth(5)
brush_zn = QtGui.QBrush(QtGui.QColor(100, 100, 100))

pen_ag = QtGui.QPen(QtGui.QColor(200, 200, 200))
pen_ag.setWidth(5)
brush_ag = QtGui.QBrush(QtGui.QColor(200, 200, 200))

pen_un = QtGui.QPen(QtCore.Qt.blue)
pen_un.setWidth(5)
brush_un = QtGui.QBrush(QtCore.Qt.blue)

pen_selected_1 = QtGui.QPen(QtCore.Qt.darkCyan)
pen_selected_1.setWidth(6)
brush_selected_1 = QtGui.QBrush(QtCore.Qt.darkCyan)

pen_selected_2 = QtGui.QPen(QtCore.Qt.yellow)
pen_selected_2.setWidth(3)
brush_selected_2 = QtGui.QBrush(QtCore.Qt.transparent)

pen_atom_pos = QtGui.QPen(QtCore.Qt.red)
pen_atom_pos.setWidth(1)
brush_atom_pos = QtGui.QBrush(QtCore.Qt.transparent)

pen_atom_pos_hidden = QtGui.QPen(QtCore.Qt.red)
pen_atom_pos_hidden.setWidth(1)
brush_atom_pos_hidden = QtGui.QBrush(QtCore.Qt.transparent)

pen_graph = QtGui.QPen(QtCore.Qt.black)
pen_graph.setWidth(1)
brush_graph_0 = QtGui.QBrush(QtCore.Qt.black)
brush_graph_1 = QtGui.QBrush(QtCore.Qt.white)

pen_edge = QtGui.QPen(QtCore.Qt.black)
pen_edge.setWidth(1)
pen_inconsistent_edge = QtGui.QPen(QtCore.Qt.red)
pen_inconsistent_edge.setWidth(3)
pen_dislocation_edge = QtGui.QPen(QtCore.Qt.blue)
pen_dislocation_edge.setWidth(3)

default_config_string = '[theme]\n' \
                        'theme: dark\n\n' \
                        '[tooltips]\n' \
                        'tooltips: True\n\n' \
                        '[colors]\n'

