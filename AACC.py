"""Entry module for the GUI. Reads 'config.ini' and launches the PyQt5 application loop."""

#
# AutomAl 6000.
# ----------------------------------------
# By Haakon Tvedt.
# ----------------------------------------
# Master project in technical physics at NTNU. Supervised by Prof. Randi Holmestad. Co-supervised by Calin Maroiara
#

# Program imports:
import GUI
# External imports:
from PyQt5 import QtWidgets, QtGui, QtCore
import sys
import os
import configparser


def dump_log(program):
    string = program.terminal_window.handler.widget.toPlainText()
    with open('temp/log.txt', 'w') as f_:
        for line in iter(string.splitlines()):
            f_.write('{}\n'.format(line))


if __name__ == '__main__':

    app = QtWidgets.QApplication(sys.argv)
    app.setStyle("Fusion")

    default_config_string = '[theme]\n' \
                            'theme: dark\n\n' \
                            '[tooltips]\n' \
                            'tooltips: True\n\n' \
                            '[colors]\n'

    # Check for existence of config file:
    if not os.path.isfile('config.ini'):
        with open('config.ini', 'w') as f:
            f.write(default_config_string)

    # Import configurations from config file
    config = configparser.ConfigParser()
    config.read('config.ini')

    # Set theme
    if config.get('theme', 'theme') == 'dark':
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
        app.setPalette(dark_palette)
    else:
        pass

    # Start app
    program_session = GUI.MainUI(settings_file=config)
    app.exec_()
    dump_log(program_session)


