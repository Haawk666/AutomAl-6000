from PyQt5 import QtWidgets, QtCore
from .. import advanced_plugin_example as my_api


class MyGUI(QtWidgets.QDialog):

    def __init__(self, *args, ui_obj=None):
        super().__init__(*args)

        self.ui_obj = ui_obj

        self.btn_print_num_mg = QtWidgets.QPushButton('Print number of Mg in image')
        self.btn_print_num_mg.clicked.connect(self.print_mg_trigger)

        self.chb_toggle_all_mg_in_overlay = QtWidgets.QCheckBox('Toggle overlay Mg')
        self.chb_toggle_all_mg_in_overlay.toggled.connect(self.toggle_mg_trigger)

        self.btn_layout = QtWidgets.QHBoxLayout()
        self.btn_layout.addStretch()
        self.btn_layout.addWidget(self.btn_print_num_mg)
        self.btn_layout.addWidget(self.chb_toggle_all_mg_in_overlay)
        self.btn_layout.addStretch()

        self.setLayout(self.btn_layout)

        self.ui_obj.plugin_dock = QtWidgets.QDockWidget()
        self.ui_obj.plugin_dock.setWidget(self)
        self.ui_obj.plugin_dock.setWindowTitle('My advanced plugin')

        self.ui_obj.addDockWidget(QtCore.Qt.BottomDockWidgetArea, self.ui_obj.plugin_dock)

    def show(self):
        super().show()

        self.chb_toggle_all_mg_in_overlay.blockSignals(True)
        self.chb_toggle_all_mg_in_overlay.setChecked(self.ui_obj.control_window.chb_mg_columns.isChecked())
        self.chb_toggle_all_mg_in_overlay.blockSignals(False)

    def print_mg_trigger(self):

        if self.ui_obj.project_instance is not None:
            num_mg = 0
            for vertex in self.ui_obj.project_instance.graph.vertices:
                if vertex.species() == 'Mg':
                    num_mg += 1
            my_api.logger.info('Number of Mg: {}'.format(num_mg))
        else:
            my_api.logger.warning('There is no active project!')

    def toggle_mg_trigger(self, state):

        if self.ui_obj.project_instance is not None:
            self.ui_obj.control_window.chb_mg_columns.setChecked(state)
        else:
            my_api.logger.warning('There is no active project!')
            self.chb_toggle_all_mg_in_overlay.blockSignals(True)
            self.chb_toggle_all_mg_in_overlay.setChecked(self.ui_obj.control_window.chb_mg_columns.isChecked())
            self.chb_toggle_all_mg_in_overlay.blockSignals(False)

