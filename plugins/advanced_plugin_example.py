import logging
import plugins.advanced_plugin_example_files.advanced_plugin_example_GUI as my_GUI

# Instantiate logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Bridge:

    def __init__(self, ui_obj):
        self.ui_obj = ui_obj
        self.my_gui = None

    def trigger(self):
        if self.my_gui is None:
            self.my_gui = my_GUI.MyGUI(ui_obj=self.ui_obj)
        self.my_gui.show()

