# By Haakon Tvedt @ NTNU
# Contributors:
"""API for automal 6000"""

import core
import logging
# Instantiate logger
logger = logging.getLogger(__name__)
logger.name = 'AutomAl 6000'


class Project:

    def __init__(self):
        self.project = None

    def set_image(self, filename, extension):
        self.project = core.Project.import_from_file(filename, extension)

    def set_model(self, model):
        if self.project is not None:
            self.project.graph.active_model = model
        else:
            logger.info('Project is None!')

    def find_columns(self, threshold):
        if self.project is not None:
            self.project.threshold = threshold
            self.project.column_detection(search_type='t', plot=False)
        else:
            logger.warning('No project set..')

    def column_characterization(self):
        if self.project.num_columns > 0:
            self.project.column_characterization(search_type=0, starting_index=0)
        else:
            logger.info('No columns detected!')

    def save(self, filename):
        if self.project is not None:
            self.project.save(filename)
        else:
            logger.info('No data to save!')

    @ staticmethod
    def load(filename):
        return core.Project.load(filename)



